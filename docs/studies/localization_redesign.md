# Compartment assignment — flux-free placement + materialised-FBA certification

The design of `assign_compartments` (in `src/raven_toolbox/localization/certify.py`): place reactions
into subcellular compartments by localisation score, then certify functionality with a real FBA on the
materialised model. It replaced a legacy approach that coupled a flux model into the placement
optimisation and could, at genome scale, certify a placement the materialised model cannot actually
grow; keeping the flux model out of the optimisation removes that failure mode at the source.

/ Certainty is flagged: **CERTAIN** = backed by data or by construction. **DESIGN** = an intended
property, validated where noted. /

---

## 1. The idea

**Take the flux model out of the placement problem.** Place reactions by localisation score in a
*flux-free* MILP (sound by construction — no flux variable exists to leak), then make functionality a
*separate, tractable* step whose certificate is a real `cobra` FBA on the model `apply_assignment`
actually builds. Because the certificate is the same FBA that is the authoritative acceptance test, the
optimiser can never certify a placement the materialised model cannot grow.

The concrete method is **Place-then-Repair with a functionality feedback loop**:

1. a flux-free **placement master** (score maximisation, like the existing `predict_localization`);
2. a **structural confinement repair** that co-locates reactions sharing a non-transportable metabolite
   and adds transports for the transportable pools a placement splits;
3. **materialised-FBA certification** over the primary medium and every `GrowthCondition`;
4. a **feedback loop** that, on a real growth failure, tightens placement (co-location / pin) or
   gap-fills — so functionality still drives placement, but through real flux.

This reuses primitives the repository already has: a sound score-only placement (`predict_localization`,
in `src/raven_toolbox/localization/predict.py`) and the repair/materialisation primitives
(`apply_assignment`, `manipulation.transport.add_transport_reactions`,
`cobra.flux_analysis.find_blocked_reactions` / `pfba` / `gapfill`). `assign_compartments` returns an
`AssignmentProposal` marked `certified=True` **iff** a real FBA on the materialised model reaches the
growth floor on every medium.

---

## 2. Architecture

### 2.1 Shared scope setup

Movable-vs-pinned split (boundary reactions, existing multi-compartment reactions and the biomass
reaction are always pinned); `genes_in_scope` and the per-gene per-compartment score lookup;
`base_metabolite` keying (default strips `_<compartment>`; yeast-GEM uses `lambda m: m.name`); the
`transportable` set of base metabolites; the biomass reaction and `min_growth` (default 10 % of the
draft FBA optimum). Nothing about scoring intent changes.

### 2.2 Placement master — flux-free (CERTAIN sound)

A MILP over only `x[r,c]` (mono-localisation `Σ_c x[r,c] = 1` by default), `y[g,c]` with the existing
`couple_`/`has_`/`gene1_` coupling, and the score objective:

```
max Σ score·y − Σ multi_compartment_penalty·y_extra
```

There is **no** flux variable, no mass-balance row, and no growth constraint — a clean assignment
problem with 0/±1 coefficients and an `O(1)` score objective. `predict_localization` already solves a
problem of exactly this shape at genome scale, direct evidence the master is tractable. A leak is
**unrepresentable**: there is nothing for a tolerance-rounded binary to leak into.

The master also accepts two feedback inputs (empty on the first round):

- **co-location groups** — sets of movable reactions constrained to share one compartment
  (`x[a,c] = x[b,c] ∀c`); the score objective then chooses *which* compartment, so the higher-weight
  gene anchors the group;
- **forced pins** — `x[r, c] = 1` when a reaction must align with a *pinned* reaction's compartment.

### 2.3 Structural confinement repair (CERTAIN)

Given a placement, a base metabolite is **stranded** if it is *non-transportable* yet touched by
placed/pinned reactions in more than one compartment (it cannot be reconnected by a transport). For
each such metabolite, tie its reactions together:

- all-movable → a **co-location group** (let the master's score pick the compartment);
- involves a pinned reaction (e.g. a biomass precursor consumed by the pinned biomass reaction) →
  a **forced pin** to that pinned compartment.

Re-solve the master with the new constraints; iterate until no confined metabolite is split. This is a
purely structural fixed point (no FBA, deterministic) and it is exactly the mechanism by which
"functionality drives placement against a reaction's top score": a confined precursor forces its
producer to where it is consumed.

For **transportable** pools a placement splits (present in ≥2 compartments), add star-topology
transports through the default compartment — the connectivity `apply_assignment` materialises.

### 2.4 Materialised-FBA certification (CERTAIN — the certificate is the acceptance test)

Build the `AssignmentProposal`, call `apply_assignment`, then for the primary medium and every
`GrowthCondition`: set `model.medium` and `slim_optimize`. This FBA has no placement binary — a
reaction is in the model or not — so its optimum is the *true* biomass flux. If all media clear their
floor → **certified**. Certification and acceptance are the same computation, so they can never
disagree.

### 2.5 Feedback on real failure (DESIGN — preserves "functionality is emergent")

If a medium fails even after §2.3, the failure is a genuine gap, not a confinement artifact:

- **gap-fill** (only if `universal` supplied): the minimum reactions to restore growth — strictly
  optional, parsimonious;
- **placement gap**: diagnose the unproducible precursor pool (`find_blocked_reactions` + a solver-free
  producibility check) and re-open the reactions feeding it in the master with a functionality
  incentive toward the compartment where the pool *is* producible.

Cap the rounds; return the best **certified** proposal, or — if the cap is hit — the best placement
*plus the explicit media it fails and the transports/reactions it still needs*. Never a false
certificate.

### 2.6 Multi-localisation and multi-medium

- **Multi-localisation** is expressed as certification, not as an in-optimisation activity-binary
  family (which does not scale): a reaction is multi-placed only if a loopless FVA on the *materialised*
  model shows both placements carry real flux ≥ `eps`. This keeps the master tractable while a dead
  duplicate — which would harvest a compartment's score for free — cannot pass the FVA.
- **`GrowthCondition`** media *are* the certification FBA (§2.4) and the source of feedback (§2.5) over
  the shared placement — so a transport or placement a respiration / β-oxidation medium needs is added
  because that real medium's FBA failed, not because of any in-optimisation proxy.

---

## 3. Why it scales (DESIGN)

The optimisation has **zero continuous variables and zero flux gates**:

- **Master size:** ~28k `x` + ~16k `y` binaries with only the cheap coupling rows; it carries none of
  the mass-balance equalities, placement/transport gates, or continuous fluxes a fused formulation
  would. Coefficients are 0/±1 or `O(1)` scores, so the LP relaxation is tight and Gurobi
  presolve/clique-cover handle the set-cover shape well.
- **Oracle cost per round:** one `apply_assignment` (deep copy + `O(reactions)` edits) plus 1–14
  `slim_optimize` calls (~0.1–0.5 s each on a ~2600-reaction model). Dominated by the master solve.
- **Rounds:** confinement repair is a deterministic structural fixed point; the common genome-scale
  case (everything transportable via evidence-aware cost) triggers **no** co-location and certifies in
  a single round after adding transports for the split pools. Feedback rounds fire only for genuine
  gaps, over a low-dimensional deficit space (~20 recurring biomass precursor pools).

---

## 4. How it preserves intent

- **Functionality drives placement.** A reaction lands against its top DeepLoc score exactly when a
  confined precursor (or a certified growth failure) forces it there — enforced by real flux.
  `test_functionality_overrides_score` and `test_pathway_coherence_overrides_individual_score` remain
  meaningful: a confined metabolite fires §2.3 and forces co-location / pinning.
- **Evidence-aware transports.** `evidence_aware_transport_cost` still supplies the `{base_met: cost}`
  mapping weighting which split pool is worth a transport and which reconnection is preferred.
- **Gap-fill stays strictly optional and parsimonious** — generated only when a real FBA proves
  reconnection alone cannot restore growth; with `universal=None` it never fires.
- **The carbon panel is the authoritative test**, and it is *the same object* the algorithm optimises
  against.

---

## 5. Soundness guarantee (CERTAIN, by construction)

The returned proposal is marked **certified** *if and only if* `apply_assignment` + FBA reaches the
growth floor on every medium — because that FBA *is* the certification. There is no separate optimistic
certificate that can disagree with the materialised model. A placement that cannot grow is returned
marked **uncertified**, with the failing media listed — never a silent false positive. This is the
single property the whole design exists to guarantee, and it holds by construction rather than by
solver tolerance.

---

## 6. Results (CERTAIN — measured on the yeast-GEM replication task)

`scripts/benchmark_certified_yeast.py` on the flattened yeast-GEM draft (2569 rxns / 1337 mets,
draft growth 0.1426; `min_growth = 0.5 × 0.0809 = 0.0405`; scalar `transport_cost = 0.5`):

| Method | wall | certified? | materialised growth | reaction agr | gene agr | transports |
|---|---|--:|--:|--:|--:|--|
| **`assign_compartments`** | **~23 s** | **yes** | **0.1426** (floor 0.0405) | **72.5 %** (1408/1943) | **88.7 %** (716/807) | 1261 full-connect / **983 FVA-usable** (665 base mets) |

Reading:

- **The certified path gets high agreement *and* real functionality, in ~23 s.** 72.5 % reaction /
  88.7 % gene agreement with a materialised model that genuinely grows at 0.1426. It **certified in a
  single round**: every base metabolite was transportable, so no confinement co-location fired; the
  placement was score-optimal and the repair added transports for the pools it split. 0 reactions were
  left unplaced.
- **The honest agreement-vs-transport frontier is now measurable.** 72.5 % agreement costs ~1261
  reconnecting transports (983 FVA-usable), against curated yeast-GEM's own 1467 inter-compartment
  reactions — high agreement honestly requires many transports (there is no free lunch).

Every number is authoritative because the growth column *is* the certificate: there is no separate MILP
floor to disagree with it. Raw JSON in `.research_tmp/certify_yeast_gem.json`.

_Validation: the unit tests in `tests/test_localization_certify.py` cover the
certified-implies-materialised-growth invariant, functionality-overrides-score, pathway coherence,
gap-fill, per-medium growth conditions, transport minimisation, and multi-localisation, plus two
adversarially-found repair bugs (non-default-hub transport drop; diagnoser oscillation)._
