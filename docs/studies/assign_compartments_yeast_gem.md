# Benchmark: reproducing yeast-GEM compartmentalisation from UniProt evidence

**Question.** Given only *agnostic* subcellular-localisation evidence (UniProt annotations), how
well does `assign_compartments` reproduce the **curated** compartmentalisation of a real fungal
genome-scale model — and is the result a model that still grows?

This is the end-to-end test of the whole pipeline:

```
UniProt "Subcellular location"  ──fetch_uniprot_localization──▶  gene × compartment scores
                                ──assign_compartments──────────▶  functional, compartmentalised model
```

## Setup

* **Model:** [yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM) (~4100 reactions, ~2750
  metabolites, 14 compartments, ~1140 genes). Each reaction's curated compartment is the ground
  truth.
* **Evidence:** `fetch_uniprot_localization(559292)` — UniProt's reviewed *S. cerevisiae* entries,
  the `Subcellular location [CC]` annotation mapped to model compartment ids via
  `DEFAULT_COMPARTMENT_MAP`, keyed by ordered-locus (ORF) id so it lines up with yeast-GEM genes.
  Each annotated compartment scores 1.0 (UniProt's curated location is qualitative).
* **Task:** relocate every internal, single-compartment reaction whose gene UniProt annotates,
  in a compartment the evidence can speak to; ask the MILP to re-place them from the evidence while
  keeping biomass producible. `base_metabolite = name` (yeast-GEM ids are compartment-specific),
  `transportable = []` (no new transporters — pure relocation, tractable at genome scale).

Reproduce with [`scripts/assign_compartments_yeast_uniprot.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/scripts/assign_compartments_yeast_uniprot.py).

## Result

Agnostic UniProt evidence, reconciled with functionality, reproduces **~74–83 %** of the curated
compartmentalisation (~4450 UniProt-annotated genes; *S. cerevisiae* taxon 559292):

| relocate set | reactions | agreement with curation | applied model grows |
|---|---|---|---|
| targeted subset (default `--max-reactions`) | ~250–500 | **74–83 %** | **yes — baseline** |
| entire UniProt-covered genome (`--all`) | ~1170 | **77.5 %** | **no** — see *scale limit* below |

* **~3/4 agreement** is encouraging for an agnostic prior. The remaining gap is where evidence and
  curation legitimately differ: UniProt coverage gaps, multi-located proteins, membrane
  sub-compartments the evidence can't resolve, and reactions the **functionality constraint
  deliberately places against the score**.
* **For a targeted relocate set the result is *functional*** — not just a maximum-likelihood label
  assignment. An essential mitochondrial reaction whose gene UniProt annotates as cytosolic (e.g.
  ERG10 / r_0104) is kept in the mitochondrion because moving it breaks biomass; the localisation
  prior yields to the network. This is the point of the tool.

## Functionality soundness — and a scale limit

The functionality guarantee rests on Big-M flux gating (`|v[r,c]| ≤ ub·x[r,c]`), which is only as
tight as the solver's **integer-feasibility tolerance**. A placement binary rounded to `~tol` still
lets `ub·tol` of *ghost flux* pass through a reaction "placed elsewhere".

* **Per-reaction leak (fixed).** With `ub=1000` and the default tolerance `1e-5`, a single reaction
  could leak `0.01` — more than a typical growth floor — so the MILP would certify a single
  essential reaction in the wrong compartment (it reported optimal while the materialized model
  grew 0). Tightening the tolerance (`integrality_tol`, default `1e-9`) closes this: a targeted
  relocate set now yields a genuinely functional model.
* **Accumulated leak at full-genome scale (open).** Gurobi's `IntFeasTol` floors at `1e-9`, so the
  residual per-binary leak (`ub·1e-9 = 1e-6`) cannot be tightened further. Relocating *every*
  reaction at once (~1170 reactions × 13 alternative compartments) sums to `~0.015` of admissible
  ghost flux — above `min_growth ≈ 0.008` — so the genome-wide stress run reports 77.5 % agreement
  but a model that does not grow. optlang does not expose indicator constraints (the exact,
  leak-free gating), so the robust remedy is **post-solve verify-and-repair** (materialise, FBA,
  re-pin the essential reactions whose move broke growth, re-solve) — tracked as future work.

**Practical guidance:** relocate a *targeted* set of reactions (the normal use), where the guarantee
holds; treat `--all` as a stress benchmark, not a production setting, until verify-and-repair lands.

## Multi-localization (opt-in) at genome scale

With `multi_localization=True` (the sound flux-activity-coupling formulation; a reaction may occupy
several compartments but every extra placement must carry flux), yeast-GEM was used to check three
things: does it stay **tractable**, does it stay **sound** at scale, and does genuine dual-targeting
emerge? Relocate sets are batches of single-compartment cytosol/mito/peroxisome reactions, each gene
scored 1.0 for its true compartment (so the solver has no incentive to multi-localize *unless*
functionality warrants it); the *dual-lure* run additionally scores every cytosolic gene 0.9 in the
mitochondrion to tempt spurious dual-targeting.

| relocate set | reactions | status | time (Gurobi) | recovery | multi-localized | dead placements |
|---|---|---|---|---|---|---|
| true-comp prior | 120 | optimal | ~10 s | 118/120 | 0 | 0 |
| true-comp prior | 300 | optimal | ~37 s | 296/300 | 1 | 0 |
| dual-lure | 36 | optimal | ~14 s | 33/36 | **2 (genuine)** | 0 |

* **Sound at scale.** Across every run, *zero* dead placements: each materialized multi-placement
  carries ≥ `eps_flux` flux (FVA-checked). When scores favour one compartment the solver stays mono
  (no spurious multi-localization); under the dual-lure it dual-localizes only the reactions that can
  carry flux in both, and the model still grows.
* **Tractable.** The larger MILP (activity + home/used binaries) solves to optimality up to ~300
  relocated reactions in well under a minute — a modest constant-factor over mono. Far past that the
  same accumulated-leak scale limit as mono applies; use `time_limit`/`mip_gap` for large batches.
* **`eps_flux` must track the model's flux scale.** A *fixed* activity threshold creates a deadzone
  `(0, eps_flux)`: a reaction forced to carry a tiny but nonzero flux in a compartment cannot clear
  the threshold, and — being *used* — is allowed no inactive home, so the model becomes **infeasible**
  rather than merely mono. A fixed `1e-4` made yeast-GEM (biomass ≈ 0.08, many sub-`1e-4` fluxes)
  infeasible above ~80 relocated reactions. The default now scales to `10 · big_m · integrality_tol`
  (`1e-5`) — just above the ghost-flux floor, below yeast-GEM's meaningful fluxes — and is feasible
  and sound across all the sets above. Lower `eps_flux` further (toward the floor) for models with
  even smaller fluxes.

## Caveats

* UniProt locations are **qualitative** (presence, not ranked probabilities); DeepLoc 2 / MULocDeep
  give graded scores and would rank ambiguous compartments better.
* Membrane sub-compartments (ER membrane, mitochondrial membrane, …) are not in the default map, so
  reactions there are out of scope for this run.
* `transportable=[]` forbids adding transporters; allowing them (or coupling to gap-filling) is a
  different, richer experiment.
* Agreement is **not** an accuracy ceiling — curation itself is imperfect, and a reaction placed
  against its score for functional reasons is a feature, not an error.
