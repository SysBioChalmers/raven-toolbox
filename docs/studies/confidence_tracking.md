# Per-reaction confidence tracking

**Status: the three facets are shipped.** The data model, the notes round-trip (YAML + SBML), and the
`localization`, `equation` and `gene_association` scorers live in
[`raven_toolbox/confidence.py`](../../src/raven_toolbox/confidence.py) (tests in `tests/test_confidence.py`).
The facet set is **closed** — see [§10](#10-what-is-left) for the remaining work.

Every reaction carries a small structured record scoring how well-supported each of its **facets** is,
persisted in the model file, computed from evidence, updated by curation, and consumed by the raven-toolbox
tools. It generalises the localisation-curation work
([curation_priority_signals.md](curation_priority_signals.md)) from "which placements to review" to "how
well-supported is each facet of each reaction".

## 1. The two rules every score obeys

`ReactionConfidence.overall` is the **weakest** facet (`min`), and the record is dropped from the model when
no facet is left. Those two facts, taken together, force the design:

1. **Abstain rather than guess.** A facet that *does not apply* to a reaction is **not written at all**. An
   absent facet is neutral under `min` — exactly as neutral as `1.0` would be — but it claims nothing. An
   exchange reaction is imbalanced by construction, so it gets *no* `equation` facet rather than a perfect
   one. Ask *why* with `equation_exempt()` / `gene_association_exempt()`, which return the reason string.

   This matters more than it sounds. Writing `1.0` for the 469 exempt yeast-GEM reactions would make them
   indistinguishable from the 3617 genuinely verified ones, and "the fraction of this model whose chemistry
   is verified balanced" would read 4086/4102 instead of the true 3617/3633.

2. **A zero is a measurement, never ignorance.** `score == 0.0` means the evidence *contradicts* the model —
   a proven mass imbalance, or zero localisation support at the assigned compartment. Where evidence is
   merely missing or uninterpretable, the score is low but positive. Hence the invariant the tests assert:

   ```
   max(defect scores) < min(ignorance scores)          # 0.1 < 0.2
   ```

   so a reaction with a proven defect always outranks one that is merely unverifiable, and `overall == 0.0`
   is a usable filter for "this model is provably wrong here".

## 2. Facets

| facet | what it scores | evidence source |
|---|---|---|
| `localization` | the compartment assignment | DeepLoc support at the assigned compartment + FBA certification |
| `equation` | mass & charge balance, formula completeness | `get_elemental_balance` + a recomputed charge sum |
| `gene_association` | is there gene evidence, and is it corroborated | GPR presence + a `pubmed` annotation |

Each facet is scored independently, so a model can be annotated one facet at a time and the record grows
incrementally. A `ConfidenceEntry` is a continuous `score` in [0, 1] plus optional provenance: a categorical
`level`, the `basis` evidence, and `method` / `source` / `note`.

## 3. Storage & round-trip

Grounded by direct checks against cobrapy:

- **YAML / JSON models:** a value under `reaction.notes["raven_confidence"]` round-trips losslessly and is
  ignored by cobrapy's own logic. *Verified.*
- **Arbitrary top-level reaction keys are dropped** by cobrapy's serialiser, so `notes` — not custom
  top-level fields — is the channel. *Verified.*
- **SBML models** store `<notes>` as strings only, and a JSON string picks up HTML entities (`&quot;`).
  *Verified.*

So the record is stored as one **JSON string** under `reaction.notes["raven_confidence"]` (with a
`schema_version`) — the same on write for every format. Reading does `json.loads(html.unescape(...))`, which
handles the clean YAML string and the HTML-escaped SBML string alike (and tolerates a raw dict from a
hand-edited YAML model). A confidence-annotated model still loads and solves in **plain cobra unchanged** —
this is a test invariant.

## 4. The `equation` facet

### Exemptions (SBO-driven, no name heuristics)

| SBO / property | label | why exempt |
|---|---|---|
| `reaction.boundary` | exchange / demand / sink | exchanges mass with the environment |
| `SBO:0000629` | biomass production | lumped; "some reactants are used in unrepresented processes" |
| `SBO:0000395` | encapsulating process | "an aggregation of interactions and entities into a single process" (SLIME, pool reactions) |

`SBO:0000630` (ATP maintenance) is deliberately **not** exempt: `ATP + H2O -> ADP + Pi + H` is real chemistry
and must balance. Detecting biomass by *name* is likewise refused — `\bgrowth\b` matches yeast-GEM's `r_4046`,
"non-growth associated maintenance reaction", and a name regex must never silence a chemistry check. When a
model carries **no** reaction SBO terms at all, the scorers warn: they then cannot tell a pseudo-reaction from
a defect.

Note also that `reaction.boundary` is `len(metabolites) == 1`, independent of id and of bounds — so it
catches an exchange reaction that is not named `EX_`, and a blocked `(0, 0)` one.

### Bands

| score | `basis` | meaning |
|---:|---|---|
| 0.0 | `mass-imbalanced` | **defect**: atoms do not conserve |
| 0.1 | `charge-imbalanced` | **defect**, one rung up: usually a protonation-state convention |
| 0.3 | `formula-missing` / `formula-unparseable` / `formula-generic` | **ignorance**: the verdict is impossible, not negative |
| 0.6 | `charge-unknown` | mass proven balanced; a metabolite has no charge |
| 1.0 | `balanced` | **proven**: mass and charge conserve |

Three cobra behaviours the bands exist to survive, each verified by running it:

- **`check_mass_balance() == {}` is ambiguous.** A metabolite with no formula contributes zero atoms
  *silently*, so a reaction can read "balanced" when the truth is that the data is missing — or read
  "imbalanced" purely because a formula is absent. `get_elemental_balance` checks formula coverage first and
  returns `unknown`.
- **`check_mass_balance()` raises on a parenthesised polymer** such as `(C5H8)n` (glycogen, starch), because
  `Metabolite.elements` returns `None` there. Shipped `get_elemental_balance` used to propagate that
  exception; it now reports `unknown`. A present-but-uninterpretable formula is not a crash.
- **`check_mass_balance()["charge"]` is fabricated** when some metabolites have `charge=None`: cobra
  accumulates the sum over only the non-`None` ones, inventing a residual out of half the reaction. The
  charge verdict is therefore recomputed behind an all-charges-present guard.

Generic symbols (`R` groups) get two different treatments, because they mean two different things. A residual
that *lands in* an R group is uninterpretable → `formula-generic`, 0.3. An R group that merely *appears* and
then cancels is recorded in `note` and does not move the score.

Yeast-GEM shows why the distinction has to be drawn that way rather than by demoting every reaction that
touches an R group: 131 of its metabolites carry one, 107 scored reactions touch such a metabolite, and in
**every one of the 107 the R cancels** — not a single residual lands in an R group. 102 balance outright; the
other five fail for reasons that have nothing to do with the R group (two H imbalances, two charge
imbalances, one missing formula). Demoting all 107 would have dragged 2.9% of the model down for a symbol
that never affected a verdict. Cobra's own 114-symbol element table is used for the test, since it contains
`U` and `F`: a hand-rolled element set would read `FULLR2`'s letters as uranium and fluorine.

## 5. The `gene_association` facet

Exempt: everything `equation_exempt` excludes, plus the reactions that legitimately have no catalyst —
`SBO:0000630` (ATP maintenance) and `SBO:0000672` (spontaneous reaction, "no catalyst … is needed to
proceed"). Note the deliberate **asymmetry**: ATP maintenance is `equation`-scored and `gene`-exempt. Real
chemistry, no catalyst.

Transport (`SBO:0000655`) is **not** exempt: the term is defined as movement "mediated by a transporter
protein", so a transport reaction with no transporter gene is a genuine curation gap. (816 of yeast-GEM's 950
gene-less scored reactions are transporters.)

| score | `basis` | ≈ Thiele-Palsson |
|---:|---|---|
| 0.2 | `no-gpr` | 0–1 |
| 0.6 | `gpr` | 2 |
| 0.9 | `gpr+literature` (a `pubmed` annotation) | 3 |

`1.0` is reserved for `mark_curated()`, so an inferred score never ties a curator's call. GPR *shape* —
isozyme count, complex size — is recorded in `note` and never scored from: nothing justifies a number there.

## 6. Measured on yeast-GEM (4102 reactions)

Every number below is the output of `scripts/measure_confidence_facets.py`, which runs the shipped scorers
over the model — not prose arithmetic:

```
           facet             basis  score    level     n
        equation   mass-imbalanced    0.0     none     2
        equation charge-imbalanced    0.1     weak     8
        equation   formula-missing    0.3     weak     6
        equation          balanced    1.0   strong  3617
gene_association            no-gpr    0.2     weak   950
gene_association               gpr    0.6 moderate   749
gene_association    gpr+literature    0.9   strong  1933
```

469 reactions are `equation`-exempt (273 boundary + 195 encapsulating + 1 biomass) and 470 are
`gene_association`-exempt (those plus 1 ATP maintenance), so 3633 and 3632 are scored. **`overall == 0.0`
selects exactly two reactions** — `r_0438` and `r_0439`, the model's only proven mass imbalances — followed by
the eight charge defects. The review queue's top ten is ten real chemistry problems.

### An independent check of the gene rubric

The rubric reads only the GPR and the `pubmed` annotation. It never reads yeast-GEM's own curator-assigned
`Confidence Level` note — which leaves that note free to serve as an *independent* check rather than as an
input. Cross-tabulated over the scored reactions:

| our `basis` | recorded Thiele-Palsson | agreement |
|---|---|---|
| `gpr+literature` | 3 | 1932 / 1933 (99.9%) |
| `gpr` | 2 | 685 / 749 (91.5%) |
| `no-gpr` | 0 or 1 | 905 / 950 (95.3%) |

A rubric derived from two observable facts recovers a human curator's confidence assignment. The
disagreements are informative rather than embarrassing: 43 gene-less reactions carry a recorded confidence of
2 or 3, i.e. the model claims genetic evidence it does not encode. The scorer flags exactly these in `note`
("recorded Confidence Level 3 but no GPR").

### What yeast-GEM cannot test

Three bands never fire on it, and are therefore exercised only by synthetic fixtures: `formula-unparseable`
(no parenthesised formulas), `formula-generic` (no residual lands in an R group), and `charge-unknown` (all
2748 metabolites carry a charge). This is stated rather than hidden — a distribution measured on one model is
not a guarantee about another.

## 7. Integration with the existing tools

- **After `assign_compartments`**, `score_localization_confidence(model, proposal, scores)` writes a
  `localization` confidence per placement from DeepLoc support (+ FBA certification). A reaction whose genes
  are absent from the score table is **not scored**: no measurement was possible, and a `0.0` there would veto
  the reaction's `overall` on the strength of a missing input rather than of evidence. Such reactions still
  surface — through their `gene_association` facet and through `curation_priority`'s `no_evidence` signal.
- **After a curator relocation** (`relocate_reactions`), `mark_curated(reaction)` stamps
  `{score: 1, level: "curated", source: "curator"}` so the decision persists in the model file and the next
  automated pass leaves it untouched. `mark_curated(reaction, facet=...)` pins any facet.
- **`curation_priority`** is the inverse view: high localisation confidence here == low review priority there.
  `confidence_report(model)` is the lowest-confidence-first review queue; `facet_summary(model)` groups by
  `facet` × `basis` and is the audit trail abstention leaves behind.

## 8. API

**Shipped:** `ConfidenceEntry`, `ReactionConfidence`, `get_confidence` / `set_confidence` /
`clear_confidence`, `read_confidence`, `mark_curated`, `equation_exempt`, `gene_association_exempt`,
`score_localization_confidence`, `score_equation_confidence`, `score_gene_association_confidence`,
`confidence_report`, `facet_summary`. Storage lives in `reaction.notes["raven_confidence"]`; there is no
separate save step — the record serialises with the model.

**Planned:** an umbrella `annotate_confidence(model, facets=[...])` — see [§10](#10-what-is-left).

## 9. Standards alignment (for the paper)

Map the categorical `level` onto the established **Thiele & Palsson reconstruction confidence score (0–4)** so
it is familiar to modellers and reviewers, and reference **ECO** (Evidence & Conclusion Ontology) terms where a
facet maps to an evidence class. Two cautions carried forward from the design work:

- Thiele & Palsson's table assigns score **2 to two different evidence classes** (physiological data *and*
  sequence data), so a recorded 2 does not determine an evidence class and must not be mapped to one.
- Assertion-method ECO terms (e.g. "inferred by curator") are not evidence classes and must not be used as
  such. Any ECO id must be checked against the ontology before it reaches the paper — a wrong ontology id is
  worse than an omitted one. The SBO ids used above were each verified against EBI's SBO: `SBO:0000629`
  biomass production, `SBO:0000395` encapsulating process, `SBO:0000630` ATP maintenance, `SBO:0000672`
  spontaneous reaction, `SBO:0000655` transport reaction.

## 10. What is left

The **facet set is closed**: `localization`, `equation` and `gene_association` are shipped, and no further
facet is planned. What remains is finishing the work *around* them.

### 10.1 Wire the facets together

- **`annotate_confidence(model, facets=[...])`** — one call that runs every applicable scorer, instead of
  making a caller know the three scorer names and their argument shapes. `localization` needs a proposal and
  a score table while the other two need only the model, so the umbrella must skip a facet whose inputs are
  absent rather than fail — the same abstain-rather-than-guess rule the scores themselves follow.
- **Let `curation_priority` read the record.** Today it re-derives localisation evidence from scratch and
  cannot see that a curator already settled a placement, so a `mark_curated` reaction keeps surfacing in the
  review queue. Skipping facets at `level == "curated"` closes the loop between the two tools: score → review
  → curate → *stop being asked about it*. This is the single change that makes curation feel finished.
- **Point the SBO precondition at its remedy.** The scorers warn on a model with no SBO terms, but do not say
  that `raven_toolbox.annotation.add_sbo_terms(model)` is the fix. The warning should name it.

### 10.2 Validate beyond one model

Every number in §6 comes from yeast-GEM, and §6 already flags that three bands (`formula-unparseable`,
`formula-generic`, `charge-unknown`) never fire there — they are covered only by synthetic fixtures. A
distribution measured on one model is not a guarantee about another. Running
`scripts/measure_confidence_facets.py` over Human-GEM and a non-curated draft would show whether the bands
are calibrated or merely yeast-shaped, and would exercise the branches yeast cannot reach. The
gene-rubric-vs-`Confidence Level` check only replicates on a model that records that note, so its absence
elsewhere is itself worth reporting.

### 10.3 Standards alignment for the paper

§9 above: map `level` onto Thiele & Palsson 0–4, and attach ECO terms where a facet maps to an evidence
class — with the two cautions recorded there. The SBO half is already done and verified.
