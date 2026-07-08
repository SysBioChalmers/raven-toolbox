# Per-reaction confidence tracking

**Status: P1 shipped** — the data model, notes round-trip (YAML + SBML), and the ``localization`` scorer
live in [`raven_toolbox/confidence.py`](../../src/raven_toolbox/confidence.py) (tests in
`tests/test_confidence.py`). The ``equation`` / ``gene_association`` / ``reversibility`` facets (P2-P3)
below remain planned. A design for attaching structured, multi-dimensional confidence scores to every
reaction in a genome-scale model, persisted in the model file, computed from evidence, updated by
curation, and consumed by the raven-toolbox tools. It generalises the localisation-curation work
([curation_priority_signals.md](curation_priority_signals.md)) from "which placements to review" to
"how well-supported is each facet of each reaction".

## 1. Concept

Every reaction carries a small structured record scoring how well-supported each of its **facets** is:

| facet | what it scores | cheap evidence source |
|---|---|---|
| `localization` | the compartment assignment | DeepLoc support + FBA certification + `curation_priority` (already built) |
| `equation` | mass & charge balance, formula completeness | `cobra.Reaction.check_mass_balance()` |
| `gene_association` | is there gene evidence; experimental vs inferred/orthology | GPR presence & provenance, DeepLoc coverage |
| `reversibility` | are the bounds thermodynamically justified | ΔG hook / FVA-attainable direction / database directionality |
| *(extensible)* | `subsystem`, `EC/annotation`, `presence` (should the reaction exist) | — |

Each facet is scored independently, so a model can be annotated one facet at a time and the record grows
incrementally.

## 2. Data model

```
ConfidenceEntry:
  score:   float 0-1        # continuous confidence, for ranking
  level:   str (optional)   # categorical: "curated" | "strong" | "weak" | "none"
  basis:   str              # evidence: "deeploc" | "fba-certified" | "mass-balanced" | "orthology" | ...
  method:  str (optional)   # the function/version that produced it
  source:  str (optional)   # "auto" | "curator:<id>" | "database:<name>"
  note:    str (optional)   # e.g. which curation flags fired, or the imbalance
  updated: str (optional)   # ISO date (passed in; not generated inside pure code)

ReactionConfidence:
  { facet_name: ConfidenceEntry, ... }
  overall: float            # derived aggregate (e.g. min, or a weighted mean)
```

Represented as Python dataclasses; serialised to a plain nested dict.

## 3. Storage & round-trip

Grounded by direct checks against cobrapy:

- **YAML / JSON models:** a nested dict under `reaction.notes["raven_confidence"]` round-trips
  **losslessly** through cobrapy and is ignored by cobrapy's own logic. *Verified.*
- **Arbitrary top-level reaction keys are dropped** by cobrapy's serialiser (only its fixed schema is
  written), so `notes` — not custom top-level fields — is the channel. *Verified.*
- **SBML models** (e.g. yeast-GEM) store `<notes>` as strings only: a nested dict `str()`s to
  Python-repr (single-quoted, not JSON) and a JSON string picks up HTML entities (`&quot;`). *Verified.*
  So SBML needs a **canonical JSON string** under one key plus HTML-entity-aware read/write helpers.

Design (as shipped): store the whole record as one **JSON string** under
`reaction.notes["raven_confidence"]` (with a `schema_version`) — the same on write for every format. Read
with `json.loads(html.unescape(...))`, which handles the clean YAML string and the HTML-escaped SBML
string alike (and also tolerates a raw dict from a hand-edited YAML model). A confidence-annotated model
still loads and solves in **plain cobra unchanged** — this is a test invariant.

## 4. Scorers (each returns a `ConfidenceEntry`; all cheap)

- **localization** — directly from the existing work: `1 - normalised(curation_priority)` (or
  evidence-support × certified). A curator relocation ([relocate_reactions](../../src/raven_toolbox/localization/relocate.py))
  stamps `level="curated", score=1.0`. This ties the whole curation pipeline into a persisted score.
- **equation** — `check_mass_balance()` empty ⇒ high; an imbalance is recorded as `basis`; missing
  metabolite formulas are flagged.
- **gene_association** — no gene ⇒ low; gene present ⇒ medium; experimentally evidenced / high-orthology
  ⇒ high.
- **reversibility** — v1 heuristic (assigned bounds vs FVA-attainable direction), with a ΔG hook for
  later.

## 5. Integration with the existing tools (the payoff)

P1 ships the hooks; a caller composes them with the assignment pipeline:

- **After `assign_compartments`**, `score_localization_confidence(model, proposal, scores)` writes an
  initial `localization` confidence per placement from DeepLoc support (+ FBA certification), so a
  freshly-assigned model is already annotated.
- **After a curator relocation** (`relocate_reactions`), `mark_curated(reaction)` stamps a `localization`
  confidence of `{score:1, level:"curated", source:"curator"}`, so the decision persists in the model
  file and the next automated scoring pass leaves it untouched.
- **`curation_priority`** is the inverse view: high localisation confidence here == low review priority
  there. (Having `curation_priority` *read* confidence to skip curator-verified placements is P2 wiring.)

## 6. API — the `raven_toolbox.confidence` module

**Shipped (P1):** `ConfidenceEntry`, `ReactionConfidence`, `get_confidence`/`set_confidence`/
`clear_confidence`, `read_confidence(model)`, `mark_curated`, `score_localization_confidence`, and
`confidence_report(model) -> DataFrame` (reaction × facet matrix + `overall`, lowest-confidence first —
a natural companion to `curation_priority`). Storage lives in `reaction.notes["raven_confidence"]`; there
is no separate save step — the record serialises with the model.

**Planned (P2-P3):** `score_mass_balance_confidence`, `score_gene_association_confidence`,
`score_reversibility_confidence`, and an umbrella `annotate_confidence(model, types=[...])`.

## 7. Standards alignment (for the paper)

Map the categorical `level` to the established **Thiele & Palsson reconstruction confidence score
(0-4)** so it is familiar to modellers and reviewers, and reference **ECO** (Evidence & Conclusion
Ontology) / **SBO** terms where a facet maps to an evidence class — interoperable rather than bespoke,
while keeping the continuous 0-1 for ranking.

## 8. Phasing

- **P1 — foundation + the facet we already have:** data model + storage/round-trip helpers (YAML **and**
  SBML) + the `localization` scorer wired to `relocate`/`assign`/`curation` + tests. Ships immediate
  value (persists the curation work).
- **P2 — cheap structural facets:** `equation` (mass/charge balance) + `gene_association` (model-only).
- **P3 — reversibility + aggregation + reporting:** reversibility heuristic (ΔG hook), the aggregate
  `overall`, `confidence_report`, lowest-confidence ranking.
- **P4 — paper:** Thiele-Palsson / ECO / SBO mapping and documentation.

## 9. Open decisions

1. **Scale** — recommended: continuous 0-1 per facet **plus** an optional categorical `level` mapped to
   Thiele-Palsson 0-4. (Alt: 0-4 only — simpler, coarser for ranking.)
2. **v1 breadth** — recommended: P1 first (localisation, leveraging what exists), then structural facets.
   (Alt: build the full multi-facet skeleton up front.)
3. **Storage layout** — recommended: one JSON blob under `notes["raven_confidence"]` (robust across
   formats). (Alt: flat scalar keys per facet — SBML-native but verbose and structure-less.)
