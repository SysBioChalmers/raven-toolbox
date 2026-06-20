# Localization parameter benchmarks

Function: `raven_toolbox.localization.predict.predict_localization`

Date: 2026-06-20.

---

## `time_limit` — MILP wall-clock cap

**Parameters tested:** `None` (Python default), `900` (MATLAB default, 15 min)

`predict_localization` solves a MILP to optimally assign reactions to compartments
given localisation prediction scores. MATLAB caps this at 900 seconds (15 min).
Python has no cap.

**Timing measurements (2026-06-20, Gurobi, yeast-GEM):**

| Scenario | Reactions | Wall time |
|---|---|---|
| yeast-GEM pilot (200 reactions) | 200 | 11.6 s |
| yeast-GEM full (linear extrapolation) | 2,682 | ~155 s |
| Human-GEM (linear extrapolation) | ~13,000 | ~750 s |

Note: MILPs do not scale linearly — hard instances with ambiguous or uniform
localization scores can be dramatically slower than these estimates. The linear
extrapolation from 200 reactions should be treated as a lower bound.

**Decision: keep `time_limit=None`.** For yeast-GEM (the primary raven-toolbox
development model), the solve completes in ~2.5 minutes with no cap. For Human-GEM
scale or noisy scores, users should pass `time_limit=900` explicitly. Add a
docstring note:

> For genome-scale models with >5000 gene-associated reactions, or when localization
> scores are ambiguous (many reactions with similar scores across compartments),
> consider setting `time_limit=900` (15 minutes, matching MATLAB) to prevent runaway
> solves.

---

## `transport_cost` — penalty for assigning a reaction to a non-default compartment

**Parameter:** `transport_cost=0.5` (Python and MATLAB)

Not benchmarked. This parameter was introduced in the MATLAB RAVEN implementation
without a published sensitivity analysis. Both implementations use 0.5.

Proposed future test: run `predict_localization` on yeast-GEM at `transport_cost`
∈ {0.1, 0.25, 0.5, 1.0, 2.0} and compare the number of reactions relocated
and whether the result matches the known yeast-GEM compartment assignments.

---

## `default_compartment` — compartment assigned to reactions with no score

**Parameter:** `default_compartment='c'` (Python); required argument in MATLAB

Python provides a better UX by defaulting to cytosol (`'c'`), which is correct
for the vast majority of metabolic reactions in well-curated GEMs. No change
needed; the MATLAB version requires this to be specified, which is a regression
in usability.
