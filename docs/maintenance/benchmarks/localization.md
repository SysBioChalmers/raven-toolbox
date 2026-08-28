# Localization parameter benchmarks

Function: `raven_toolbox.localization.predict.predict_localization`

Date: 2026-06-20.

---

## `time_limit` — MILP wall-clock cap

**Parameters tested:** `None`, `900` (both Python; the MATLAB side was not run for
this comparison — see the correction below)

`predict_localization` solves a MILP to optimally assign reactions to compartments
given localisation prediction scores. Python has no cap by default.

**Correction (2026-08-28):** this section originally described MATLAB's
`predictLocalization` as capping "this" (the same MILP) at 900 seconds, implying a
direct like-for-like comparison. Reading `core/predictLocalization.m` on
`origin/develop` directly shows that's not accurate: MATLAB solves the placement
problem with **simulated annealing**, not a MILP, and its `maxTime` parameter
(default 15 minutes = 900 s — same number, different meaning) is how long that
stochastic search runs, not a solver cutoff on a bounded optimum. More time in
simulated annealing generally improves the answer; more time on a MILP up to
`time_limit` tightens a proven optimality gap. They aren't the same kind of
parameter, so "matching MATLAB" below is a numeric coincidence worth keeping as a
starting point, not a claim that the two now behave equivalently.

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
scale or noisy scores, users should pass `time_limit=900` explicitly — 900 s is kept
as the suggested starting point since it's a real number from practice, but the
docstring note should no longer call it "matching MATLAB", since MATLAB's number
governs a different kind of process (see the correction above):

> For genome-scale models with >5000 gene-associated reactions, or when localization
> scores are ambiguous (many reactions with similar scores across compartments),
> consider setting `time_limit=900` (15 minutes) to prevent runaway solves.

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
