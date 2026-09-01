# Context-specific modeling (tINIT / ftINIT)

Extract a tissue- or condition-specific model from a reference GEM plus gene scores derived
from omics data. Two algorithms are provided in {mod}`raven_toolbox.init`.

## Scoring

Gene scores drive both algorithms. Build them from expression with
{func}`raven_toolbox.init.gene_scores_from_expression` and turn them into reaction scores via
{func}`raven_toolbox.init.score_reactions_from_genes` (a GPR walk shared with the omics
adapters — see the [omics guide](omics.md)).

## tINIT

- {func}`raven_toolbox.init.run_init` — the classic INIT MILP (rewritten in optlang).
- {func}`raven_toolbox.init.get_init_model` — the full tINIT pipeline (dead-end removal →
  `run_init`).

## ftINIT (faster, staged)

- {func}`raven_toolbox.init.run_ftinit` — the single-step ftINIT MILP (continuous indicators
  for positive-score reactions; binaries only on negatives — the speedup over `run_init`).
- {func}`raven_toolbox.init.ftinit` — the full pipeline:
  {func}`raven_toolbox.init.prep_init_model` → staged `run_ftinit` →
  {func}`raven_toolbox.init.fill_tasks` → {func}`raven_toolbox.init.remove_low_score_genes`.

## Tasks and defaults

ftINIT's task layer keeps essential metabolic tasks feasible; define and check tasks with the
[tasks guide](tasks_and_gapfilling.md). The parameter defaults (`mip_gap`, `big_m`,
`force_on`, `eps`, `prod_weight`, scaling) and their robustness to noisy input are calibrated
in the
[INIT parameter calibration study](https://github.com/edkerk/raven-docs/blob/main/docs/parameter-tuning/studies/init-param-calibration.md),
and the equivalence to MATLAB RAVEN is established in the
[Human-GEM validation](https://github.com/edkerk/raven-docs/blob/main/docs/parameter-tuning/studies/humangem-validation.md)
(both on raven-docs).

`ftinit()` also takes two opt-in reproducibility parameters, both default off → exact RAVEN
behaviour. `resolve_ties=True` pins *which* of the degenerate MILP's many equal-score optima is
returned, making repeated builds identical and halving the seed-to-seed spread in predicted
essential genes. `prove_abs_gap=1.0` proves each step to a fixed absolute MIP gap; it is worth
setting because the default gap escalation returns a measurably **suboptimal** extraction. Do
not set it tighter — below ~1.0 it stops being provable at genome scale and returns the same
model anyway. Neither makes the extraction *stable* under a curated template: re-extraction can
move genes unrelated to the edit, so compare against the edit applied to the extracted model.
Measurements in the
[ftINIT reproducibility study](https://github.com/edkerk/raven-docs/blob/main/docs/parameter-tuning/studies/ftinit-determinism.md)
(raven-docs). Pinning the solver stack (raven-toolbox commit + `gurobipy` version) remains the
zero-cost lever for run-to-run identity.

:::{important}
Genome-scale (f)tINIT MILPs currently require **Gurobi** for tractable solve times; toy and
unit-test problems run on GLPK. See the
[INIT solver benchmark](https://github.com/edkerk/raven-docs/blob/main/docs/parameter-tuning/studies/init-solver-benchmark.md)
(raven-docs). Metabolomics-based scoring is the one piece not yet implemented (raises
`NotImplementedError`).
:::
