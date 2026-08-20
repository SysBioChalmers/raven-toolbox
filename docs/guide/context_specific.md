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
in the [parameter-calibration study](../studies/init_param_calibration.md), and the
equivalence to MATLAB RAVEN is established in the
[Human-GEM validation](../studies/humangem_validation.md).

`ftinit()` also takes two opt-in determinism flags, `strict_gap` and `canonical` (both default
off → exact RAVEN behaviour), which pin *which* of the degenerate MILP's many equal-score optima
is returned, giving a more parsimonious and more reproducible extraction. They do not make the
model biologically more accurate, they cost 3–7× build time, and they are **not** a fix for
gene-essentiality reproducibility — measurements and the trade-off are in the
[extraction-determinism study](../studies/ftinit_determinism.md). For reproducible results
across runs, pin the solver stack (raven-toolbox commit + `gurobipy` version) instead.

:::{important}
Genome-scale (f)tINIT MILPs currently require **Gurobi** for tractable solve times; toy and
unit-test problems run on GLPK. See the
[solver benchmark](../studies/init_solver_benchmark.md). Metabolomics-based scoring is the one
piece not yet implemented (raises `NotImplementedError`).
:::
