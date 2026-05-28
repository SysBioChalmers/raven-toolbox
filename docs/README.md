# Documentation

Start with the [top-level README](../README.md). The docs are organised as:

## RAVEN ↔ ravengem reference

* **[raven_migration.md](raven_migration.md)** — function-by-function map from MATLAB RAVEN
  to ravengem (and cobrapy where appropriate). Read this if you're coming from RAVEN.
* **[../IMPROVEMENTS.md](../IMPROVEMENTS.md)** — improvements ravengem makes over RAVEN that
  are also candidates to back-port into the MATLAB toolbox.

## Open work

* **[todo.md](todo.md)** — what's still on the books.
* **[known_issues.md](known_issues.md)** — low-priority backlog from the full-codebase
  review. None affects correctness on well-formed inputs.

## Empirical studies & calibrations

* **[humangem_validation.md](humangem_validation.md)** — ravengem ftINIT vs MATLAB RAVEN on
  5 Hart2015 cell lines (Jaccard 0.975–0.980).
* **[init_param_calibration.md](init_param_calibration.md)** — clean-data calibration +
  input-robustness study for (f)tINIT (sweeps of `mip_gap` / `big_m` / `force_on` / `eps` /
  `prod_weight` / scaling; dropout / noise / downsample robustness).
* **[init_solver_benchmark.md](init_solver_benchmark.md)** — Gurobi vs HiGHS vs GLPK on
  genome-scale ftINIT.
* **[kegg_hmm_cutoff_calibration.md](kegg_hmm_cutoff_calibration.md)** — HMM E-value /
  score-ratio sensitivity for the KEGG HMM-query reconstruction path.

## Data formats & maintenance

* **[kegg_data_format.md](kegg_data_format.md)** — layout of the KEGG artefact bundle.
* **[maintaining_kegg_data.md](maintaining_kegg_data.md)** — building and publishing the
  KEGG artefact releases.
* **[maintaining_binaries.md](maintaining_binaries.md)** — building and publishing the
  external-binary (BLAST/DIAMOND/HMMER) ZIP releases.

## Archive

Historical design notes (preserved for reference but no longer part of the user-facing
documentation):

* **[archive/ftinit_review_and_plan.md](archive/ftinit_review_and_plan.md)** — the
  pre-implementation critical review of `ftINIT` that drove the Phase 4d port plan.
* **[archive/localization_design.md](archive/localization_design.md)** — the
  pre-implementation critical review of `predictLocalization` that drove the Phase 7
  redesign (caller-passed relocate set, MILP instead of SA, partial-update mode).
* **[archive/plan_get_model_from_homology.md](archive/plan_get_model_from_homology.md)** —
  early planning notes for the homology-based reconstruction.
