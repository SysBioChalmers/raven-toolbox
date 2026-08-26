# Studies & validation

Empirical validation runs and parameter calibrations that back raven-toolbox's defaults and
its equivalence claims against MATLAB RAVEN.

- **[Human-GEM validation](humangem_validation.md)** — raven-toolbox ftINIT vs MATLAB RAVEN
  on 5 Hart2015 cell lines (Jaccard 0.975–0.980).
- **[(f)tINIT parameter calibration](init_param_calibration.md)** — clean-data calibration
  plus input-robustness study (`mip_gap` / `big_m` / `force_on` / `eps` / `prod_weight` /
  scaling sweeps; dropout / noise / downsample robustness).
- **[(f)tINIT solver benchmark](init_solver_benchmark.md)** — Gurobi vs HiGHS vs GLPK on
  genome-scale ftINIT.
- **[Yeast localization benchmark](yeast_localization_benchmark.md)** —
  `predict_localization` against curated yeast-GEM, with a predictor-noise sweep.
- **[KEGG HMM cut-off calibration](kegg_hmm_cutoff_calibration.md)** — HMM E-value /
  score-ratio sensitivity for the KEGG HMM-query reconstruction path.
- **[Homology cut-off calibration](homology_cutoff_calibration.md)** — `get_model_from_homology`
  thresholds against independent KEGG/OMA ortholog references, across a relatedness series;
  copied from `develop` PR #92 (not yet ported to this branch's code).
- **[Sampling convergence calibration](sampling_convergence_calibration.md)** — between-chain
  Gelman-Rubin R-hat for ACHR sampling, complementing the existing single-chain ESS result;
  in progress.

```{toctree}
:hidden:

humangem_validation
init_param_calibration
init_solver_benchmark
yeast_localization_benchmark
kegg_hmm_cutoff_calibration
homology_cutoff_calibration
sampling_convergence_calibration
```
