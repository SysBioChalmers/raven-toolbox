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
- **[DeepLoc 2.1 yeast-GEM benchmark](deeploc_yeast_benchmark.md)** — real DeepLoc 2.1
  predictions vs curated yeast-GEM compartments; whether membrane-type recovers the
  lumen/membrane split, and an organelle-collapsed run (slow model: 64.6% collapsed).
- **[DeepLoc normalisation benchmark](deeploc_normalisation_benchmark.md)** — normalised
  (top→1.0) vs raw DeepLoc probabilities for compartment assignment on the whole yeast-GEM;
  accuracy-neutral, so normalisation stays the default and `normalise=False` is opt-in.
- **[Localisation finetuning](localization_finetuning.md)** — tuning the DeepLoc-loading
  hyperparameters (`membrane_threshold`, `min_confidence`, triage trust) on the slow yeast run.
- **[KEGG HMM cut-off calibration](kegg_hmm_cutoff_calibration.md)** — HMM E-value /
  score-ratio sensitivity for the KEGG HMM-query reconstruction path.

```{toctree}
:hidden:

humangem_validation
init_param_calibration
init_solver_benchmark
yeast_localization_benchmark
deeploc_yeast_benchmark
deeploc_normalisation_benchmark
localization_finetuning
kegg_hmm_cutoff_calibration
```
