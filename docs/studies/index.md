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
- **[ftINIT extraction determinism](ftinit_determinism.md)** — what the opt-in `strict_gap` /
  `canonical` flags buy on Human-GEM/DLD1: a 3.7× smaller reaction seed-swing, but *worse*
  gene-essentiality determinism (5 → 19 flips) and a 3–7× build-time cost.
- **[Yeast localization benchmark](yeast_localization_benchmark.md)** —
  `predict_localization` against curated yeast-GEM, with a predictor-noise sweep.
- **[Compartment-assignment redesign](localization_redesign.md)** — the design of
  `assign_compartments` (flux-free placement master + real-FBA certification): keeping the flux model
  out of the placement optimisation so the result stays sound and functional at genome scale.
- **[Curation-priority signals](curation_priority_signals.md)** — the signal catalogue behind
  `curation_priority`: which placements and transports to review by hand, and how the evidence-gated
  score ranks them.
- **[Yeast-GEM validation](yeast_validation.md)** — `assign_compartments` on curated yeast-GEM:
  recovery of the compartmentalisation, a CarveFungi head-to-head (McNemar), and a biological check
  (transporter connectivity, pathway localisation, dual-localised enzymes).
- **[Multi-organism validation](multiorganism_validation.md)** — the same method across four kingdoms
  (yeast, Human-GEM, AraCore, iCre1355), including the chloroplast, with no per-organism changes.
- **[Confidence tracking](confidence_tracking.md)** — per-reaction, multi-facet confidence persisted in
  the model (`raven_toolbox.confidence`): the localisation / equation / gene-association scorers and the
  abstain-vs-zero discipline that keeps `overall == 0.0` meaning "provably wrong", validated against
  yeast-GEM's own curator-assigned confidence scores.
- **[DeepLoc 2.1 yeast-GEM benchmark](deeploc_yeast_benchmark.md)** — real DeepLoc 2.1
  predictions vs curated yeast-GEM compartments; whether membrane-type recovers the
  lumen/membrane split, and an organelle-collapsed run (slow model: 64.6% collapsed).
- **[DeepLoc 2.1 AraCore benchmark](deeploc_aracore_benchmark.md)** — cross-kingdom
  generalisation on an independent *Arabidopsis* plant model (80.3% overall, chloroplast 89.9%).
- **[DeepLoc 2.1 iCre1355 benchmark](deeploc_icre1355_benchmark.md)** — the most training-distant
  test, the green alga *Chlamydomonas* (chloroplast 78%, but cytosol/mito poor on an auto-model).
- **[DeepLoc 2.1 Human-GEM benchmark](deeploc_humangem_benchmark.md)** — human positive control
  (84.7% addressable) and a circularity lesson: 15% of its compartments are DeepLoc-derived.
- **[Localisation finetuning](localization_finetuning.md)** — tuning the DeepLoc-loading
  hyperparameters (`membrane_threshold`, `min_confidence`, triage trust) on the slow yeast run.
- **[predictLocalization head-to-head](predictlocalization_comparison.md)** — the deterministic MILP
  vs RAVEN's stochastic `predictLocalization` on identical inputs (determinism + accuracy + runtime).
- **[CarveFungi analysis](carvefungi_analysis.md)** — how the contemporary carve-a-universal-model
  method works and how our transport-minimising assignment differs.
- **[CarveFungi assignment head-to-head](carvefungi_milp_benchmark.md)** — our transport cost added to
  CarveFungi's *own* carve MILP (CPLEX): ~1.6× fewer transports per reaction (41% fewer) at no
  detectable assignment-accuracy cost; deterministic near-optimal incumbents (the big-M carve does not
  prove optimality).
- **[DeepLoc normalisation benchmark](deeploc_normalisation_benchmark.md)** — normalised
  (top→1.0) vs raw DeepLoc probabilities for compartment assignment on the whole yeast-GEM;
  accuracy-neutral, so normalisation stays the default and `normalise=False` is opt-in.
- **[KEGG HMM cut-off calibration](kegg_hmm_cutoff_calibration.md)** — HMM E-value /
  score-ratio sensitivity for the KEGG HMM-query reconstruction path.
- **[Homology cut-off calibration](homology_cutoff_calibration.md)** — `max_evalue` is
  inert across five orders of magnitude; `min_identity` is the only real lever, and 40 is
  confirmed once precision is weighted above recall (β = 0.5). Defaults unchanged, and the
  curated-GEM ground truth is shown to be circular.

```{toctree}
:hidden:

humangem_validation
init_param_calibration
init_solver_benchmark
ftinit_determinism
yeast_localization_benchmark
deeploc_yeast_benchmark
deeploc_aracore_benchmark
deeploc_icre1355_benchmark
deeploc_humangem_benchmark
localization_finetuning
predictlocalization_comparison
carvefungi_analysis
carvefungi_milp_benchmark
deeploc_normalisation_benchmark
kegg_hmm_cutoff_calibration
homology_cutoff_calibration
localization_redesign
curation_priority_signals
yeast_validation
multiorganism_validation
confidence_tracking
```
