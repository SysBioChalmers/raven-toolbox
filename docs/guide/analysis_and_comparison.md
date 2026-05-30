# Analysis & comparison

## Analyses — {mod}`raven_python.analysis`

RAVEN analyses that are not in cobrapy's core:

- {func}`raven_python.analysis.reporter_metabolites` — Reporter Metabolites, an
  around-metabolite gene-score test. raven-python uses an exact closed-form background in
  place of RAVEN's Monte-Carlo.
- {func}`raven_python.analysis.fseof` — Flux Scanning based on Enforced Objective Flux:
  regression slope + correlation, with amplify / knockdown / knockout classes and gene
  aggregation.
- {func}`raven_python.analysis.random_sampling` — random-objective flux sampling (wraps
  cobra's samplers); {func}`raven_python.analysis.find_good_reactions` is the companion
  screen.

## Comparison — {mod}`raven_python.comparison`

{func}`raven_python.comparison.compare_models` compares any number of models and returns tidy
DataFrames (reaction / metabolite / gene / subsystem presence, pairwise Jaccard, and an
optional `check_tasks` pass/fail matrix). Plotting and tSNE/MDS are deliberately left out —
they are one-liners in seaborn / scikit-learn on the returned frames.
