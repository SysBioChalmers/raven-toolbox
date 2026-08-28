# Sub-cellular localisation

{mod}`raven_toolbox.localization` assigns reactions to compartments by MILP — deterministic
(not simulated annealing), predictor-agnostic, and partial-update friendly.

1. **Load predictor scores** into the `gene × compartment`
   {class}`raven_toolbox.localization.LocalizationScores` frame:
   {func}`raven_toolbox.localization.load_wolfpsort` (WoLF PSORT summary output) or
   {func}`raven_toolbox.localization.load_deeploc` (DeepLoc 2 per-protein CSV). raven-toolbox
   does **not** shell out to the predictor — run it separately and feed in its output.
2. **Predict / apply:** {func}`raven_toolbox.localization.predict_localization` is the MILP
   entry point. Pass the set of reactions to relocate (everything else is pinned); extra
   compartments beyond a reaction's primary one pay a `multi_compartment_penalty`. With
   `apply=False` you get a {class}`raven_toolbox.localization.LocalizationProposal` diff to
   inspect before committing; {func}`raven_toolbox.localization.apply_localization` applies a
   result.

The defaults and accuracy (including a predictor-noise sweep) are validated against curated
yeast-GEM in the
[yeast-GEM localisation benchmark](https://github.com/edkerk/raven-docs/blob/main/docs/parameter-tuning/studies/yeast-localization-benchmark.md)
(raven-docs).
