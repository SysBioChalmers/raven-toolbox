# Sub-cellular localisation

{mod}`raven_toolbox.localization` assigns reactions to compartments by MILP — deterministic
(not simulated annealing), predictor-agnostic, and partial-update friendly.

1. **Load predictor/database scores** into the `gene × compartment`
   {class}`raven_toolbox.localization.LocalizationScores` frame:
   {func}`raven_toolbox.localization.load_deeploc` (DeepLoc 2 per-protein CSV),
   {func}`raven_toolbox.localization.load_mulocdeep` (MULocDeep wide table),
   {func}`raven_toolbox.localization.load_compartments` (COMPARTMENTS evidence database), or
   {func}`raven_toolbox.localization.load_uniprot` (curated UniProtKB `Subcellular location`).
   {func}`raven_toolbox.localization.fetch_uniprot_localization` queries the UniProt REST API by
   organism id directly, no manual export needed.
   Pass {data}`raven_toolbox.localization.DEFAULT_COMPARTMENT_MAP` to rename predictor labels
   to your model's compartment ids. raven-toolbox does **not** shell out to the predictor — run
   it separately and feed in its output.

   **Fuse and tune the evidence.** Since no single source is authoritative (two curated sources
   agree only ~90% on yeast-GEM — see the [DeepLoc benchmark](../studies/deeploc_yeast_benchmark.md)),
   {func}`raven_toolbox.localization.combine_scores` merges several `LocalizationScores` into a
   consensus (agreement reinforced). `load_deeploc` also takes `min_confidence=` (drop unreliable
   low-confidence genes — DeepLoc's probability is well calibrated) and `membrane_split={"m":"mm"}`
   (route mitochondrion to its membrane sub-compartment using the transmembrane signal; mito only).
2. **Predict / apply:** {func}`raven_toolbox.localization.predict_localization` is the MILP
   entry point. Pass the set of reactions to relocate (everything else is pinned); extra
   compartments beyond a reaction's primary one pay a `multi_compartment_penalty`. With
   `apply=False` you get a {class}`raven_toolbox.localization.LocalizationProposal` diff to
   inspect before committing; {func}`raven_toolbox.localization.apply_localization` applies a
   result.
3. **Triage for curation (optional):** {func}`raven_toolbox.localization.triage_localization` takes
   the proposal + scores and returns a {class}`raven_toolbox.localization.ReviewReport` ranking the
   genes/reactions whose localisation is shakiest (low confidence, near-ties, source disagreement,
   no evidence, low-trust compartment, multi-localisation), each with a reason — so a curator knows
   where to look. Pass `load_deeploc(..., keep_raw_confidence=True)` so the (strongest) confidence
   signal survives normalisation.

The defaults and accuracy (including a predictor-noise sweep) are validated against curated
yeast-GEM in the
[yeast localization benchmark](../studies/yeast_localization_benchmark.md).

## Preparing input for a sequence predictor (DeepLoc 2.1)

Sequence predictors take a protein **FASTA** and emit the per-protein table the loaders read back.
{func}`raven_toolbox.localization.prepare_deeploc_input` writes that FASTA for your model's genes —
sequences fetched from UniProtKB, headers set to the **gene ids** so the predictor's output lines up
with the model (and with {func}`raven_toolbox.localization.load_deeploc`) without remapping:

```python
from raven_toolbox.localization import prepare_deeploc_input, load_deeploc, DEFAULT_COMPARTMENT_MAP

# yeast-GEM gene ids are ORF names (e.g. YNR001C) = UniProt ordered-locus, the default key
res = prepare_deeploc_input(model, 559292, "deeploc_yeast.fasta")  # 559292 = S. cerevisiae
print(res)                 # "1143/1143 sequences -> deeploc_yeast_001.fasta, ... (0 missing)"
print(res.missing)         # genes with no reviewed UniProt sequence — chase these separately

# ... run DeepLoc 2.1 on res.paths yourself, then:
scores = load_deeploc("deeploc_results.csv", compartment_map=DEFAULT_COMPARTMENT_MAP)
```

**DeepLoc 2.1 has no public batch API**, so run it yourself on the FASTA:
the [web server](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/) (max **500 sequences** per
submission — `prepare_deeploc_input` chunks the FASTA at that limit into `…_001.fasta`, `…_002.fasta`,
…) or the downloadable standalone (no limit; pass `max_records_per_file=None`). Sequences must be
≥ 10 aa (enforced). {func}`raven_toolbox.localization.fetch_protein_sequences` and
{func}`raven_toolbox.localization.write_fasta` are the underlying building blocks if you need finer
control. A ready-to-run script is `scripts/prepare_deeploc_yeast.py`.
