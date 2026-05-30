# De-novo reconstruction

Two independent routes build a draft model for an organism that has no curated GEM yet.

## Homology — {mod}`raven_python.reconstruction.homology`

Transfer reactions from a curated **template** model to the target organism via an ortholog
search.

1. Run the search (or bring your own table):
   {func}`raven_python.reconstruction.homology.run_blast` /
   {func}`raven_python.reconstruction.homology.run_diamond`, then
   {func}`raven_python.reconstruction.homology.make_ortholog_hits` to get the canonical
   `gene × gene` hits DataFrame (bidirectional / best-hits-only policies supported).
2. Draft the model:
   {func}`raven_python.reconstruction.homology.get_model_from_homology` — AST-based GPR
   rewrite, configurable complex policy, returns a {class}`cobra.Model` (plus a
   `HomologyResult`).

The external aligners are resolved from the pinned binary registry
({mod}`raven_python.binaries`).

## KEGG — {mod}`raven_python.reconstruction.kegg`

Draft directly from KEGG orthology, either for a KEGG-listed species or from your own protein
FASTA via HMM search.

- **KEGG species (no FASTA):**
  {func}`raven_python.reconstruction.kegg.get_kegg_model_for_organism`.
- **Your sequences (HMM search):** {func}`raven_python.reconstruction.kegg.assign_kos` →
  {func}`raven_python.reconstruction.kegg.get_kegg_model_from_sequences`. The HMM cut-off
  defaults are calibrated in the
  [KEGG HMM cut-off study](../studies/kegg_hmm_cutoff_calibration.md).

The KEGG artefact bundle (KO tables, reference model, HMM libraries) is fetched by
{mod}`raven_python.data`; building and publishing it is a maintainer task — see
[Maintaining KEGG data](../maintenance/maintaining_kegg_data.md) and the
[KEGG data format](../maintenance/kegg_data_format.md).

After drafting, fill connectivity gaps with the
[gap-filling guide](tasks_and_gapfilling.md).
