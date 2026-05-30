# Omics integration

{mod}`raven_python.omics` ingests **Human Protein Atlas** data and turns it into the gene
scores that drive context-specific extraction.

- **Proteomics:** {func}`raven_python.omics.parse_hpa` →
  {func}`raven_python.omics.hpa_gene_scores`.
- **RNA-seq:** {func}`raven_python.omics.parse_hpa_rna` →
  {func}`raven_python.omics.rna_gene_scores`.

Both return tidy pandas DataFrames, and the scoring adapters reuse
{func}`raven_python.init.score_reactions_from_genes` (a single source of truth for the GPR
walk), so omics-derived scores plug straight into
{func}`raven_python.init.ftinit` / {func}`raven_python.init.get_init_model` — see the
[context-specific modeling guide](context_specific.md).

`HPA_LEVEL_SCORES` exposes the categorical-level → score mapping used for the proteomics
expression levels.
