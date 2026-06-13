# Omics integration

{mod}`raven_toolbox.omics` ingests **Human Protein Atlas** data and turns it into the gene
scores that drive context-specific extraction.

- **Proteomics:** {func}`raven_toolbox.omics.parse_hpa` →
  {func}`raven_toolbox.omics.hpa_gene_scores`.
- **RNA-seq:** {func}`raven_toolbox.omics.parse_hpa_rna` →
  {func}`raven_toolbox.omics.rna_gene_scores`.

Both return tidy pandas DataFrames, and the scoring adapters reuse
{func}`raven_toolbox.init.score_reactions_from_genes` (a single source of truth for the GPR
walk), so omics-derived scores plug straight into
{func}`raven_toolbox.init.ftinit` / {func}`raven_toolbox.init.get_init_model` — see the
[context-specific modeling guide](context_specific.md).

`HPA_LEVEL_SCORES` exposes the categorical-level → score mapping used for the proteomics
expression levels.
