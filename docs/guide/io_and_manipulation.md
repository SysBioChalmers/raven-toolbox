# I/O & model manipulation

## I/O — {mod}`raven_toolbox.io`

raven-toolbox keeps everything as a {class}`cobra.Model`, so cobra's SBML I/O works
unchanged. On top of that it adds the RAVEN-specific formats:

- {func}`raven_toolbox.io.read_yaml_model` / {func}`raven_toolbox.io.write_yaml_model` —
  cobra-standard YAML (the `!!omap` layout), transparently handling `.yml.gz`. RAVEN-only and
  GECKO `ec-*` side-fields are preserved on each entry's `notes` so a read→write round-trip is
  lossless. The full schema (top-level layout, field order, quoting rules, the GECKO
  `ec-*` and `metaData` extensions) is documented in
  [the YAML model format reference](https://github.com/edkerk/raven-docs/blob/main/docs/yaml-format.md)
  (raven-docs — it also covers RAVEN MATLAB's writer/reader).
- {func}`raven_toolbox.io.export_to_excel` — the RAVEN 5-sheet workbook (RXNS / METS / COMPS /
  GENES / MODEL). Requires the `excel` extra. Excel **import** is intentionally not provided.
- {func}`raven_toolbox.io.export_for_git` — the Standard-GEM repository layout
  (`model/<fmt>/…`), for version-controlled model repos.

## Manipulation — {mod}`raven_toolbox.manipulation`

Structural transforms that cobra does not cover cleanly:

- **Build / copy reactions:** {func}`raven_toolbox.manipulation.add_reactions_from_equations`
  (equation-string → reactions, matching metabolites by id / name / `name[comp]`),
  {func}`raven_toolbox.manipulation.add_transport_reactions`,
  {func}`raven_toolbox.manipulation.add_reactions_from_model`.
- **Merge:** {func}`raven_toolbox.manipulation.merge_models` (N-model merge, unify metabolites
  by `name[comp]`).
- **GPR / bounds:** {func}`raven_toolbox.manipulation.change_gene_reaction_rules`,
  {func}`raven_toolbox.manipulation.change_reaction_equations`,
  {func}`raven_toolbox.manipulation.set_variance_bounds`.
- **Simplify:** {func}`raven_toolbox.manipulation.remove_dead_end_reactions`,
  `remove_duplicate_reactions`, `constrain_reversible_reactions`, `group_linear_reactions`.
- **Topology:** {func}`raven_toolbox.manipulation.convert_to_irreversible`,
  {func}`raven_toolbox.manipulation.expand_model` (split isozyme OR-GPRs),
  and the compartment helpers `merge_compartments` / `copy_to_compartment`.

See the [migration map](../reference/migration.md) for which RAVEN functions these
correspond to (and which became cobra one-liners), and the full
[`io`](../reference/api/io.md) / [`manipulation`](../reference/api/manipulation.md) API.
