# Maintenance

Maintainer-facing documentation: the layout of the published data bundles and how to
rebuild and release them.

- **[Artefact hosting & publishing](artefact_hosting.md)** — **start here.** Where all
  downloadables live (the `raven-data` repo), the versioning model, and the exact
  build → publish → sync → PR workflow.
- **[KEGG data format](kegg_data_format.md)** — layout of the KEGG artefact bundle.
- **[Maintaining KEGG data](maintaining_kegg_data.md)** — building the KEGG artefacts.
- **[Maintaining binaries](maintaining_binaries.md)** — building the external-binary
  (BLAST / DIAMOND / HMMER) ZIPs and the per-platform / licensing matrix.
- **[Data & binary manifest](data_manifest.md)** — the shared manifest format
  (consumed by raven-toolbox and MATLAB RAVEN) and the optional GitHub→Zenodo mirror.
- **[Parameter defaults — inventory and evaluation plan](parameter_defaults.md)** — full
  inventory of optional parameters with current defaults and a methodology for deciding
  whether each default is well-chosen.

```{toctree}
:hidden:

artefact_hosting
kegg_data_format
maintaining_kegg_data
maintaining_binaries
data_manifest
parameter_defaults
```
