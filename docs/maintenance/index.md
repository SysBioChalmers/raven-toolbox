# Maintenance

Maintainer-facing documentation: the layout of the published data bundles and how to
rebuild and release them.

- **[KEGG data format](kegg_data_format.md)** — layout of the KEGG artefact bundle.
- **[Maintaining KEGG data](maintaining_kegg_data.md)** — building and publishing the KEGG
  artefact releases.
- **[Maintaining binaries](maintaining_binaries.md)** — building and publishing the
  external-binary (BLAST / DIAMOND / HMMER) ZIP releases.
- **[Data & binary manifest](data_manifest.md)** — the shared manifest that lists every
  published artefact / binary (consumed by raven-toolbox and MATLAB RAVEN), where to host
  assets (GitHub Releases vs Zenodo), and the GitHub→Zenodo auto-publish setup.

```{toctree}
:hidden:

kegg_data_format
maintaining_kegg_data
maintaining_binaries
data_manifest
```
