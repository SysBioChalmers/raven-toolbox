# Maintainer scripts

Release-time tooling. Not part of the installed package — run them from a checkout
with raven-python installed (`pip install -e .`). End users never need these.

## `build_kegg_artefacts.py`

Build the publishable KEGG artefact set from an arranged KEGG dump (see
`download_kegg_dump`): the gzipped-YAML reference model, the gzipped-TSV tables,
and (with `--hmms`) the per-domain pressed HMM libraries. Output is laid out ready
to upload as release assets. See [docs/maintaining_kegg_data.md](../docs/maintaining_kegg_data.md).

```bash
python scripts/build_kegg_artefacts.py --keggdb keggdb --out artefacts          # tables + model
python scripts/build_kegg_artefacts.py --keggdb keggdb --out artefacts --hmms --threads 8
```

## `make_registry_snippet.py`

After uploading the files to a release, compute their SHA256 and print the entry
to merge into the runtime registry — `raven_python.data._DATA_REGISTRY` (data) or
`raven_python.binaries._REGISTRY` (binary ZIP bundles). The checksum helper is shared
with the resolvers, so published checksums always match what `ensure_data` /
`ensure_binary` verify.

```bash
# Data artefacts:
python scripts/make_registry_snippet.py data --dataset kegg --version kegg116 \
    --dir artefacts --base-url https://github.com/ORG/raven-python/releases/download/kegg-data-kegg116

# Binary bundle (ZIPs named <bundle>-<version>-<os>-<arch>.zip):
python scripts/make_registry_snippet.py binary --bundle blast --version 2.16.0 \
    --provides blastp makeblastdb --dir zips \
    --base-url https://github.com/ORG/raven-python/releases/download/blast-2.16.0
```
