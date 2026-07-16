# Maintainer scripts

Release-time tooling. Not part of the installed package — run them from a checkout
with raven-toolbox installed (`pip install -e .`). End users never need these.

The full publish workflow (build → upload → manifest → sync) is documented in
[docs/maintenance/artefact_hosting.md](../docs/maintenance/artefact_hosting.md).

## `build_binary_bundles.py`

Build the per-platform binary ZIPs (BLAST+/DIAMOND/HMMER) from RAVEN's vetted
`software/` binaries (pinned commits) into `dist/binaries/`, with checksums and
provenance. See [docs/maintenance/maintaining_binaries.md](../docs/maintenance/maintaining_binaries.md).

```bash
python scripts/build_binary_bundles.py        # -> dist/binaries/*.zip (+ checksums, PROVENANCE)
```

## `publish_to_raven_data.py`

Upload release assets to the [`raven-data`](https://github.com/SysBioChalmers/raven-data)
repo with `gh`, idempotently (immutable per-version tags; skips assets already present).

```bash
python scripts/publish_to_raven_data.py binaries --dir dist/binaries
python scripts/publish_to_raven_data.py release --tag kegg118 --dir artefacts
python scripts/publish_to_raven_data.py --dry-run release --tag manifest-v1 data/manifest.json
```

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
to merge into the runtime registry — `raven_toolbox.data._DATA_REGISTRY` (data) or
`raven_toolbox.binaries._REGISTRY` (binary ZIP bundles). The checksum helper is shared
with the resolvers, so published checksums always match what `ensure_data` /
`ensure_binary` verify.

```bash
# Data artefacts (--tag is the release tag; the asset URL is built from it):
python scripts/make_registry_snippet.py data --dataset kegg --version kegg116 \
    --dir artefacts --tag v0.3.0

# Binary bundle (ZIPs named <bundle>-<version>-<os>-<arch>.zip):
python scripts/make_registry_snippet.py binary --bundle blast --version 2.16.0 \
    --provides blastp makeblastdb --dir zips --tag blast-2.16.0
```
