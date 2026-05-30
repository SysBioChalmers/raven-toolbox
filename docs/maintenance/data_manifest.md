# Data & binary manifest

Large artefacts (KEGG tables / HMMs, template models) and external-binary bundles
(BLAST / DIAMOND / HMMER) are **not** committed to the code repository. They are published
as downloadable assets and described by a single, language-agnostic **manifest** that both
raven-python and MATLAB RAVEN read. Every file carries a **SHA256**, so consumers verify
integrity after download.

- Format: [`data/manifest.schema.json`](https://github.com/SysBioChalmers/raven-python/blob/develop/data/manifest.schema.json) (JSON Schema)
- Worked example: [`data/manifest.example.json`](https://github.com/SysBioChalmers/raven-python/blob/develop/data/manifest.example.json)
- Live manifest: [`data/manifest.json`](https://github.com/SysBioChalmers/raven-python/blob/develop/data/manifest.json) (empty until assets are published)

The manifest is a superset of the two runtime registries:

| Manifest section | Runtime registry |
| --- | --- |
| `data` | {data}`raven_python.data._DATA_REGISTRY` |
| `binaries` | `raven_python.binaries._REGISTRY` |

```json
{
  "manifest_version": 1,
  "data":     { "<dataset>": { "version": "...", "doi": "...", "files":     { "<name>": {"url": "...", "sha256": "...", "bytes": 0} } } },
  "binaries": { "<bundle>":  { "version": "...", "provides": ["..."], "platforms": { "<os>-<arch>": {"url": "...", "sha256": "...", "bytes": 0} } } }
}
```

## Consuming it — Python

Point raven-python at a manifest and the resolvers populate themselves on first use,
verifying each download's checksum:

```bash
export RAVEN_PYTHON_MANIFEST=https://github.com/SysBioChalmers/raven-data/releases/download/manifest-v1/manifest.json
```

```python
from raven_python import manifest
manifest.load_into_registries()           # or load_into_registries("/path/or/url")
# now data.ensure_kegg_data() / binaries.ensure_binary("diamond") resolve from the manifest
```

If `RAVEN_PYTHON_MANIFEST` is set, `data.ensure_*` and `binaries.ensure_binary` load it
lazily — no explicit call needed.

## Consuming it — MATLAB

The same JSON is trivial to read from MATLAB (`webread` + `jsondecode`), download
(`websave`), and verify (Java's `MessageDigest`, always available in MATLAB):

```matlab
function file = ensureDataFile(manifestUrl, dataset, name, cacheDir)
    m = jsondecode(webread(manifestUrl, weboptions('ContentType','text')));
    entry = m.data.(dataset).files.(matlab.lang.makeValidName(name));
    file = fullfile(cacheDir, name);
    if ~isfile(file)
        websave(file, entry.url);
    end
    assert(strcmp(sha256(file), entry.sha256), 'SHA256 mismatch for %s', name);
end

function hex = sha256(file)
    fid = fopen(file, 'r'); raw = fread(fid, Inf, '*uint8'); fclose(fid);
    md = java.security.MessageDigest.getInstance('SHA-256');
    md.update(raw);
    hex = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'))', 1, []));
end
```

## Publishing — generating manifest entries

After uploading a release's files, add/update an entry with the maintainer script
([`scripts/make_registry_snippet.py`](https://github.com/SysBioChalmers/raven-python/blob/develop/scripts/make_registry_snippet.py)),
which computes each SHA256 and byte size:

```bash
python scripts/make_registry_snippet.py manifest --manifest data/manifest.json \
    --target data --dataset kegg --version kegg116 --dir artefacts \
    --base-url https://github.com/SysBioChalmers/raven-data/releases/download/kegg-kegg116 \
    --doi 10.5281/zenodo.0000000 --source https://zenodo.org/records/0000000

python scripts/make_registry_snippet.py manifest --manifest data/manifest.json \
    --target binary --bundle diamond --version 2.1.9 --provides diamond --dir zips \
    --base-url https://github.com/SysBioChalmers/raven-data/releases/download/diamond-2.1.9 \
    --license GPL-3.0-only
```

## Where to host: GitHub Releases vs Zenodo

Both are just URLs in the manifest, so consumers don't care — choose per asset:

- **GitHub Releases** — simplest, free, language-agnostic, up to ~2 GB per file. Good default,
  and you're already on GitHub for the code.
- **Zenodo** — adds a citable **DOI**, long-term archival, and handles files larger than 2 GB
  (up to 50 GB/record). Right for the KEGG HMM bundle and anything you want citable.

### Auto-publishing to Zenodo from GitHub

:::{important}
The **native GitHub↔Zenodo integration** (flip a switch, publish a Release → DOI) archives
the **repository source zipball** at the tag — it does **not** capture files attached to the
Release. So it only works for assets *committed into the repo*, which defeats the purpose for
multi-GB binaries. Use it for a *software* DOI, not for the data assets.
:::

For the data assets, keep everything GitHub-driven with a small **GitHub Action** that, on
release, uploads the assets to Zenodo via its REST API (e.g. [`zenodraft`](https://github.com/zenodraft/zenodraft)).
You cut a normal GitHub Release with the files attached; the Action mirrors them to Zenodo and
mints a new version DOI. Drop this in the data repo as `.github/workflows/zenodo.yml`:

```yaml
name: Mirror release assets to Zenodo
on:
  release:
    types: [published]
jobs:
  zenodo:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-node@v4
        with: { node-version: "20" }
      - name: Download this release's assets
        run: gh release download "${{ github.event.release.tag_name }}" --dir assets
        env: { GH_TOKEN: "${{ github.token }}" }
      - name: Deposit a new version on Zenodo
        run: npx zenodraft@latest version create --publish ${{ vars.ZENODO_CONCEPT_DOI }} assets/*
        env: { ZENODO_ACCESS_TOKEN: "${{ secrets.ZENODO_TOKEN }}" }
```

Then record the resulting DOI in the manifest via the `--doi` flag above. Net result: you only
ever interact with GitHub Releases; Zenodo archiving + DOIs happen automatically.

## Per-asset recommendations

| Asset | Home | Notes |
| --- | --- | --- |
| **Software binaries** (BLAST / DIAMOND / HMMER) | **bioconda** preferred; or release ZIPs via the resolver | DIAMOND is **GPL-3.0** — ship its license text in the ZIP; keep it as a separate asset, never bundled into the MIT wheel. |
| **KEGG HMMs / tables** | **Zenodo** (DOI, >2 GB, archival) | ⚠️ Derived from the subscription-licensed KEGG dump — **confirm redistribution rights with KEGG before publishing publicly**. If not permitted, keep access-gated and have users build from their own dump (the resolver supports a local dir). |
| **Template models** (Human-GEM, yeast-GEM) | **Don't re-host** | Fetch from their canonical repos by pinned release tag — respects their licenses and avoids stale copies. |
