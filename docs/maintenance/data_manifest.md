# Data & binary manifest

Large artefacts (KEGG tables / HMMs, template models) and external-binary bundles
(BLAST / DIAMOND / HMMER) are **not** committed to the code repository. They are published
as downloadable assets and described by a single, language-agnostic **manifest** that both
raven-toolbox and MATLAB RAVEN read. Every file carries a **SHA256**, so consumers verify
integrity after download.

- Format: [`data/manifest.schema.json`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/data/manifest.schema.json) (JSON Schema)
- Worked example: [`data/manifest.example.json`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/data/manifest.example.json)
- Live manifest: [`data/manifest.json`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/data/manifest.json) (empty until assets are published)

The manifest is a superset of the two runtime registries:

| Manifest section | Runtime registry |
| --- | --- |
| `data` | {data}`raven_toolbox.data._DATA_REGISTRY` |
| `binaries` | `raven_toolbox.binaries._REGISTRY` |

```json
{
  "manifest_version": 1,
  "data":     { "<dataset>": { "version": "...", "doi": "...", "files":     { "<name>": {"url": "...", "sha256": "...", "bytes": 0} } } },
  "binaries": { "<bundle>":  { "version": "...", "provides": ["..."], "platforms": { "<os>-<arch>": {"url": "...", "sha256": "...", "bytes": 0} } } }
}
```

## Consuming it — Python

Point raven-toolbox at a manifest and the resolvers populate themselves on first use,
verifying each download's checksum:

```bash
export RAVEN_PYTHON_MANIFEST=https://github.com/SysBioChalmers/raven-toolbox/releases/download/manifest-v1/manifest.json
```

```python
from raven_toolbox import manifest
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
([`scripts/make_registry_snippet.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/scripts/make_registry_snippet.py)),
which computes each SHA256 and byte size. Pass `--tag` (the GitHub release tag); the
script builds the `https://github.com/SysBioChalmers/raven-toolbox/releases/download/<tag>`
asset URL itself:

```bash
python scripts/make_registry_snippet.py manifest --manifest data/manifest.json \
    --target data --dataset kegg --version kegg116 --dir artefacts --tag kegg-kegg116 \
    --doi 10.5281/zenodo.0000000 --source https://zenodo.org/records/0000000

python scripts/make_registry_snippet.py manifest --manifest data/manifest.json \
    --target binary --bundle diamond --version 2.1.9 --provides diamond --dir zips \
    --tag diamond-2.1.9 --license GPL-3.0-only
```

## Where to host

Release **assets are stored separately from the git tree** (GitHub keeps them in a blob
store), so attaching them to a release does **not** bloat the repository. A dedicated assets
repository is therefore **optional** — attach the assets to releases on an existing RAVEN
repo (this one, or MATLAB [RAVEN](https://github.com/SysBioChalmers/RAVEN)) and have **both
packages reuse the same release-asset URLs** via this manifest.

Use **dedicated tags** for the assets — e.g. `kegg-kegg116`, `diamond-2.1.9` — rather than
attaching them to code-milestone releases like `v0.1.0a1`. KEGG data updates roughly yearly
while the code changes often; dedicated tags keep the two cadences decoupled while still
living in one repository. The manifest's per-dataset `version` does the rest (it namespaces
the download cache).

Both GitHub Releases and Zenodo are just URLs in the manifest, so consumers don't care —
mix them per file:

- **GitHub Releases** — simplest, free, language-agnostic, up to **~2 GB per file**. The
  default home for the manifest and most assets.
- **Zenodo** — adds a citable **DOI**, long-term archival, and handles files **larger than
  2 GB** (up to 50 GB/record). Use it for individual large HMM libraries or anything you want
  citable; point just that file's `url` at the Zenodo record.

### Auto-publishing to Zenodo from GitHub (only if you need DOIs / >2 GB files)

:::{important}
The **native GitHub↔Zenodo integration** (flip a switch, publish a Release → DOI) archives
the **repository source zipball** at the tag — it does **not** capture files attached to the
Release. So it only works for assets *committed into the repo*, which defeats the purpose for
multi-GB binaries. Use it for a *software* DOI, not for the data assets.
:::

If you do want Zenodo DOIs (or need to host files >2 GB), keep it GitHub-driven with a small
**GitHub Action** that, on release, uploads the assets to Zenodo via its REST API (e.g.
[`zenodraft`](https://github.com/zenodraft/zenodraft)). You cut a normal GitHub Release with
the files attached; the Action mirrors them to Zenodo and mints a new version DOI. Drop this
into whichever repo hosts the asset releases as `.github/workflows/zenodo.yml`:

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
| **KEGG HMMs / tables** | GitHub release (dedicated `kegg-*` tag); Zenodo for libraries >2 GB | Derived from the KEGG dump and **redistributed with permission from KEGG**. Note the provenance in the release notes / manifest `license`. |
| **Template models** (Human-GEM, yeast-GEM) | **Don't re-host** | Fetch from their canonical repos by pinned release tag — respects their licenses and avoids stale copies. |
