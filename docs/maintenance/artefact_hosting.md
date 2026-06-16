# Artefact hosting & publishing (the raven-data repo)

Audience: **raven-toolbox / RAVEN maintainers.** This is the one-stop explanation of
where the large downloadable artefacts live, how they're versioned, and the exact
commands to publish a new one. The per-artefact build details are in
[Maintaining KEGG data](maintaining_kegg_data.md) and
[Maintaining binaries](maintaining_binaries.md); the file format is in
[Data & binary manifest](data_manifest.md).

## Where everything lives

All large downloadables — KEGG reference data, prebuilt HMM libraries, and the
BLAST+/DIAMOND/HMMER binaries — are hosted as **GitHub Release assets in a dedicated
repository, [`SysBioChalmers/raven-data`](https://github.com/SysBioChalmers/raven-data)**,
indexed by a single `manifest.json`. They are **not** attached to the
`raven-toolbox` or `RAVEN` code releases.

Why a separate repo:

- The HMM libraries are ~135–155 MB each — past GitHub's **100 MB per-file git
  limit** — so they can only be release assets, never committed into a repo's tree.
- It keeps the code repos' release pages and licenses clean (KEGG terms, GPL/PD/BSD
  binaries stay on the data repo).
- One shared, versioned source of truth serves **both** raven-toolbox (Python) and
  MATLAB RAVEN.

## Versioning model (no duplication)

There is **no single "raven-data version".** Each artefact has its own **immutable
release tag**, versioned by its *upstream* version, holding only that artefact's
assets:

| Tag | Assets |
|---|---|
| `kegg118` (`kegg119`, …) | `kegg<NNN>_core.tar.gz`, `_taxonomy.gz`, `_prokaryotes.hmm.gz`, `_eukaryotes.hmm.gz` |
| `blast-2.17.0` | `blast-2.17.0-{linux-x86_64,macos-arm64,windows-x86_64}.zip` |
| `diamond-2.1.17` | the three DIAMOND ZIPs |
| `hmmer-3.4.0` | Linux + macOS `hmmsearch` ZIPs |
| `hmmer-3.3.2` | the native-Windows `hmmsearch` ZIP (Cygwin) |
| `manifest-v1` (`-v2`, …) | `manifest.json` (a pinned snapshot) |

An asset is uploaded **once** under its tag and **never re-uploaded**. Bumping one
tool (e.g. DIAMOND → 2.1.18) means a *new* tag `diamond-2.1.18` with only its ZIPs;
every other tag is untouched. The "snapshot of the whole set" is the tiny
`manifest.json` (versioned `manifest-v1`, `-v2`, …), not a release that re-bundles
binaries — so binaries are stored once and merely re-referenced. (`manifest_version`
*inside* the JSON is the **schema** version, independent of the snapshot tag.)

## How consumers use it

Both tools read the same `manifest.json` and **verify every file's SHA256** after
download:

- **raven-toolbox** bakes a pinned snapshot of the manifest into its default
  registries (`raven_toolbox.data._DATA_REGISTRY`, `raven_toolbox.binaries._REGISTRY`),
  so a given code release always fetches the exact, checksum-verified assets it was
  tested against. `RAVEN_PYTHON_MANIFEST=<url|path>` overrides with a newer manifest.
- **RAVEN (MATLAB)** `getKEGGModelForOrganism` downloads the HMM libraries from the
  matching `kegg<NNN>` release.

The tools coordinate **through the manifest + SHA256 pinning**, not by matching
version numbers — each chooses which manifest snapshot it targets.

## Publishing a new artefact — the whole workflow

Everything is scripted and idempotent. Re-running skips assets already present, so a
single-tool bump never re-uploads the rest.

```bash
# 0. one-time: gh auth with write access to SysBioChalmers/raven-data

# 1. BUILD the assets (from raven-toolbox)
python scripts/build_binary_bundles.py                 # -> dist/binaries/*.zip (+ checksums, PROVENANCE)
python scripts/build_kegg_artefacts.py --keggdb … --out … --version kegg118 --hmms

# 2. PUBLISH to raven-data releases (idempotent; immutable tags)
python scripts/publish_to_raven_data.py binaries --dir dist/binaries
python scripts/publish_to_raven_data.py release --tag kegg118 --dir <kegg-out-dir>

# 3. UPDATE the manifest (compute URLs + SHA256 + bytes from the uploaded files)
python scripts/make_registry_snippet.py manifest --target data --dataset kegg \
    --version kegg118 --dir <kegg-out-dir> --manifest data/manifest.json \
    --base-url https://github.com/SysBioChalmers/raven-data/releases/download/kegg118
#   (binaries: one `--target binary` call per bundle, or hand-merge the
#    dist/binaries/manifest_binaries.json that build_binary_bundles.py emits —
#    it already handles HMMER's split 3.4.0/3.3.2 versions.)

# 4. SYNC the baked Python registries from the manifest (single source of truth)
python scripts/make_registry_snippet.py sync           # rewrites _DATA_REGISTRY + _REGISTRY

# 5. PUBLISH the manifest snapshot, then open the raven-toolbox PR
python scripts/publish_to_raven_data.py release --tag manifest-v2 data/manifest.json
gh pr create -t "Adopt raven-data manifest-v2 (kegg118)" -b "…"
```

`--dry-run` on `publish_to_raven_data.py` prints the `gh` calls without running them.

### What stays manual
- Building the **KEGG dump** (needs `~/.netrc` FTP credentials + a multi-GB download
  and the KEGG license) — run by a maintainer, not CI.
- Deciding **whether** to adopt a new upstream tool / KEGG release into a code
  release — a reviewed version bump; the scripts just execute it.

## Editing rule
`data/manifest.json` is the **single source of truth**. Never hand-edit the baked
registries in `data.py`/`binaries.py` — change the manifest and run
`make_registry_snippet.py sync`. Per-asset upstream sources and checksums are
recorded in raven-data's `PROVENANCE-binaries.txt`.
