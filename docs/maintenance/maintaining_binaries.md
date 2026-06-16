# Maintaining bundled binaries (BLAST+, DIAMOND, …)

Audience: **raven-toolbox maintainers / the GitHub repo owner.** This explains how
raven-toolbox ships external command-line tools, how to update their versions, and how
to build **minimal-footprint** ZIPs to attach to a GitHub release.

> End users never read this. They get a binary automatically via `ensure_binary`,
> or use their own (system/conda) install. This doc is only for whoever publishes
> the release assets.

---

## 1. How binary provisioning works

raven-toolbox does **not** vendor binaries in the git repo or on PyPI. Instead:

1. For each tool we publish **version-pinned ZIPs as GitHub release assets**.
2. A **registry** (`src/raven_toolbox/binaries_registry.json`) maps each *bundle* to its
   version, the executables it provides, and per-platform `{asset, sha256}`.
3. At run time `raven_toolbox.binaries.ensure_binary("blastp")` resolves a tool in this
   order — and only reaches the download as a last resort:

   ```
   explicit binary= arg  →  env var (RAVEN_PYTHON_BLASTP / RAVEN_PYTHON_DIAMOND / …)
     →  shutil.which on PATH (system / conda / apt / brew)
     →  ensure_binary: download the pinned ZIP → verify SHA256 → cache → return path
     →  actionable error (with conda / manual instructions)
   ```

So a pre-installed binary always wins; the bundle is the zero-setup fallback.
Pinning the version makes reconstruction **reproducible**.

A *bundle* can provide several executables from one download (e.g. the `blast`
bundle provides both `blastp` and `makeblastdb`), so they are fetched once.

---

## 2. What raven-toolbox actually needs — ship only these

Distribute the **minimum** set of executables. Everything else (other suite
tools, docs, examples, changelogs) must be excluded.

| Bundle | Executables to include | Everything else |
|---|---|---|
| `diamond` | `diamond` | — (it is a single static binary) |
| `blast` | `blastp`, `makeblastdb` | **drop** `blastn`, `tblastn`, `psiblast`, `rpsblast`, `blast_formatter`, `*_vdb`, the `doc/`, `ChangeLog`, `README`, ~30 other tools |

(Confirmed against RAVEN `getBlast`/`getDiamond`: only `makeblastdb`+`blastp`, and
`diamond` for its `makedb`/`blastp` subcommands, are ever invoked.)

For BLAST+ this is the big win: the full NCBI suite is ~hundreds of MB; two
binaries (stripped) are a small fraction.

---

## 3. Asset & ZIP conventions

**Asset filename:** `<bundle>-<version>-<os>-<arch>.zip`

- `<os>` ∈ `linux`, `macos`, `windows`
- `<arch>` ∈ `x86_64`, `arm64`
- examples: `diamond-2.1.11-linux-x86_64.zip`, `blast-2.16.0-macos-arm64.zip`

**ZIP layout — flat, executables at the root, plus the upstream licence:**

```
diamond-2.1.11-linux-x86_64.zip
├── diamond
└── LICENSE

blast-2.16.0-linux-x86_64.zip
├── blastp
├── makeblastdb
└── LICENSE
```

No nested `bin/`, no extra files. `ensure_binary` extracts the ZIP into the cache
and expects the executable at the top level.

**Windows ZIPs** ship the `.exe` (and any runtime DLL the build needs) at the root:

```
blast-2.17.0-windows-x86_64.zip       hmmer-3.3.2-windows-x86_64.zip
├── blastp.exe                        ├── hmmsearch.exe
├── makeblastdb.exe                   ├── cygwin1.dll        # Cygwin runtime, required
├── nghttp2.dll        # required     └── LICENSE
└── LICENSE
```

`ensure_binary` looks up the executable as `<name>.exe` on Windows (falling back to
the bare name), so the resolver returns `hmmsearch.exe` and the sibling DLL sits
next to it in the cache dir where the OS loader finds it.

---

## 4. Step-by-step: add or update a version

Example: bump DIAMOND to a new version for Linux x86-64. Repeat per `(os, arch)`.

1. **Download the official upstream build** (never rebuild from source unless you
   must):
   - DIAMOND → <https://github.com/bbuchfink/diamond/releases>
     (`diamond-linux64.tar.gz`, `diamond-macos.tar.gz`)
   - BLAST+ → <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/> or a
     pinned version dir (`ncbi-blast-<ver>+-x64-linux.tar.gz`,
     `-x64-macosx.tar.gz`, `-aarch64-linux.tar.gz`, `-x64-win64.tar.gz`).
   - Record the upstream URL **and** its published checksum for provenance.
2. **Extract only the needed executables** (see §2) to a clean staging dir.
3. **Strip debug symbols** to shrink (skip on Windows / signed macOS builds):
   ```bash
   strip diamond           # or: strip blastp makeblastdb
   ```
4. **Smoke-test the stripped binaries in a clean shell** (no other tools on PATH):
   ```bash
   ./diamond --version
   ./blastp -version && ./makeblastdb -version
   ```
   If they fail for a missing shared library, add that `.so`/`.dylib` to the ZIP
   (rare — NCBI/DIAMOND release builds are largely self-contained).
5. **Add the upstream licence file** as `LICENSE` (see §6).
6. **Zip with max compression, flat layout:**
   ```bash
   zip -9 -j diamond-2.1.11-linux-x86_64.zip diamond LICENSE
   # -j junks paths so entries sit at the ZIP root
   ```
7. **Compute the SHA256:**
   ```bash
   sha256sum diamond-2.1.11-linux-x86_64.zip   # shasum -a 256 on macOS
   ```
8. **Attach the ZIP to a raven-toolbox GitHub release** (a release tagged for the binary
   set, e.g. `binaries-2024.06`, keeps them independent of code releases).
9. **Update the registry** `src/raven_toolbox/binaries_registry.json` — bump `version`
   and set the per-platform `asset` + `sha256`:
   ```json
   {
     "diamond": {
       "version": "2.1.11",
       "provides": ["diamond"],
       "platforms": {
         "linux-x86_64": {
           "asset": "diamond-2.1.11-linux-x86_64.zip",
           "url": "https://github.com/SysBioChalmers/raven-toolbox/releases/download/binaries-2024.06/diamond-2.1.11-linux-x86_64.zip",
           "sha256": "<sha256>"
         }
       }
     },
     "blast": {
       "version": "2.16.0",
       "provides": ["blastp", "makeblastdb"],
       "platforms": { "linux-x86_64": { "asset": "...", "url": "...", "sha256": "..." } }
     }
   }
   ```
10. **Commit the registry change**, run the homology tests, and (if you have the
    binary) confirm `ensure_binary("diamond", version="2.1.11")` downloads,
    verifies, and runs.

---

## 5. Keeping the footprint minimal — checklist

- ✅ Only the executables in §2 (for BLAST+, exactly `blastp` + `makeblastdb`).
- ✅ `strip` the binaries (often halves their size).
- ✅ `zip -9 -j` (max compression, flat — no `bin/`, no folders).
- ✅ Exactly one extra file: `LICENSE`.
- ❌ No docs, examples, `ChangeLog`, `README`, man pages, test data, or sibling tools.
- ❌ No `.dSYM`/debug bundles; no duplicate static `.a` libraries.
- ➕ Only add a shared library if step-4 testing proves it is required.

---

## 6. Platform / architecture matrix & licensing

**Coverage = what you build.** Start with `linux-x86_64` (CI default), then add
`macos-arm64`, `macos-x86_64`, `linux-arm64`, `windows-x86_64` as capacity allows.
For any `(os, arch)` **not** in the registry, `ensure_binary` raises an actionable
error pointing to conda (`conda install -c bioconda diamond blast`) or a manual
install — that is the documented fallback, not a failure to fix urgently.

**Licensing (must comply when redistributing):**

- **BLAST+** — produced by NCBI (US Government); **public domain**, free to
  redistribute. Include NCBI's `LICENSE` for courtesy/provenance.
- **DIAMOND** — **GPLv3**. Redistribution is allowed; you **must** include the
  GPLv3 licence text in the ZIP and keep the binary unmodified (or offer source).
- **HMMER** (future) — BSD-3-Clause; include its `LICENSE`.

Always ship the upstream licence in the ZIP, and keep a `BINARIES_PROVENANCE.md`
(or a note in the release body) recording, per asset: upstream URL, upstream
version, upstream checksum, and the SHA256 you published.

### Native OS support per tool

raven-toolbox invokes each tool through `subprocess.run([resolved_path, …])` — that
call is itself cross-platform, so the real constraint is whether a given tool has
a binary that runs natively on each OS. It varies:

| Tool | Linux | macOS (incl. arm64) | Windows (native) |
|---|---|---|---|
| BLAST+ (`blastp`, `makeblastdb`) | ✅ | ✅ | ✅ (NCBI ships Windows builds; needs `nghttp2.dll`) |
| DIAMOND | ✅ | ✅ | ✅ (RAVEN ships `diamond.exe`) |
| HMMER `hmmsearch` (query, 3b.5) | ✅ 3.4 | ✅ 3.4 | ✅ **3.3.2** (Cygwin build from RAVEN 2.10.5; needs `cygwin1.dll`) |
| HMMER `hmmbuild` (build, 3b.3) | ✅ 3.4 | ✅ 3.4 | ⚠️ 3.3.2 `.exe` exists (RAVEN 2.10.5) but the build also needs MAFFT/CD-HIT |
| MAFFT | ✅ | ✅ | ❌ no usable native build |
| CD-HIT | ✅ | ✅ | ❌ no Windows build exists |

(Sources for the bundled builds and versions: RAVEN `develop3` `software/` — BLAST+
2.17.0, DIAMOND 2.1.17, HMMER 3.4.0 hmmsearch for Linux/macOS — and RAVEN `v2.10.5`,
which is the last release to ship native-Windows HMMER 3.3.2 `hmmsearch.exe`/`hmmbuild.exe`.)

Implications:

- **Linux / macOS** — everything works. `conda install -c bioconda hmmer mafft
  cd-hit blast diamond`, or point the `RAVEN_PYTHON_*` env vars at your installs.
- **Native Windows — runtime (end users) works.** BLAST+, DIAMOND, and `hmmsearch`
  all have native Windows builds, so homology reconstruction **and** the KEGG HMM
  *query* (3b.5) run without WSL. The `hmmsearch` we bundle for Windows is HMMER
  **3.3.2** (from RAVEN 2.10.5); it reads the same `HMMER3/f` profile format that
  3.4's `hmmbuild` writes, so it searches the published 3.4-built HMM libraries
  identically. See §10.
- **Native Windows — the HMM *build* (3b.3) does not work.** Even though a 3.3.2
  `hmmbuild.exe` exists, the build pipeline also needs **MAFFT and CD-HIT**, which
  have no Windows binaries (and no bioconda Windows packages). Bundling can't fix
  that — there is nothing to bundle. Build on **WSL2** (or Linux/macOS).
- raven-toolbox does **not** replicate RAVEN's `getWSLpath`/`wsl …` path
  translation: it calls the resolved binary directly, so mixing native-Windows
  Python with WSL binaries is unsupported — for the build, keep the whole stack
  inside WSL2.
- The common end-user paths — homology reconstruction and the KEGG *species* model
  (3b.4) — need no HMMER/MAFFT/CD-HIT, so they are fully cross-platform.

---

## 7. Emitting the registry entry

After building the per-platform ZIPs (named `<bundle>-<version>-<os>-<arch>.zip`)
and uploading them to the release, generate the `_REGISTRY` entry — checksums and
URLs — with [`scripts/make_registry_snippet.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/scripts/README.md):

```bash
python scripts/make_registry_snippet.py binary --bundle blast --version 2.16.0 \
    --provides blastp makeblastdb --dir zips \
    --base-url https://github.com/ORG/raven-toolbox/releases/download/blast-2.16.0
```

It prints the ready-to-paste `_REGISTRY["blast"]` block; its SHA256 helper is the
same one `ensure_binary` verifies with, so the checksums always match. (Producing
the minimal ZIPs themselves — download upstream, `strip`, `zip -9 -j`, add
`LICENSE` per §3–§6 — is still a manual/per-tool step.)

---

## 8. Adding a new tool later (e.g. HMMER for KEGG reconstruction)

1. Decide the **minimal executable set** (e.g. HMMER → `hmmsearch`, `hmmscan`,
   maybe `hmmbuild`/`hmmpress`).
2. Add a bundle entry to the registry with `provides` listing those executables.
3. Build/attach ZIPs per §3–§4; include the tool's licence (§6).
4. The wrappers call `ensure_binary("hmmsearch", …)` with the same resolution
   order — no new provisioning code needed.

---

## 9. Binary sets, the CLI, and auto-fetch control

End users don't need every tool. raven-toolbox groups the executables into two
**sets** (`raven_toolbox.binaries.BINARY_SETS`) for the two audiences:

| Set | Executables | Audience |
|---|---|---|
| `runtime` | `blastp`, `makeblastdb`, `diamond`, `hmmsearch` | end users — homology + KEGG HMM query (3b.5) |
| `build` | `hmmbuild`, `mafft`, `cd-hit` | maintainers — KEGG HMM-library build (3b.3) |

**Why not a `pip install raven-toolbox[runtime]` extra?** pip extras can only pull
**PyPI wheels**, and none of these tools (except HMMER, via `pyhmmer`) are on PyPI;
and downloading binaries during `pip install` only runs for sdists, breaks
offline/locked-down installs, and is a known anti-pattern. So provisioning is
**decoupled from pip** into two layers:

1. **Explicit fetch (console script).** After `pip install`, run:
   ```bash
   raven-toolbox-binaries --set runtime   # blast + diamond + hmmsearch for this OS
   raven-toolbox-binaries --set build     # hmmbuild + mafft + cd-hit
   raven-toolbox-binaries --list          # show the sets for this platform
   ```
   It fetches via `ensure_binary` (SHA256-verified, cached), **skips tools already
   on PATH**, and reports any with no bundle for this OS/arch (with a conda/WSL2
   hint) instead of failing.
2. **Lazy first-use download.** Any wrapper that needs a tool calls
   `resolve_binary`, which downloads the bundle on first use if it isn't already
   resolvable. This is the zero-setup default.

**Turning auto-fetch off.** Set `RAVEN_PYTHON_AUTOFETCH=0` (also `false`/`no`/`off`)
to stop `resolve_binary` ever reaching the network: resolution then ends at
arg → env var → PATH and raises an actionable error otherwise. For air-gapped or
strictly conda/system-managed environments. (The `raven-toolbox-binaries` command
still fetches when run explicitly.)

---

## 10. Native-Windows HMMER (3.3.2 from RAVEN 2.10.5)

There is **no native-Windows HMMER 3.4** build — the project targets POSIX, and
RAVEN's `develop3` ships only Linux/macOS `hmmsearch` 3.4.0 (it runs the Linux
binary via WSL on Windows). But RAVEN **`v2.10.5`** is the last release to bundle a
**native Windows HMMER 3.3.2**, Cygwin-compiled: `hmmsearch.exe`, `hmmbuild.exe`,
and the `cygwin1.dll` they depend on.

So a native-Windows end user has two choices for the KEGG HMM *query* (3b.5):

- **Native HMMER 3.3.2** — register a `windows-x86_64` entry on the `hmmer` bundle
  pointing at a `hmmer-3.3.2-windows-x86_64.zip` (`hmmsearch.exe` + `cygwin1.dll` +
  `LICENSE`, repackaged from RAVEN 2.10.5). `ensure_binary` then provisions it
  automatically — no WSL.
- **HMMER 3.4 via WSL2** — run the whole stack inside WSL2 (Linux Python + Linux
  `hmmsearch` 3.4).

**Is searching 3.4-built HMMs with 3.3.2 a problem? No.** raven-toolbox publishes a
concatenated **ASCII** `.hmm` library and searches it directly (no `hmmpress`, so no
version-sensitive binary `.h3m/.h3f/.h3i/.h3p`). The ASCII profile format is
`HMMER3/f`, unchanged from 3.1 through 3.4, and `hmmbuild` 3.4 writes it; `hmmsearch`
3.3.2 reads it. HMMER 3.4 is a maintenance release over 3.3.2 with no change to the
protein scoring model, so bit scores / E-values match and the calibrated KEGG
cutoffs transfer. Caveats: ship `cygwin1.dll` with the `.exe`; it's an older,
unmaintained build; and **keep publishing ASCII libraries** (don't switch to
`hmmpress`-ed binaries) to preserve this cross-version compatibility. A one-time
check — search a fixed test set with 3.3.2 (Windows) and 3.4 (Linux) on the same
library and confirm identical hits — is cheap insurance worth recording in the
[KEGG HMM cutoff calibration study](../studies/kegg_hmm_cutoff_calibration.md).

> The example registry entries (including this Windows 3.3.2 asset) are in
> [`data/manifest.example.json`](../../data/manifest.example.json). The bundle's
> single `version` is the Linux/macOS version (3.4.0); the Windows asset URL names
> 3.3.2. If you ever need fully per-platform versions, add an optional `version` to
> the manifest's per-platform `file` entries and thread it into the
> `ensure_binary` cache path — a small schema change noted for the future.
