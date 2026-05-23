# Maintaining bundled binaries (BLAST+, DIAMOND, …)

Audience: **ravengem maintainers / the GitHub repo owner.** This explains how
ravengem ships external command-line tools, how to update their versions, and how
to build **minimal-footprint** ZIPs to attach to a GitHub release.

> End users never read this. They get a binary automatically via `ensure_binary`,
> or use their own (system/conda) install. This doc is only for whoever publishes
> the release assets.

---

## 1. How binary provisioning works

ravengem does **not** vendor binaries in the git repo or on PyPI. Instead:

1. For each tool we publish **version-pinned ZIPs as GitHub release assets**.
2. A **registry** (`src/ravengem/binaries_registry.json`) maps each *bundle* to its
   version, the executables it provides, and per-platform `{asset, sha256}`.
3. At run time `ravengem.binaries.ensure_binary("blastp")` resolves a tool in this
   order — and only reaches the download as a last resort:

   ```
   explicit binary= arg  →  env var (RAVENGEM_BLASTP / RAVENGEM_DIAMOND / …)
     →  shutil.which on PATH (system / conda / apt / brew)
     →  ensure_binary: download the pinned ZIP → verify SHA256 → cache → return path
     →  actionable error (with conda / manual instructions)
   ```

So a pre-installed binary always wins; the bundle is the zero-setup fallback.
Pinning the version makes reconstruction **reproducible**.

A *bundle* can provide several executables from one download (e.g. the `blast`
bundle provides both `blastp` and `makeblastdb`), so they are fetched once.

---

## 2. What ravengem actually needs — ship only these

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
8. **Attach the ZIP to a ravengem GitHub release** (a release tagged for the binary
   set, e.g. `binaries-2024.06`, keeps them independent of code releases).
9. **Update the registry** `src/ravengem/binaries_registry.json` — bump `version`
   and set the per-platform `asset` + `sha256`:
   ```json
   {
     "diamond": {
       "version": "2.1.11",
       "provides": ["diamond"],
       "platforms": {
         "linux-x86_64": {
           "asset": "diamond-2.1.11-linux-x86_64.zip",
           "url": "https://github.com/SysBioChalmers/ravengem/releases/download/binaries-2024.06/diamond-2.1.11-linux-x86_64.zip",
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

---

## 7. Optional: a build helper

A small script (`scripts/build_binary_zip.sh`, to be added) can standardise
steps 2–7: download upstream, extract only the needed executables, `strip`, add
`LICENSE`, `zip -9 -j`, and print the SHA256 + a ready-to-paste registry snippet.
Keeping it in `scripts/` makes version bumps a one-command, reproducible operation.

---

## 8. Adding a new tool later (e.g. HMMER for KEGG reconstruction)

1. Decide the **minimal executable set** (e.g. HMMER → `hmmsearch`, `hmmscan`,
   maybe `hmmbuild`/`hmmpress`).
2. Add a bundle entry to the registry with `provides` listing those executables.
3. Build/attach ZIPs per §3–§4; include the tool's licence (§6).
4. The wrappers call `ensure_binary("hmmsearch", …)` with the same resolution
   order — no new provisioning code needed.
