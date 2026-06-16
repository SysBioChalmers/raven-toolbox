# Binary provisioning plan (BLAST / DIAMOND / HMMER / MAFFT / CD-HIT)

Status: **partly implemented** (code + docs landed; release ZIPs + real manifest
URLs/SHA256 remain a maintainer step). Audience: raven-toolbox maintainers.

This records the design for *how users choose and obtain* the external command-line
tools raven-toolbox shells out to, including native-Windows HMMER. Operational
how-to (building ZIPs, the registry) lives in
[`maintaining_binaries.md`](../maintenance/maintaining_binaries.md); this doc is the
rationale and the survey behind it.

## 1. Problem

Two audiences want different tools, and OS coverage differs per tool:

- **End users (runtime):** `blastp`, `makeblastdb`, `diamond` (homology) and
  `hmmsearch` (KEGG HMM query, step 3b.5).
- **Developers (build):** `hmmbuild`, `mafft`, `cd-hit` (KEGG HMM-library build,
  step 3b.3).

Question raised: *can the choice — and whether binaries auto-download — be made at
`pip install` time?*

## 2. What RAVEN actually ships (surveyed)

| Tool | Linux | macOS (`.mac`) | Windows (`.exe`) | Version | RAVEN source |
|---|---|---|---|---|---|
| BLAST+ (`blastp`,`makeblastdb`) | ✅ | ✅ | ✅ (+`nghttp2.dll`) | 2.17.0 | `develop3` `software/` |
| DIAMOND | ✅ | ✅ | ✅ | 2.1.17 | `develop3` `software/` |
| HMMER `hmmsearch` | ✅ | ✅ | ❌ | 3.4.0 | `develop3` `software/hmmer` |
| HMMER `hmmsearch.exe` + `hmmbuild.exe` | ✅ | ✅ | ✅ (Cygwin, +`cygwin1.dll`) | **3.3.2** | `v2.10.5` `software/hmmer` |
| MAFFT, CD-HIT | — | — | — | — | not bundled by RAVEN |

Key facts:

- RAVEN's `develop3` bundles only the **runtime/query** tools (blast, diamond,
  hmmsearch) — not the build tools (mafft, cd-hit, hmmbuild). This matches the
  users-vs-developers split exactly.
- There is **no native-Windows HMMER 3.4** anywhere; `develop3` drops the Windows
  HMMER and uses WSL. But RAVEN **`v2.10.5`** is the last release with a native
  Windows HMMER **3.3.2** (Cygwin) — `hmmsearch.exe` and `hmmbuild.exe`.
- DIAMOND's Windows build is real (`diamond.exe`), correcting the earlier doc that
  called it "Linux-first".

## 3. Why not `pip install raven-toolbox[runtime]`

pip extras can only depend on **PyPI wheels**. Of these tools only HMMER has one
(`pyhmmer`, Linux/macOS only). Downloading binaries during `pip install` runs only
for sdists (not wheels), breaks offline/isolated/locked-down installs, and is a
known anti-pattern. So provisioning is **decoupled from pip**.

## 4. Design (implemented)

Three layers, plus a Windows-specific fix:

1. **Binary sets** — `binaries.BINARY_SETS = {runtime, build}`; `all` is the union.
2. **Explicit fetch CLI** — `raven-toolbox-binaries --set runtime|build|all`
   (console script → `binaries_cli:main`). Fetches via `ensure_binary`, skips tools
   already on PATH, reports per-tool `present/downloaded/unavailable/error`, exits
   non-zero only on a real download error (not on an OS that simply has no bundle).
3. **Lazy first-use download + master switch** — `resolve_binary` still downloads on
   first use (zero-setup default); `RAVEN_PYTHON_AUTOFETCH=0` disables that so
   resolution stops at arg → env var → PATH. (Decision: **on by default**, opt-out.)
4. **Windows `.exe` resolution** — `ensure_binary` now looks up `<name>.exe` on
   Windows (falling back to the bare name) so Windows bundles resolve at all; sibling
   DLLs (`nghttp2.dll`, `cygwin1.dll`) extract next to the `.exe`.

Sourcing: reuse RAVEN's vetted binaries, repackaged into raven-toolbox release ZIPs
(`<bundle>-<ver>-<os>-<arch>.zip`) and registered in `data/manifest.json` — **not**
vendored in-repo (keeps the wheel lean; consistent with the existing model). The
build tools (mafft, cd-hit) aren't in RAVEN; repackage from upstream for Linux/macOS
or point users to conda.

## 5. Native-Windows HMMER: 3.3.2 vs 3.4-via-WSL

Offer both. Register a `windows-x86_64` entry on the `hmmer` bundle pointing at a
`hmmer-3.3.2-windows-x86_64.zip` (repackaged from RAVEN 2.10.5: `hmmsearch.exe` +
`cygwin1.dll` + `LICENSE`); `ensure_binary` then auto-provisions it on native
Windows. Users who prefer 3.4 run inside WSL2.

**Cons of searching 3.4-built HMMs with 3.3.2: effectively none.** raven-toolbox
publishes a concatenated **ASCII** `.hmm` library and searches it directly (no
`hmmpress`, so no version-sensitive binary index). The ASCII format is `HMMER3/f`,
unchanged 3.1→3.4; 3.4's `hmmbuild` writes it and 3.3.2's `hmmsearch` reads it. 3.4
is a maintenance release over 3.3.2 with no protein-scoring change, so scores/E-values
match and cutoffs transfer. Caveats: ship `cygwin1.dll`; it's an older unmaintained
build; keep libraries ASCII; validate once (search a test set with both versions).
Full write-up in [`maintaining_binaries.md` §10](../maintenance/maintaining_binaries.md).

## 6. What remains (maintainer steps)

- Build/upload the per-platform ZIPs (incl. `hmmer-3.3.2-windows-x86_64.zip`) and
  fill the **real** `data/manifest.json` `binaries` block with URLs + SHA256 (use
  `scripts/make_registry_snippet.py binary`). Example entries (zeroed SHA) are in
  `data/manifest.example.json`.
- Confirm the actual arch of RAVEN's `.mac` HMMER/BLAST/DIAMOND builds (x86_64 vs
  universal) before publishing the `macos-*` keys.

## 7. Possible future enhancements (not done)

- **`pyhmmer` extra** — `raven-toolbox[hmmer]` → `pyhmmer` gives HMMER via a Python
  API on Linux/macOS with no binary download (would need an optional pyhmmer-backed
  `hmmsearch` path; doesn't help Windows).
- **Per-platform versions** — add an optional `version` to the manifest per-platform
  `file` entries and thread it into the `ensure_binary` cache path, so the `hmmer`
  bundle can carry 3.4.0 (Linux/macOS) and 3.3.2 (Windows) as first-class versions
  rather than via the asset URL.
