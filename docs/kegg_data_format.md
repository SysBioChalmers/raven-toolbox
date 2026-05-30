# KEGG relational-table storage format

This note records *why* raven-python stores its KEGG-derived relational tables as
**gzipped TSV**, and what other options we deliberately deferred. It applies to
the maintainer-built KEGG artefacts described in PLAN.md §2.3b — the `ko_reaction`,
`organism_gene_ko`, KO-name, and reaction-flag tables.

The reference GEM itself is stored as **gzipped RAVEN/cobra YAML**
(`reference_model.yml.gz`) — RAVEN-native and MATLAB-readable, gzipped to match the
tables (the YAML I/O transparently gzips on a `.gz` suffix). On the real KEGG dump
this is ~1.1 MB (vs ~30 MB as SBML) for the full 12k-reaction gene-free model.

End users do not build any of this: the published artefacts are fetched and cached
under `~/.cache/raven-python/data/kegg-<version>/` by `ensure_data` (see
`raven_python.data`), mirroring how binaries are provisioned.

## Decision (current)

- **Small tables** (`ko_reaction`, `ko_names`, `rxn_flags`): **gzipped TSV
  (`.tsv.gz`)**. Each is well under 1 MB, so compression choice is irrelevant;
  gzip keeps them MATLAB-native and dependency-free.
- **The large `organism_gene_ko` table**: **xz-compressed TSV
  (`organism_gene_ko.tsv.xz`), with rows sorted by `(organism, gene)`**.

Why the large table differs. It carries KEGG's ~9M gene↔KO associations and
dominates the artefact set (≈78 MB as unsorted gzipped TSV). Two cheap,
stdlib-only changes cut that to ≈27 MB (2.9×):

1. **Sort by `(organism, gene)`** before writing. Gene IDs from one organism
   share long common prefixes (locus tags, numeric runs); sorting makes them
   adjacent so the compressor can fold them. This alone takes 78 → 48 MB and
   happens to match the by-organism query pattern in
   `get_kegg_model_for_organism`. The sort is an external merge sort bounded to
   `chunk_rows` in memory (see `stream_organism_gene_ko`), so it stays scalable.
2. **xz instead of gzip** (Python stdlib `lzma`). Its larger dictionary captures
   cross-row redundancy gzip's 32 KB window misses: sorted + xz reaches ≈27 MB.

- **pandas reads/writes both with zero extra dependencies** — compression is
  inferred from the `.gz`/`.xz` suffix; `lzma` and `gzip` are both stdlib, so
  this works natively on Windows, macOS, and Linux with no external binary.
- **MATLAB caveat:** `readtable` reads gzipped TSV after a `gunzip`, but MATLAB
  has no built-in xz decompressor. The small tables stay MATLAB-native; the
  large table needs an external `unxz` (or Java/`7-Zip`) before `readtable` on
  the MATLAB side. The xz file is raven-python's (Python) primary artefact; this
  trades a little MATLAB convenience on the one big file for a ~3× size cut.

## Options considered

| Format | Python cost | MATLAB cost | Notes |
| --- | --- | --- | --- |
| **Gzipped TSV** ✅ | none (stdlib/pandas) | none (`readtable`) | Universal, text, types re-specified on read. Chosen. |
| Parquet | `pyarrow` or `fastparquet` (~40–60 MB wheel) as a `raven-python[kegg]` extra | needs ≥ R2019a (`parquetread`, native) | Smaller, faster, typed, columnar. Win mainly at scale / repeated random access. |
| SQLite | none (stdlib `sqlite3`) | **needs Database Toolbox** | Rejected: the MATLAB-side toolbox requirement breaks the "same files, both languages, no extra deps" goal. |

## When to revisit

Reconsider Parquet (or SQLite) if any of these become true:

- The `organism_gene_ko` table grows large enough that load *time* (not just
  size — the sort+xz change above already addresses on-disk size) becomes a real
  bottleneck. The remaining inefficiency is that building one species' model
  still loads all ~9M rows; sorted order makes a `searchsorted`/row-group
  by-organism read the natural next step before reaching for Parquet.
- We start doing repeated random-access / columnar reads rather than a single
  load-once-per-run pattern.
- A typed, self-describing schema becomes valuable (TSV loses dtypes; they are
  re-specified on read).

If revisited, prefer **Parquet** over SQLite (no MATLAB toolbox dependency; MATLAB
reads Parquet natively from R2019a). It could be offered as an optional
`raven-python[kegg]` extra (pyarrow) alongside the TSV default, rather than replacing
it — keeping the dependency-free path intact for users who don't opt in.
