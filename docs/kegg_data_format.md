# KEGG relational-table storage format

This note records *why* ravengem stores its KEGG-derived relational tables as
**gzipped TSV**, and what other options we deliberately deferred. It applies to
the maintainer-built KEGG artefacts described in PLAN.md §2.3b — the `ko_reaction`,
`organism_gene_ko`, `phyl_dist`, KO-name, and reaction-flag tables. (The reference
GEM itself is a cobra model file — YAML/SBML — not covered here.)

## Decision (current)

**Gzipped TSV (`.tsv.gz`)**, partitioned per organism for the large
`organism_gene_ko` table.

- **pandas reads/writes it with zero extra dependencies** — `pd.read_csv` /
  `DataFrame.to_csv` with `sep="\t", compression="gzip"` are built in.
- **MATLAB reads it natively** — `readtable` handles TSV with no toolbox.
- This makes it the genuinely *dependency-free, cross-language* format, which is
  exactly the requirement: the **same files** must serve both the Python
  (ravengem) and MATLAB (RAVEN) sides.

The tables are small by design (minimal columns, gene-free reference GEM) and are
read once per reconstruction, so TSV's parsing/size overhead is not a practical
concern at our current scale.

## Options considered

| Format | Python cost | MATLAB cost | Notes |
| --- | --- | --- | --- |
| **Gzipped TSV** ✅ | none (stdlib/pandas) | none (`readtable`) | Universal, text, types re-specified on read. Chosen. |
| Parquet | `pyarrow` or `fastparquet` (~40–60 MB wheel) as a `ravengem[kegg]` extra | needs ≥ R2019a (`parquetread`, native) | Smaller, faster, typed, columnar. Win mainly at scale / repeated random access. |
| SQLite | none (stdlib `sqlite3`) | **needs Database Toolbox** | Rejected: the MATLAB-side toolbox requirement breaks the "same files, both languages, no extra deps" goal. |

## When to revisit

Reconsider Parquet (or SQLite) if any of these become true:

- The `organism_gene_ko` table grows large enough that gzipped-TSV load time or
  on-disk size becomes a real bottleneck in reconstruction.
- We start doing repeated random-access / columnar reads rather than a single
  load-once-per-run pattern.
- A typed, self-describing schema becomes valuable (TSV loses dtypes; they are
  re-specified on read).

If revisited, prefer **Parquet** over SQLite (no MATLAB toolbox dependency; MATLAB
reads Parquet natively from R2019a). It could be offered as an optional
`ravengem[kegg]` extra (pyarrow) alongside the TSV default, rather than replacing
it — keeping the dependency-free path intact for users who don't opt in.
