# Homology reconstruction parameter benchmarks

Functions: `raven_toolbox.reconstruction.homology.blast.run_blast`,
`raven_toolbox.reconstruction.homology.blast.run_diamond`,
`raven_toolbox.reconstruction.homology.homology.get_model_from_homology`

Date: 2026-06-20. Binary: BLAST 2.17.0
(`~/.cache/raven_toolbox/binaries/blast-2.17.0-windows-x86_64/blastp.exe`).
FASTAs: `hanpo-GEM/data/genomes/` — sce.faa (6717 seqs), rhto.faa (8140 seqs),
hanpo.faa (5177 seqs).

---

## `evalue` — BLAST/Diamond search E-value cutoff

**Python default:** `1e-5`
**MATLAB default:** `1e-4` (10e-5)

The E-value is the expected number of random hits with that score or better in a
database of this size. A lower E-value is more stringent (fewer but more confident
hits); a higher E-value admits more hits including more false positives.

**Benchmark (2026-06-20): H. polymorpha vs S. cerevisiae (hanpo.faa vs sce.faa)**

| E-value | Total hits | Time (s) |
|---|---|---|
| `1e-4` (MATLAB) | 58,707 | 1303 |
| `1e-5` (Python) | 53,503 | 1181 |
| Marginal hits (`1e-4` only) | 5,204 (8.9% of 1e-4 hits) | — |

The 5,204 marginal hits are gene pairs with alignment E-values between 1e-5 and 1e-4.
Their identity distribution is pending (follow-up analysis running), but the
`get_model_from_homology` post-BLAST filter applies `min_identity=40` and
`max_evalue=1e-30`, which would discard virtually all of these marginal hits regardless
of the initial BLAST E-value cutoff. The 1e-4 vs 1e-5 distinction matters only for
the raw BLAST table, not the final homology model.

**For closely related organisms** (≥70% AAI, e.g. different *Saccharomyces* species):
the difference between 1e-4 and 1e-5 is negligible — all real homologs score well
above 1e-5.

**For distantly related organisms** (≤30% AAI, e.g. *H. polymorpha* vs bacteria):
neither 1e-4 nor 1e-5 is stringent enough — `get_model_from_homology` applies
additional filters (`max_evalue=1e-30`, `min_identity=40`) that dominate.

**Decision: ✓ keep `evalue=1e-5`.** Matches the `blastp` command-line default.
MATLAB's `1e-4` adds 8.9% more raw hits (5,204 gene pairs) that are filtered out
downstream by `get_model_from_homology`'s identity/evalue cutoffs. No correctness
benefit; marginally more compute.

---

## `threads` — number of parallel BLAST/Diamond processes

**Python default:** `1`
**MATLAB default:** all available cores (auto-detected)

BLAST and Diamond are documented as deterministic across thread counts — the same
hits are returned regardless of parallelism (alignment scores are computed
independently per query sequence). HMMER may show negligible E-value differences
across threads due to floating-point accumulation in background frequency
estimation, but below the threshold of any meaningful cutoff.

**Benchmark (2026-06-20): 500-query subset of hanpo.faa vs sce.faa, threads=1 vs threads=4**

| Threads | Hits | Wall time (s) | Identical results? |
|---|---|---|---|
| 1 | 2,469 | 45.2 | — |
| 4 | 2,469 | 23.7 | ✓ yes |
| Speedup | — | **1.9×** | — |

Hit counts are identical (deterministic). Speedup of 1.9× on this 8-core Windows
machine (where some cores are already allocated to other processes); on a dedicated
Linux server or with `cpu_count-1` cores reserved, the speedup is typically closer
to 3–4×. On the full proteome (5177 hanpo × 6717 sce), single-threaded takes ~20 min
per direction; 4-threaded would take ~10–13 min per direction.

**Decision: ✓ implemented — `threads` changed to `max(1, os.cpu_count()-1)`.** BLAST
is deterministic across thread counts; this is a pure performance improvement with no
correctness risk.

---

## Thresholds in `get_model_from_homology`

These parameters act as post-BLAST filters on the homology table before mapping
genes from a template model to a new organism.

| Parameter | Default | Notes |
|---|---|---|
| `max_evalue` | `1e-30` | Very stringent; only strong alignments pass. This is ~25 orders of magnitude tighter than the BLAST E-value cutoff. |
| `min_align_len` | `200` | Minimum alignment length in amino acids (~60 aa per functional domain). |
| `min_identity` | `40` | Minimum percent identity. Below 40% the structural homology is uncertain. |

These values match MATLAB RAVEN and were chosen based on standard practice in
metabolic network reconstruction. They have not been benchmarked against a
gold-standard gene essentiality or ortholog dataset.

**Status: untested.** The most important parameter is `min_identity=40` —
the 40% threshold is a well-established heuristic (Doolittle 1981) for inferring
functional equivalence from sequence homology, but the right value depends on
the target organism pair.

Proposed benchmark: run `get_model_from_homology` on a yeast-GEM template with
H. polymorpha as the target organism at `min_identity` ∈ {30, 40, 50}, then
compare the number of genes mapped and the coverage of known H. polymorpha
metabolic functions (from the published hanpo-GEM).
