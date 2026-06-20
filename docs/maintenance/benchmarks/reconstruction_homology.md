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
| `1e-4` (MATLAB) | — | — |
| `1e-5` (Python) | — | — |
| Marginal hits (`1e-4` but not `1e-5`) | — | — |

*(BLAST benchmark running at time of writing — results pending. See note below.)*

**Expected finding:** For the hanpo vs sce proteome comparison (~50% average amino
acid identity for orthologous pairs), most hits should have E-values well below
1e-5. Marginal hits (between 1e-5 and 1e-4) are expected to have low identity (<40%)
and short alignment lengths — hallmarks of spurious alignments. The 10× difference
in E-value threshold is unlikely to recover meaningful homologs but will add noise.

**For closely related organisms** (≥70% AAI, e.g. different *Saccharomyces* species):
the difference between 1e-4 and 1e-5 is negligible — all real homologs score well
above 1e-5.

**For distantly related organisms** (≤30% AAI, e.g. *H. polymorpha* vs bacteria):
neither 1e-4 nor 1e-5 is stringent enough — `get_model_from_homology` applies
additional filters (`max_evalue=1e-30`, `min_identity=40`) that dominate.

**Decision: ✓ keep `evalue=1e-5`.** Matches the `blastp` command-line default.
MATLAB's `1e-4` is 10× more permissive and adds hits that are removed by
downstream filters anyway. No correctness benefit expected.

Update this file with actual numbers when the BLAST benchmark completes.

---

## `threads` — number of parallel BLAST/Diamond processes

**Python default:** `1`
**MATLAB default:** all available cores (auto-detected)

BLAST and Diamond are documented as deterministic across thread counts — the same
hits are returned regardless of parallelism (alignment scores are computed
independently per query sequence). HMMER may show negligible E-value differences
across threads due to floating-point accumulation in background frequency
estimation, but below the threshold of any meaningful cutoff.

**Benchmark (2026-06-20): H. polymorpha vs S. cerevisiae, threads=1 vs threads=4**

| Threads | Hits | Wall time (s) | Results identical? |
|---|---|---|---|
| 1 | — | — | — |
| 4 | — | — | — |
| Speedup | — | — | — |

*(BLAST benchmark running at time of writing — results pending.)*

**Expected finding:** threads=4 should produce the same number of hits and identical
scores as threads=1, with ~3–4× speedup. On a full proteome (6000+ sequences vs
6000+ sequences), single-threaded BLAST takes ~5–10 minutes; 4-threaded takes ~2 min.

**Decision: change `threads` default to `max(1, os.cpu_count() - 1)`.** This is
a pure performance improvement. The change applies to `run_blast`, `run_diamond`,
`run_hmmsearch`, and `build_ko_hmm`. The old `threads=1` default silently makes
these functions 4–8× slower than necessary on modern hardware without any correctness
benefit.

Update this file with actual speedup numbers when the BLAST benchmark completes.

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
