# Homology reconstruction parameter benchmarks

Functions: `raven_toolbox.reconstruction.homology.blast.run_blast`,
`raven_toolbox.reconstruction.homology.blast.run_diamond`,
`raven_toolbox.reconstruction.homology.homology.get_model_from_homology`

Date: 2026-06-20 (evalue/threads sections below); `get_model_from_homology`
thresholds section updated 2026-08-26 with results from the 2026-08-25 study.
Binary: BLAST 2.17.0
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

| Parameter | Default (this branch) | MATLAB | Notes |
|---|---|---|---|
| `max_evalue` | `1e-30` | `1e-30` | Very stringent; only strong alignments pass. This is ~25 orders of magnitude tighter than the BLAST E-value cutoff. |
| `min_align_len` | `100` | `200` | Minimum alignment length in amino acids (~60 aa per functional domain). Lowered from `200` on 2026-08-26, see below. |
| `min_identity` | `40` | `40` | Minimum percent identity. Below 40% the structural homology is uncertain. |

**Status: measured, and ported to this branch's code.** See
[Homology cut-off calibration](../../studies/homology_cutoff_calibration.md) for
the full study (done on `develop`, PR #92, 2026-08-25; doc copied here 2026-08-26,
`min_align_len` default changed here the same day — see below for what's still
`develop`-only).

**Method, in brief:** reconstruct four organisms from an *S. cerevisiae* template
across a relatedness series (*K. lactis* close → *Y. lipolytica* moderate →
*A. nidulans* distant → *E. coli* very distant), and check each transferred gene
match against two references that don't share BLAST's biases with each other:
KEGG gene annotations and OMA ortholog calls. A wrong transfer is scored worse
than a missed one (β=0.5: gap-filling can add a missing reaction; a wrongly
transferred one is hard to notice and hard to remove).

**Results:**

| Parameter | Finding |
|---|---|
| `min_identity=40` | **Confirmed optimal.** Both KEGG and OMA put the best value between 35–45; 40 sits inside that band and wins outright on the closest and most distant organisms. |
| `min_align_len=200` | **Too strict — should be 100.** Performance is flat from 50–150 and drops between 150–200; 100 recovers 3–4 points of recall on every organism tested for ≤0.6 points of precision. Diverges from MATLAB (200), which was never itself measured. |
| `max_evalue=1e-30` | **Confirmed inert.** Identical model output from 1e-4 to 1e-50 — identity and length have already excluded whatever a looser e-value would admit. Kept at 1e-30 for RAVEN continuity; not worth tuning further. |

**The circularity problem was tested directly, not just avoided.** An earlier
arm of the same study scored thresholds against curated GEMs (hanpo-GEM,
rhto-GEM) built by this same homology function. It confirmed the concern raised
earlier in this doc empirically: walking thresholds loose one step at a time, the
fraction of newly-admitted reactions already present in the curated model held
between 0.62–0.85 right up to the curated model's own build settings, then
collapsed to 0.06 one step past them. That cliff is the reference recognising its
own construction — optimising against it just returns the build settings,
regardless of whether they're good. That arm was retired in favour of the
KEGG/OMA references above. See the study's "A curated-GEM reference was tried and
retired" section for the full account.

**Done on this branch (2026-08-26):** `min_align_len` in
[`homology.py`](../../../src/raven_toolbox/reconstruction/homology/homology.py)
changed `200` → `100`, matching the measured optimum. `run_blast`/`run_diamond`
`evalue` also changed `1e-5` → `1e-4` in the same pass (unifying with MATLAB;
no measured downstream effect either way — see the "Cross-toolbox parity
decisions" section of [index.md](index.md)).

**Still `develop`-only, not ported here:** `min_align_len=200` remains MATLAB
RAVEN's own default (a back-port is proposed but not this session's call to
make), and `develop`'s `review_identity` parameter (surfaces near-miss
candidates instead of silently discarding them) — a follow-on feature to this
study, out of scope for a parameter-default change.
