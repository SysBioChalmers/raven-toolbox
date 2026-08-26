# KEGG reconstruction parameter benchmarks

Functions: `raven_toolbox.reconstruction.kegg.query.assign_kos`,
`raven_toolbox.reconstruction.kegg.query.run_hmmsearch`,
`raven_toolbox.reconstruction.kegg.hmm.build_ko_hmm`,
`raven_toolbox.reconstruction.kegg.assemble.*`

Date: 2026-06-20 (threads section below); `assign_kos` score-ratio section
updated 2026-08-26 — study predates this file (pre-0.2.0 release) and was
already on this branch under `docs/studies/`. HMMER binary available at
`~/.cache/raven_toolbox/binaries/hmmer-3.4.0-windows-x86_64/hmmsearch.exe`.

---

## `threads` in `run_hmmsearch` and `build_ko_hmm`

**Python default:** `1`
**MATLAB default:** all available cores

`hmmsearch` and `hmmbuild` (called by `build_ko_hmm`) support `-cpu N` for
multi-core parallelism. HMMER is documented as deterministic across thread counts
for the Viterbi algorithm used in `hmmsearch --cut_tc`; small floating-point
differences can appear in E-value estimation across threads, but are below
the significance of the score cutoffs used here.

**Status: threads performance test not yet run.** All current HMMER benchmarks
are single-threaded.

**Decision: change to `max(1, os.cpu_count() - 1)`.** This is a pure performance
fix. On a modern 8-core laptop, single-threaded hmmsearch on the full KEGG KO
library (>26,000 HMMs) can take 30–60 minutes per proteome; multi-threaded cuts
this to ~5 minutes.

---

## `seq_identity` in `build_ko_hmm`

**Parameter:** `seq_identity=0.9` (Python and MATLAB)

Used by CD-HIT to cluster sequences within each KO before building the HMM.
At 90% identity, highly similar sequences are collapsed to one representative,
reducing HMM overfitting. The CD-HIT documentation recommends 0.9 as the
default for protein sequences.

**Decision: ✓ keep `0.9`.** Matches CD-HIT recommendation and MATLAB default.

---

## Score ratio cutoffs in `assign_kos`

**`min_score_ratio_ko=0.3`** — a gene is assigned to a KO if its hmmsearch
bit score is at least 30% of the best score for that KO. This is a relative
threshold that accounts for KO-specific variation in HMM length.

**`min_score_ratio_g=0.9`** — within a gene, only KOs scoring ≥90% of the best
KO-level score are kept as candidate assignments.

**`cutoff=1e-30`** — minimum E-value for any hmmsearch hit to be considered.

**These values already diverge from MATLAB RAVEN**, deliberately:
`cutoff` (RAVEN `1e-50` vs here `1e-30`) and `min_score_ratio_g` (RAVEN `0.8` vs
here `0.9`) were both changed based on the measurements below; only
`min_score_ratio_ko=0.3` matches RAVEN, because it was found to have no effect.

**Status: measured.** See
[KEGG HMM cut-off calibration](../../studies/kegg_hmm_cutoff_calibration.md) for
the full study — already present on this branch, reproducible via
`scripts/analyze_hmm_cutoffs.py`.

**Method, in brief:** run `assign_kos` on four organisms of varying study depth
(*S. cerevisiae*, *C. merolae*, *E. coli* K-12, and the minimal genome
*M. genitalium*) against the real KEGG release-118 KO HMM libraries, and compare
predicted gene→KO assignments to each organism's own curated KEGG annotations
(precision/recall/F1 at the gene→KO level, plus reaction-level recovery via
KO→reaction mapping).

**Results:**

| Parameter | Finding |
|---|---|
| `cutoff` (RAVEN `1e-50` → **`1e-30`**) | RAVEN's `1e-50` sits *inside* the tail of real (matched) hits, not at the noise boundary — matched E-values cluster at 1e-100…1e-155 (even the weakest 1% sit at 1e-15…1e-36), while noise clusters at ~1e-8, so there is a ~20-order-of-magnitude gap between them. `1e-50` discards real annotations for no gain against noise. Effect is worst on the divergent minimal genome (`mge`): reaction recall 0.87 at `1e-50` vs 0.98 at `1e-30`. |
| `min_score_ratio_g` (RAVEN `0.8` → **`0.9`**) | The real precision lever: 0.80→0.95 lifts precision ~0.07–0.10 for ~0.02 recall loss, consistently across all four organisms. 0.9 was chosen to offset the precision cost of loosening `cutoff`. |
| `min_score_ratio_ko=0.3` | **Confirmed empirically inert** — varying 0.0/0.3/0.5 changes precision/recall by ≤0.02 (mostly 0.00) on every organism. Kept at RAVEN's value since it does nothing either way. |

**Net effect of the new defaults** `(1e-30, 0.3, 0.9)` vs RAVEN's `(1e-50, 0.3, 0.8)`:
model organisms (`sce`, `eco`) see small precision/recall trade-offs in both
directions; the divergent minimal genome `mge` — the case this parameter path
exists for — gains ~10 points of reaction recall (0.87 → 0.97) for a similar
precision cost. Full per-organism tables are in the study.

**Caveat carried over from the study:** all four organisms are in the HMM
libraries' own training set, so recall is an upper bound — this calibrates how
the parameters trade off relative to each other and relative to RAVEN's default,
not an absolute accuracy figure on a genome KEGG has never seen. `rxn_novel`
(predicted reactions absent from the reference) is also a lower bound on real
precision, since some of it is legitimate homology KEGG simply hasn't
curated for that organism yet.

---

## Model assembly flags

| Parameter | Default | Notes |
|---|---|---|
| `keep_spontaneous` | `True` | Spontaneous reactions (marked `COMMENT: This reaction is spontaneous`) have no gene rule and are always included. Excluding them would block many real metabolic routes. ✓ keep |
| `keep_undefined_stoich` | `True` | Reactions with variable stoichiometry (e.g. `n` subunits). Excluding them loses pathways; including them requires manual curation. ✓ keep for draft reconstruction |
| `keep_incomplete` | `True` | Reactions where not all enzymes are known. Same reasoning as above. ✓ keep |
| `keep_general` | `False` | Overview-map reactions that aggregate many specific reactions into one lumped step. Including them produces double-counting. ✓ keep `False` |

These match MATLAB RAVEN and reflect well-established reconstruction practices.
No empirical benchmark required.
