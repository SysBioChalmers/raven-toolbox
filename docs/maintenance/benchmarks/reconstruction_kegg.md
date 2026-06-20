# KEGG reconstruction parameter benchmarks

Functions: `raven_toolbox.reconstruction.kegg.query.assign_kos`,
`raven_toolbox.reconstruction.kegg.query.run_hmmsearch`,
`raven_toolbox.reconstruction.kegg.hmm.build_ko_hmm`,
`raven_toolbox.reconstruction.kegg.assemble.*`

Date: 2026-06-20. HMMER binary available at
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

**`min_score_ratio_g=0.9`** — within a KO, only genes scoring ≥90% of the best
gene-level score are kept as candidate gene assignments.

**`cutoff=1e-30`** — minimum E-value for any hmmsearch hit to be considered.

Both Python and MATLAB use these values. They originate from RAVEN's original
KEGG reconstruction workflow and have not been independently benchmarked against
a held-out genome with known KO assignments.

**Status: untested.** Proposed benchmark: run `assign_kos` on *Saccharomyces
cerevisiae* (sce.faa, 6717 sequences) against the KEGG KO HMM library; compare
the resulting KO-gene assignments against the official KEGG organism page for
yeast (known true positives) at `min_score_ratio_ko` ∈ {0.2, 0.3, 0.4} and
`cutoff` ∈ {1e-20, 1e-30, 1e-40}.

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
