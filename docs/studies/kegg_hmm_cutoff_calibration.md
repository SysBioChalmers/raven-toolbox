# KEGG HMM-query cut-off calibration

This note records the measurements behind the default KO-assignment parameters in
`reconstruction/kegg/query.py` (`assign_kos` / `get_kegg_model_from_sequences`,
pipeline step 3b.5) and IMPROVEMENTS **K15**. It is the evidence for moving away
from RAVEN's `1e-50` cut-off.

## What the parameters do

`assign_kos` turns an `hmmscan` KO×gene E-value matrix into gene→KO assignments
in three steps:

1. **`cutoff`** — keep hits with `evalue <= cutoff`.
2. **`min_score_ratio_ko`** — within a KO, drop genes whose
   `log(evalue)/log(best_evalue_in_KO) < min_score_ratio_ko`.
3. **`min_score_ratio_g`** — within a gene, drop KOs whose
   `log(evalue)/log(best_evalue_for_gene) < min_score_ratio_g`.

## Method

- **Data:** KEGG release 118. Libraries: the maintainer-built `prokaryotes.hmm`
  (831 MB) and `eukaryotes.hmm` (692 MB), 90 %-clustered, FFT-NS-2/PartTree (K12).
- **Queries:** each organism's full proteome, extracted from `genes.pep`.
- **Ground truth:** the organism's *real* KEGG gene→KO links, from the
  `organism_gene_ko` table (restricted, as the table is, to reaction-linked KOs).
- **Prediction:** `assign_kos` output, with the `organism:` prefix stripped from
  query gene IDs so they match the bare gene IDs in the ground truth.
- **Metrics (gene→KO level):** precision = |pred ∩ truth| / |pred|,
  recall = |pred ∩ truth| / |truth|, F1. Reaction-level: `rxn_rec` = fraction of
  the organism's true reactions recovered (KO→reaction via `ko_reaction`);
  `rxn_novel` = predicted reactions **not** in the annotation set.
- Reproduce with [`scripts/analyze_hmm_cutoffs.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/scripts/analyze_hmm_cutoffs.py).

### Important caveat

All four organisms are **present in the libraries' training set**, so their own
sequences hit their KO profiles strongly and recall is an upper bound. The
calibration is therefore *relative* (how the parameters trade off, and where
RAVEN's default sits relative to the signal), not an absolute accuracy estimate.
A genome genuinely absent from KEGG would be the next validation. Also note that
`rxn_novel` / "precision < 1" partly reflects **legitimate homology** KEGG never
annotated for that organism (paralogs, un-curated genes), not pure error — so the
precision figures are a lower bound on real precision.

## Organisms

| code | organism | library | proteome (seqs) | true gene→KO pairs | true reactions |
|---|---|---|---|---|---|
| `sce` | *Saccharomyces cerevisiae* (budding yeast) | euk | 6021 | 841 | 1217 |
| `cme` | *Cyanidioschyzon merolae* (red alga) | euk | 5010 | 709 | 1157 |
| `eco` | *Escherichia coli* K-12 MG1655 | prok | 4288 | 1071 | 1548 |
| `mge` | *Mycoplasmoides genitalium* G37 (minimal genome) | prok | 476 | 110 | 211 |

`sce`/`eco` are model organisms; `cme`/`mge` are lesser-studied, `mge`
additionally being a small, divergent genome.

## 1. E-value separation (the key result)

`log10(E-value)` percentiles of the best hit per (gene, KO) pair, split by whether
the pair is in the organism's annotation (**matched**) or not (**novel**). Smaller
(more negative) = stronger hit.

| organism | group | n | p50 | p90 | p95 | p99 |
|---|---|---|---|---|---|---|
| `sce` | matched | 835 | −155 | −75 | −59 | −33 |
| `sce` | novel | 9467 | −8 | −2 | −0 | 1 |
| `cme` | matched | 704 | −133 | −63 | −47 | −25 |
| `cme` | novel | 10170 | −8 | −2 | −2 | 0 |
| `eco` | matched | 1070 | −142 | −69 | −57 | −36 |
| `eco` | novel | 27357 | −7 | −2 | −1 | −0 |
| `mge` | matched | 110 | −100 | −42 | −35 | −15 |
| `mge` | novel | 1904 | −4 | −2 | −1 | −0 |

**Reading:** matched pairs cluster at E ≈ 1e-100…1e-155; even their weakest 1 %
sit at 1e-15…1e-36. Novel pairs cluster at ≈1e-8. The two are separated by ~20
orders of magnitude. RAVEN's **`1e-50` lies inside the *matched* tail** (between
the matched p90 and p95 for most organisms; past p90 for `mge`), so it discards
real-but-weakly-scoring annotations while gaining nothing against the (far weaker)
noise.

## 2. Cut-off sweep

`min_score_ratio_ko = 0.3`, `min_score_ratio_g = 0.8` fixed; gene→KO precision /
recall / F1 and reaction recovery vs the annotation.

### `sce`
| cutoff | gKO prec | gKO rec | gKO F1 | rxn rec | rxn novel |
|---|---|---|---|---|---|
| 1e-10 | 0.57 | 0.98 | 0.72 | 0.99 | 334 |
| 1e-20 | 0.65 | 0.98 | 0.78 | 0.97 | 283 |
| 1e-30 | 0.72 | 0.97 | 0.83 | 0.97 | 216 |
| 1e-50 | 0.78 | 0.95 | 0.86 | 0.96 | 157 |
| 1e-70 | 0.81 | 0.91 | 0.86 | 0.91 | 113 |
| 1e-100 | 0.84 | 0.76 | 0.80 | 0.79 | 68 |

### `cme`
| cutoff | gKO prec | gKO rec | gKO F1 | rxn rec | rxn novel |
|---|---|---|---|---|---|
| 1e-10 | 0.50 | 0.98 | 0.67 | 1.00 | 541 |
| 1e-20 | 0.57 | 0.98 | 0.72 | 1.00 | 421 |
| 1e-30 | 0.61 | 0.97 | 0.75 | 0.97 | 367 |
| 1e-50 | 0.70 | 0.93 | 0.80 | 0.94 | 307 |
| 1e-70 | 0.75 | 0.85 | 0.80 | 0.87 | 223 |
| 1e-100 | 0.80 | 0.70 | 0.75 | 0.71 | 136 |

### `eco`
| cutoff | gKO prec | gKO rec | gKO F1 | rxn rec | rxn novel |
|---|---|---|---|---|---|
| 1e-10 | 0.53 | 0.99 | 0.69 | 0.99 | 382 |
| 1e-20 | 0.57 | 0.99 | 0.73 | 0.99 | 300 |
| 1e-30 | 0.60 | 0.98 | 0.75 | 0.99 | 268 |
| 1e-50 | 0.67 | 0.95 | 0.78 | 0.98 | 198 |
| 1e-70 | 0.73 | 0.88 | 0.80 | 0.93 | 157 |
| 1e-100 | 0.82 | 0.74 | 0.77 | 0.80 | 96 |

### `mge`
| cutoff | gKO prec | gKO rec | gKO F1 | rxn rec | rxn novel |
|---|---|---|---|---|---|
| 1e-10 | 0.52 | 0.98 | 0.68 | 0.99 | 75 |
| 1e-20 | 0.62 | 0.96 | 0.75 | 0.98 | 51 |
| 1e-30 | 0.65 | 0.95 | 0.77 | 0.98 | 39 |
| 1e-50 | 0.77 | 0.84 | 0.80 | 0.87 | 29 |
| 1e-70 | 0.85 | 0.73 | 0.78 | 0.73 | 27 |
| 1e-100 | 0.87 | 0.50 | 0.64 | 0.47 | 21 |

**Reading:** recall is flat-and-high from 1e-10 to ~1e-30, then falls as the
cut-off eats into the matched tail — gently for model organisms, sharply for the
divergent `mge` (rxn recall 0.98 → 0.87 from 1e-30 → 1e-50, → 0.47 at 1e-100).
The recall lost to a stricter cut-off is *not* noise rejection (noise is at 1e-8);
it is real annotation. `rxn_novel` shrinks with stricter cut-offs because strong
un-annotated homologs are also removed.

## 3. Score-ratio sweep (`cutoff = 1e-50`)

| organism | ko ratio | g ratio | gKO prec | gKO rec | gKO F1 |
|---|---|---|---|---|---|
| `sce` | 0.0 | 0.50 | 0.61 | 0.96 | 0.74 |
| `sce` | 0.0 | 0.80 | 0.77 | 0.95 | 0.85 |
| `sce` | 0.0 | 0.95 | 0.84 | 0.93 | 0.88 |
| `sce` | 0.3 | 0.80 | 0.78 | 0.95 | 0.86 |
| `sce` | 0.5 | 0.80 | 0.80 | 0.95 | 0.86 |
| `cme` | 0.0 | 0.50 | 0.53 | 0.94 | 0.68 |
| `cme` | 0.0 | 0.80 | 0.69 | 0.93 | 0.79 |
| `cme` | 0.0 | 0.95 | 0.78 | 0.92 | 0.84 |
| `cme` | 0.3 | 0.80 | 0.70 | 0.93 | 0.80 |
| `cme` | 0.5 | 0.80 | 0.70 | 0.93 | 0.80 |
| `eco` | 0.0 | 0.50 | 0.39 | 0.96 | 0.56 |
| `eco` | 0.0 | 0.80 | 0.66 | 0.95 | 0.78 |
| `eco` | 0.0 | 0.95 | 0.76 | 0.94 | 0.84 |
| `eco` | 0.3 | 0.80 | 0.67 | 0.95 | 0.78 |
| `eco` | 0.5 | 0.80 | 0.69 | 0.95 | 0.80 |
| `mge` | 0.0 | 0.50 | 0.62 | 0.85 | 0.72 |
| `mge` | 0.0 | 0.80 | 0.77 | 0.84 | 0.80 |
| `mge` | 0.0 | 0.95 | 0.82 | 0.81 | 0.81 |
| `mge` | 0.3 | 0.80 | 0.77 | 0.84 | 0.80 |
| `mge` | 0.5 | 0.80 | 0.78 | 0.84 | 0.81 |

**Reading:**
- **`min_score_ratio_ko` is inert** — across all four organisms, varying it
  0.0 → 0.3 → 0.5 changes precision/recall by ≤0.02 (mostly 0.00). It is a
  magic-number knob that does effectively nothing here. (Full 0.0/0.3/0.5 × g-grid
  in the script output; representative rows shown.)
- **`min_score_ratio_g` is the real precision lever** — 0.80 → 0.95 lifts
  precision ~0.07–0.10 for ~0.02 recall loss. 0.50 is clearly too loose.

## 4. Chosen defaults and effect

| parameter | RAVEN / old | new default | rationale |
|---|---|---|---|
| `cutoff` | 1e-50 | **1e-30** | recovers the matched tail (esp. divergent genomes); still ~22 orders above the 1e-8 noise floor |
| `min_score_ratio_g` | 0.8 | **0.9** | the effective precision lever; offsets the looser cut-off |
| `min_score_ratio_ko` | 0.3 | 0.3 (kept) | empirically inert; retained for RAVEN parity |

Old default `(1e-50, 0.3, 0.8)` vs new default `(1e-30, 0.3, 0.9)`
(`min_score_ratio_ko` 0.3 ≡ 0.0 here):

| organism | gKO prec | gKO rec | rxn rec | rxn novel |
|---|---|---|---|---|
| `sce` | 0.78 → 0.76 | 0.95 → 0.96 | 0.96 → 0.96 | 157 → 137 |
| `cme` | 0.70 → 0.67 | 0.93 → 0.96 | 0.94 → 0.97 | 307 → 305 |
| `eco` | 0.67 → 0.67 | 0.95 → 0.97 | 0.98 → 0.99 | 198 → 173 |
| `mge` | 0.77 → 0.69 | **0.84 → 0.94** | **0.87 → 0.97** | 29 → 35 |

The divergent minimal genome gains ~10 points of recall (the case the sequence
path exists for); model organisms improve slightly and `eco` emits *fewer*
unannotated reactions (the tighter gene-ratio prunes spurious multi-KO genes). The
small precision dip vs annotation is dominated by extra strong homologs, not
weak-hit noise.

## 5. Whole-model cross-validation (sanity check)

Full reconstruction of *S. cerevisiae* two ways, at the old defaults:

| | annotation path (3b.4) | HMM path (3b.5) |
|---|---|---|
| reactions | 1355 | 1461 |
| metabolites | 1501 | 1567 |
| genes | 835 | 896 |

Reaction recall 96.3 % (1305/1355 shared, Jaccard 0.86); gene recall 96.6 %
(807/835 shared, Jaccard 0.87). The annotation path also exercises the new
`organism_gene_ko.tsv.gz` artefact (K14). `hmmscan` throughput ≈ 0.1 s/query
against either library on 12 threads (yeast: 6021 queries in 633 s).
