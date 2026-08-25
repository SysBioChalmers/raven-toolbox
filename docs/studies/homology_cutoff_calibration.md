# Homology cut-off calibration

:::{admonition} Run 2026-08-24 — `min_align_len` drops to 100
:class: important

Three filters decide which template reactions transfer to a new organism, and
none of them had ever been measured. Scored against KEGG orthology with precision
weighted above recall:

* **`min_identity` 40** — confirmed. It is the filter that matters and 40 is the
  measured optimum.
* **`min_align_len` 200 → 100** — changed. 200 was discarding real orthologs for
  no gain in precision.
* **`max_evalue` 1e-30** — unchanged, but it does nothing: any value from 1e-4 to
  1e-50 gives the identical model.
:::

## What the parameters do

`get_model_from_homology` filters the hit table before deciding what transfers:

```python
keep = (hits.evalue <= max_evalue) & (hits.align_len >= min_align_len) & (hits.identity >= min_identity)
```

The aligner has its own e-value (1e-4, matching RAVEN's `getBlast`), but that is
a collection threshold twenty-six orders looser than the filter above, so it
never binds. These three are the decision.

## How correctness was judged

Reconstruct a target organism from an *S. cerevisiae* template and ask, for every
hit that survives the filters, whether KEGG assigns the two genes a shared KO.
That is the question the filters actually answer — is this pair orthologous —
and it is measured before any of it reaches a model.

Judging at the hit level rather than the reaction level matters. A reaction
transfers if *any one* of its template genes has a surviving hit, so
reaction-level agreement is buffered: it would look flat across a wide range of
thresholds for reasons that have nothing to do with the thresholds being right.

Only hits where **both** genes carry a KO annotation are judged. KEGG's table
covers reaction-linked KOs only — 843 of ~6,000 yeast genes — so a hit involving
an unannotated gene is unjudgeable rather than wrong, and scoring it either way
would be an invention.

### The weighting, and why it decides the answer

A wrongly transferred reaction is worse than a missing one: gap-filling can
recover something absent, while something wrong is hard to spot and harder to
remove once it is in the model and its gene rules. Precision therefore counts for
more than recall, and results are scored at **β = 0.5**.

This is not a detail. At β = 1, the same measurements say identity 40 is costing
10 points at distance and should drop to 25. At β = 0.5 they say 40 is correct.
Same data, opposite conclusion — so a cut-off recommendation that does not state
its weighting is not a recommendation.

## Organisms

A distance series from the template, because the right answer depends on how far
the target sits from it. A single target would give a number that generalises to
nothing.

| Organism | Distance from yeast | Proteome | Sequences | KEGG-mapped |
|---|---|---|---|---|
| `sce` *S. cerevisiae* | template | UP000002311 | 6,021 | 99.2 % |
| `kla` *K. lactis* | close | UP000000598 | 5,045 | 99.9 % |
| `yli` *Y. lipolytica* | medium | UP000001300 | 5,991 | 92.8 % |
| `ani` *A. nidulans* | distant | UP000000560 | 10,281 | 97.3 % |
| `eco` *E. coli* | very distant | UP000000625 | 4,141 | 94.0 % |

Proteomes come from UniProt and are relabelled with KEGG gene ids using UniProt's
own cross-reference, which is what makes the KEGG annotations usable as truth.

:::{admonition} A KEGG organism code names a strain
:class: note
E. coli K-12 entries carry `eco:` (MG1655) **and** `ecj:` (W3110) cross-references
in equal number, so taking the most frequent prefix would have been a coin flip
between two genomes. The mapping filters on the exact code.
:::

## Results

At the current defaults, the filters are precision-heavy and increasingly
recall-poor with distance:

| Target | Hits | Truth pairs | Precision | Recall | F0.5 |
|---|---|---|---|---|---|
| `kla` | 54,478 | 1,100 | 0.976 | 0.761 | 0.923 |
| `yli` | 57,492 | 1,091 | 0.967 | 0.596 | 0.860 |
| `ani` | 76,973 | 1,131 | 0.936 | 0.514 | 0.804 |
| `eco` | 11,458 | 355 | 0.788 | 0.304 | 0.598 |

### Sequence identity: 40 is right

| identity | `kla` | `yli` | `ani` | `eco` |
|---|---|---|---|---|
| 25 | 0.861 | 0.833 | 0.771 | 0.460 |
| 30 | 0.896 | 0.851 | 0.796 | 0.518 |
| 35 | 0.917 | **0.868** | **0.814** | 0.558 |
| **40** *(default)* | **0.923** | 0.860 | 0.804 | **0.598** |
| 45 | 0.907 | 0.813 | 0.735 | 0.520 |
| 50 | 0.882 | 0.736 | 0.640 | 0.323 |

40 wins outright on the closest and most distant organisms and trails by under
0.01 on the two in between. Loosening to 35 buys recall at a precision cost that
the weighting does not accept; tightening to 45 loses real orthologs fast.

### Alignment length: 200 was too strict

| length | `kla` | `yli` | `ani` | `eco` | precision (`yli`) |
|---|---|---|---|---|---|
| 50 | 0.933 | 0.871 | 0.815 | 0.618 | 0.961 |
| **100** *(new default)* | **0.933** | **0.871** | **0.815** | **0.618** | 0.961 |
| 150 | 0.932 | 0.869 | 0.813 | 0.615 | 0.962 |
| 200 *(old default)* | 0.923 | 0.860 | 0.804 | 0.598 | 0.967 |
| 300 | 0.885 | 0.816 | 0.764 | 0.571 | 0.970 |

Everything at or below 150 measures the same; the loss appears between 150 and
200. Dropping to 100 recovers 3–4 points of recall on every organism while
precision moves by at most 0.006 — real orthologs that 200 was throwing away for
nothing.

50 and 100 are identical to three decimals, so 100 is chosen as the less
permissive of two equivalent values.

### E-value: inert

| `max_evalue` | `kla` | `ani` |
|---|---|---|
| 1e-4 … 1e-50 | 0.923 | 0.804 (identical throughout) |
| 1e-100 | 0.877 | 0.703 |

Five orders of magnitude, one answer. Identity and length have already excluded
whatever a looser e-value would admit, so the parameter has nothing left to
decide. Only at 1e-100 does it start discarding good hits. It stays at 1e-30 for
continuity with RAVEN, and is not worth tuning.

## What changed

| Parameter | Was | Now | Why |
|---|---|---|---|
| `min_identity` | 40 | **40** | Measured optimum; the filter that matters |
| `min_align_len` | 200 | **100** | 200 discarded real orthologs at no gain in precision |
| `max_evalue` | 1e-30 | **1e-30** | Makes no difference between 1e-4 and 1e-50 |

`min_align_len` now diverges from MATLAB RAVEN's 200. The evidence is here, so
the back-port is worth proposing rather than leaving the two toolboxes silently
different.

## Limitations

- **KEGG is not fully independent of BLAST.** KO assignments come from KOFAM HMMs
  and SSDB, and SSDB is built from all-against-all BLAST. This measures agreement
  with KEGG's orthology calls rather than ground truth. A curated ortholog set
  such as OMA would be a genuinely independent check.
- **Only annotated genes could be judged** — a few hundred per organism. The
  absolute precision figures are therefore flattering; the comparison *between*
  settings is still fair, since every setting is judged on the same genes.
- **One template.** Everything is measured outward from *S. cerevisiae*.
- **Hit level only.** Reaction-level effects are buffered, so the consequences
  for a finished model are smaller than these gaps suggest.

## Reproducing

```bash
python scripts/homology_cutoff_kegg.py fetch --out work/
python scripts/homology_cutoff_kegg.py align --out work/
python scripts/homology_cutoff_kegg.py score --out work/ --beta 0.5 \
    --gene-ko kegg118_organism_gene_ko.tsv.gz
```

Alignment is the only slow step — 6 to 13 minutes per organism pair — and is
cached, so re-scoring under a different weighting takes about two minutes and
needs no realignment.
