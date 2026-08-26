# Homology cut-off calibration

*Authored on `develop` (PR #92, 2026-08-25); copied here verbatim while that branch's*
*`get_model_from_homology`/`run_blast` changes have not yet been ported to this branch.*
*The reproduction script (`scripts/homology_cutoff_kegg.py`) lives on `develop` only —*
*see [reconstruction_homology.md](../maintenance/benchmarks/reconstruction_homology.md)*
*for how this bears on the current branch's (still-`200`) `min_align_len` default.*

Three settings decide which template reactions transfer to a new organism. None
of them had ever been measured, in either toolbox. This is what they are worth.

| Setting | Was | Now | Why |
|---|---|---|---|
| `min_identity` | 40 | **40** | The one that matters. 40 is right. |
| `min_align_len` | 200 | **100** | 200 was throwing away real matches for nothing. |
| `max_evalue` | 1e-30 | **1e-30** | Makes no difference anywhere between 1e-4 and 1e-50. |

`min_align_len` now differs from MATLAB RAVEN, which still uses 200. Worth
back-porting.

## How it was tested

Reconstruct four organisms from an *S. cerevisiae* template, at a range of
settings, and ask of every match that survives: are these two genes really
counterparts?

Two independent sources answer that — KEGG's gene annotations and OMA's ortholog
assignments — and wrong matches count double, because a wrong reaction is worse
than a missing one: gap-filling can add what is absent, while something wrong is
hard to notice and harder to remove.

The organisms span a range of relatedness, because a setting that suits a close
relative need not suit a distant one:

| | Organism | Relation to yeast |
|---|---|---|
| `kla` | *K. lactis* | close |
| `yli` | *Y. lipolytica* | moderate |
| `ani` | *A. nidulans* | distant |
| `eco` | *E. coli* | very distant |

Proteomes come from UniProt, relabelled with KEGG gene names so both references
can be applied to the same matches.

## Sequence identity: keep 40

Scored against KEGG annotations (higher is better):

| identity | `kla` | `yli` | `ani` | `eco` |
|---|---|---|---|---|
| 25 | 0.861 | 0.833 | 0.771 | 0.460 |
| 30 | 0.896 | 0.851 | 0.796 | 0.518 |
| 35 | 0.917 | **0.868** | **0.814** | 0.558 |
| **40** | **0.923** | 0.860 | 0.804 | **0.598** |
| 45 | 0.907 | 0.813 | 0.735 | 0.520 |
| 50 | 0.882 | 0.736 | 0.640 | 0.323 |

Against OMA the best value is 45 for the three fungi and 35 for *E. coli*. So the
two sources put the answer between 35 and 45, and 40 sits comfortably inside
that. Neither supports anything looser.

That last point is worth stating plainly, because counting a missing match as
equally bad as a wrong one moves the recommendation to 25 — a completely
different answer from the same measurements. Any recommendation about these
settings is meaningless unless it says how it weighed the two kinds of mistake.

## Alignment length: 200 was too strict

| length | `kla` | `yli` | `ani` | `eco` |
|---|---|---|---|---|
| 50 | 0.933 | 0.871 | 0.815 | 0.618 |
| **100** | **0.933** | **0.871** | **0.815** | **0.618** |
| 150 | 0.932 | 0.869 | 0.813 | 0.615 |
| 200 *(old)* | 0.923 | 0.860 | 0.804 | 0.598 |
| 300 | 0.885 | 0.816 | 0.764 | 0.571 |

Anything at or below 150 performs the same; the loss appears between 150 and
200. Dropping to 100 recovers 3–4% more real matches on every organism while the
wrong-match rate moves by at most 0.6 of a percentage point. 50 and 100 are
identical, so 100 was chosen as the less permissive of the two.

## E-value: not a setting worth touching

| `max_evalue` | `kla` | `ani` |
|---|---|---|
| 1e-4 … 1e-50 | 0.923 | 0.804 (identical throughout) |
| 1e-100 | 0.877 | 0.703 |

Five orders of magnitude, one answer. Identity and length have already removed
whatever a looser e-value would have removed, so it has nothing left to decide.
Only at 1e-100 does it start discarding good matches. Kept at 1e-30 to match
RAVEN.

## DIAMOND behaves the same as BLAST

Both toolboxes offer DIAMOND as the fast alternative and apply these same three
settings to its output, which is only safe if the two aligners behave alike.
They do:

| Organism | best with BLAST | best with DIAMOND | score at that point |
|---|---|---|---|
| `kla` | 40 | 40 | 0.923 vs 0.921 |
| `yli` | 35 | 35 | 0.868 vs 0.865 |
| `ani` | 35 | 35 | 0.814 vs 0.811 |
| `eco` | 40 | 40 | 0.598 vs 0.593 |

DIAMOND finds about half as many matches — 29,785 against 54,478 on `kla` — but
almost all the missing ones are weak ones these settings discard anyway. After
filtering, the two agree on 87–92% of what survives, and DIAMOND runs 10–20×
faster (34 s against 358 s on `kla`; 40 s against 788 s on `ani`).

So one shared set of defaults is justified, and DIAMOND is a reasonable choice
for large jobs.

## What the numbers do not cover

- **KEGG's annotations lean on BLAST comparisons**, so on their own they would
  partly measure agreement with the method being tested. OMA infers counterparts
  independently, which is why both were used; they agree.
- **The two sources are not directly comparable in absolute terms.** KEGG only
  covers genes it has annotated — a few hundred well-studied ones per organism —
  while OMA covers whole proteomes and lists strict counterparts, so ordinary
  gene duplicates count against us. Compare the shape of each column, not the
  heights between them.
- **The *Y. lipolytica* OMA row is the weakest**: 38% of its OMA pairs could not
  be matched to KEGG gene names, so that reference is incomplete.
- **One template.** Everything was measured outward from *S. cerevisiae*.
- **Measured on gene matches, not finished models.** A reaction transfers if any
  one of its genes matches, so the effect on a finished model is smaller than
  these numbers suggest.

## A curated-GEM reference was tried and retired

An earlier arm of this study also scored candidate thresholds against curated
non-model-organism GEMs (hanpo-GEM, rhto-GEM) built by RAVEN's own
homology-based reconstruction. It was dropped, and the reason is worth recording
because it generalises beyond this study.

The diagnostic: walk the thresholds loose one step at a time and ask what
fraction of newly admitted reactions are already in the curated model. The rate
held between 0.62 and 0.85 up to hanpo-GEM's actual build settings
(`max_evalue=1e-30, min_align_len=150, min_identity=35`) and collapsed to 0.06
one step past them. That cliff is the reference remembering its own
construction: the curated model agrees with predictions made at *its own* build
settings and disagrees with everything else, regardless of whether those other
settings are better or worse. Optimising against it returns the build settings
whether or not they are good — it cannot measure correctness, only
self-consistency. Most curated non-model fungal GEMs are RAVEN drafts built the
same way, so this applies to the class, not just these two models.

## Reproducing

```bash
python scripts/homology_cutoff_kegg.py fetch --out work/
python scripts/homology_cutoff_kegg.py align --out work/
python scripts/homology_cutoff_kegg.py score --out work/ --beta 0.5 \
    --gene-ko kegg118_organism_gene_ko.tsv.gz
```

Aligning is the only slow part and is saved, so re-scoring — under a different
weighting, or against a different reference — takes about two minutes.
