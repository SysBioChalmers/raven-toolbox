# Homology cut-off calibration

:::{admonition} Run 2026-08-24 — the defaults do not change
:class: important

The measurements are in. Under the success criterion fixed before the data
existed, **no change to `max_evalue` / `min_align_len` / `min_identity` is
justified**, and they stay as they are.

That is not the whole story: the evidence supports a looser identity on every
organism measured, and it was blocked by a constraint that this run suggests was
mis-specified. Rewriting the criterion after seeing the curves is exactly what
fixing it in advance was meant to prevent, so the criterion stands and the
decision it hinges on is stated openly at the end.
:::

## What this study found

1. **`max_evalue` does nothing.** From 1e-4 to 1e-50 the extracted model is
   byte-identical, on both ground truths and every organism. Only 1e-100 changes
   anything. Identity and alignment length have already excluded everything a
   looser e-value would admit.
2. **`min_align_len` barely matters below 150.** Values of 50 and 100 give
   identical results; the default of 200 costs ~2 points of F1.
3. **`min_identity` is the only real lever**, and its best value *moves with
   phylogenetic distance* — so no single global default is right for every
   reconstruction.
4. **One of the two planned ground truths is unusable**, and the study can
   demonstrate why rather than merely suspect it.

## What ran, and what did not

| Arm | Status |
|---|---|
| KEGG, cross-organism | **Ran.** Four targets on a distance ladder |
| Curated GEM pairs | **Ran, then invalidated** — see below |
| Model-level metrics | Not run. The hit-level result is decisive and the model level is buffered; see *Next* |

The KEGG arm was nearly blocked: the local KEGG 118 dump has the KO tables but no
proteomes, and KEGG's bulk FASTA is subscription-only. UniProt supplies both the
sequences and a KEGG cross-reference (`sce:YAL001C`), which makes the gene ids
line up exactly, so the arm ran after all.

## Inputs

| | |
|---|---|
| KEGG release | 118 (`organism_gene_ko`, 11,772 organisms) |
| Proteomes | UniProt reference proteomes, relabelled with KEGG gene ids |
| Aligner | BLAST+ 2.17.0 (raven-toolbox binaries, from raven-data) |
| raven-toolbox | 0.3.0 |
| Curated arm | hanpo-GEM v1.0.1, yeast-GEM and rhto-GEM as shipped in its `data/templateModels/` |

Proteome mapping coverage — the fraction of each proteome carrying a KEGG gene
id, since unmapped sequences are dropped:

| Organism | Proteome | Sequences | Coverage |
|---|---|---|---|
| `sce` (template) | UP000002311 | 6,021 | 99.2 % |
| `kla` close | UP000000598 | 5,045 | 99.9 % |
| `yli` medium | UP000001300 | 5,991 | 92.8 % |
| `ani` distant | UP000000560 | 10,281 | 97.3 % |
| `eco` very distant | UP000000625 | 4,141 | 94.0 % |

:::{admonition} Strain matters, and "most common" is not a strain
:class: note
A KEGG organism code names a strain. E. coli K-12 entries in UniProt carry
`eco:` (MG1655) **and** `ecj:` (W3110) cross-references in equal number, so
picking the most frequent prefix would have been a coin flip between two
genomes. The mapping filters on the exact code.
:::

## Arm 1: the curated GEM, and why it cannot be used

Reconstructing *H. polymorpha* from yeast-GEM and rhto-GEM and comparing against
the curated hanpo-GEM looked like the non-circular check. It is not.

Recall is capped first: of hanpo-GEM's 2,370 reactions, 2,239 exist in the
templates at all, so the **ceiling is 0.945**; the remaining 131 are manual
curation (methanol pathway, SLIME lipids) that no threshold can reach.

Sweeping thresholds against it, the best combination was
`(1e-30, 150, 35)` — recall 0.736, Jaccard 0.583. That is *exactly* the setting
hanpo-GEM was built with (`code/reconstructionProtocol.m:116`), which is either a
strong result or a worthless one, and the sweep alone cannot say which.

The discriminating test: loosen the thresholds one step at a time and ask what
fraction of the newly admitted reactions land in the curated model.

| thresholds | reactions gained | in curated | rate |
|---|---|---|---|
| ide 45 | 631 | 437 | 0.693 |
| ide 40 *(default)* | 308 | 250 | 0.812 |
| len 175, ide 37 | 255 | 159 | 0.624 |
| **len 150, ide 35** *(build settings)* | 75 | 64 | **0.853** |
| len 140, ide 33 | 32 | 2 | **0.062** |
| len 125, ide 30 | 60 | 1 | 0.017 |
| len 100, ide 25 | 19 | 1 | 0.053 |

The rate holds between 0.62 and 0.85 up to and including the build settings, then
**collapses by more than an order of magnitude immediately past them**. A smooth
decline would indicate genuine signal. A cliff at precisely the settings used to
build the reference is the reference remembering its own construction: reactions
reachable at or above those thresholds were placed there by the original draft
and survived curation; reactions only reachable beyond them were never in that
draft.

**Any optimisation against hanpo-GEM returns (1e-30, 150, 35) whether or not
those values are good.** The arm is discarded.

This generalises uncomfortably: most curated non-model fungal GEMs are themselves
RAVEN homology drafts from yeast-GEM, so the whole class of "curated GEM as
ground truth" is largely downstream of the tool being calibrated.

## Arm 2: KEGG, across a distance ladder

For every surviving hit between an *S. cerevisiae* gene and a target gene, does
KEGG assign them a shared KO? Only hits where **both** genes carry a KO
annotation are judged — KEGG's table covers reaction-linked KOs only (843 of
~6,000 `sce` genes), so a hit involving an unannotated gene is unjudgeable rather
than wrong, and counting it either way would be an invention.

| Target | Distance | Hits | Truth pairs | Default (1e-30/200/40) | Best F1 |
|---|---|---|---|---|---|
| `kla` | close | 54,478 | 1,100 | P 0.976 R 0.761 **F1 0.855** | ide 30 → **0.878** |
| `yli` | medium | 57,492 | 1,091 | P 0.967 R 0.596 **F1 0.737** | ide 25 → **0.815** |
| `ani` | distant | 76,973 | 1,131 | P 0.936 R 0.514 **F1 0.663** | ide 25 → **0.760** |
| `eco` | very distant | 11,458 | 355 | P 0.788 R 0.304 **F1 0.439** | ide 30 → **0.530** |

The default is precision-heavy and recall-poor, and it gets worse with distance:
recall falls 0.76 → 0.60 → 0.51 → 0.30 across the ladder while precision stays
above 0.93 until the very-distant pair.

Both other parameters are near-inert, confirming arm 1 on independent data:

| Dimension (`kla` / `ani`) | Values | F1 |
|---|---|---|
| `max_evalue` | 1e-4 … 1e-50 | identical (0.855 / 0.663) |
| `max_evalue` | 1e-100 | 0.811 / 0.563 |
| `min_align_len` | 50, 100 | 0.874 / 0.682 |
| `min_align_len` | 200 *(default)* | 0.855 / 0.663 |
| `min_identity` | see below | dominant |

## Applying the criterion

The criterion, fixed before any data existed: *maximise hit-level F1 against KO
sharing, averaged over the medium and distant bands, subject to the very-distant
pairs transferring no more than the current defaults do.*

| ide | objective F1 (`yli`, `ani`) | `kla` F1 | `eco` calls | `eco` precision | `eco` F1 |
|---|---|---|---|---|---|
| 25 | **0.788** | 0.866 | 487 (3.6×) | 0.435 | 0.504 |
| 30 | 0.785 | **0.878** | 384 (2.8×) | 0.510 | **0.530** |
| 35 | 0.758 | 0.873 | 256 (1.9×) | 0.602 | 0.504 |
| **40** *(default)* | 0.700 | 0.855 | **137** | **0.788** | 0.439 |
| 45 | 0.597 | 0.820 | 82 | 0.866 | 0.325 |

The objective improves by 8.8 points at ide 25 — and **every** loosening violates
the constraint, ide 35 included. Under the criterion as written, the defaults
stand.

### The criterion looks mis-specified

Stated plainly rather than quietly fixed. At ide 35, F1 improves on *all four*
organisms, the constraint organism included (`eco` 0.439 → 0.504), because the
extra calls at distance are substantially correct: recall on `eco` rises from
0.304 to 0.552. The constraint penalises finding true orthologs, which was not
its intent — it was meant to stop recall being maximised until everything
transfers, and F1 on the very-distant pair already does that job.

### But there is a real argument for the strict default

F1 weights precision and recall equally, and reconstruction does not. At `eco`,
precision falls from 0.788 to 0.602 at ide 35. A **missing** reaction can be
recovered by gap-filling; a **wrongly transferred** one silently pollutes the
model and its gene associations, and tends to survive into everything downstream.
If a false transfer costs more than a missed one, today's strict default is
defensible on exactly this evidence.

## What this leaves open

The next decision is not a measurement but a value judgement, and it belongs to
the maintainers: **what is the relative cost of a wrongly transferred reaction
versus a missed one?** Answer that and the loss function follows — an
F-beta with beta < 1, say — and the sweep above can be re-scored in minutes
without realigning anything.

That answer should be recorded *before* the re-scoring, as this criterion was.

## Recommendations that do not depend on that decision

1. **Stop treating `max_evalue` as a tuning knob.** It is inert across five
   orders of magnitude in the default regime. Document it; do not spend
   calibration effort on it.
2. **`min_identity` is the lever**, and its optimum moves with phylogenetic
   distance: ~30 for a close relative, ~25 for medium and distant. A single
   global default is a compromise, and the docs should say so, so users
   reconstructing from a distant template know which number to reach for.
3. **`min_align_len` 200 is slightly conservative** — 150 costs nothing measured
   here — but the effect is small enough not to justify a change on its own.

## Limitations

- **KEGG is not independent of BLAST.** KO assignments come from KOFAM HMMs and
  SSDB, and SSDB is all-against-all BLAST. This measures agreement with KEGG's
  homology calls, not correctness. A genuinely independent benchmark
  (OMA, OrthoDB, eggNOG — Quest-for-Orthologs style) would be the next
  improvement, and OMA's API is reachable.
- **Only reaction-linked KOs are judged**, a few hundred genes per organism, and
  the unjudgeable fraction is large (3,143 hits for `ani` at the default). The
  measured precision is precision *among annotated genes*, which is optimistic.
- **One template only.** Everything is measured from *S. cerevisiae* outward;
  a different template may behave differently.
- **Hit level only.** Reaction-level effects are buffered — a reaction transfers
  if any one of its genes hits — so the model-level consequences of these
  differences are smaller than the F1 gaps suggest, and were not measured.

## Reproducing

Scripts are under `scripts/` (see the study driver). One alignment per organism
pair is cached, and every threshold combination is post-processing of it:

| Pair | Hits | Alignment time |
|---|---|---|
| `sce`+`rhto` → hanpo | 118,939 | 872 s |
| `sce` → `kla` | 54,478 | 358 s |
| `sce` → `yli` | 57,492 | 431 s |
| `sce` → `ani` | 76,973 | 788 s |
| `sce` → `eco` | 11,458 | 360 s |

The whole study is about half an hour of alignment plus a few minutes of
sweeping — cheap enough to repeat whenever the loss function is decided.
