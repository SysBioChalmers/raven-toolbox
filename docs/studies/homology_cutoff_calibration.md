# Homology cut-off calibration

:::{admonition} Run 2026-08-24 — the defaults do not change
:class: important

Scored at **β = 0.5** — precision weighted above recall, because a wrongly
transferred reaction is worse than a missed one — `(1e-30, 200, 40)` stands.
`min_identity` 40 is optimal or within 0.01 of optimal on every organism tested.

One candidate change misses narrowly: `min_align_len` 50 improves all four
organisms at unchanged precision, but fails both pre-registered gates. It is
recorded rather than adopted.

Scored at β = 1 the same data says something quite different — a 10-point deficit
at distance — which is worth knowing about any threshold study that does not
state its weighting.
:::

## What this study found

1. **`max_evalue` does nothing.** From 1e-4 to 1e-50 the extracted model is
   byte-identical, on both ground truths and every organism. Only 1e-100 changes
   anything. Identity and alignment length have already excluded everything a
   looser e-value would admit.
2. **`min_align_len` barely matters below 150.** Values of 50 and 100 give
   identical results, and 50 scores best on every organism — a nearly free recall
   gain that the pre-registered gates nonetheless reject.
3. **`min_identity` is the only real lever, and 40 is the right value** under the
   agreed weighting. Which value looks best depends entirely on that weighting:
   at β = 1 the optimum moves to 25–30 and drifts with phylogenetic distance; at
   β = 0.5 it sits at 40 for everything.
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

| Target | Distance | Hits | Truth pairs | Default: precision / recall | F0.5 | F1 |
|---|---|---|---|---|---|---|
| `kla` | close | 54,478 | 1,100 | 0.976 / 0.761 | **0.923** | 0.855 |
| `yli` | medium | 57,492 | 1,091 | 0.967 / 0.596 | **0.860** | 0.737 |
| `ani` | distant | 76,973 | 1,131 | 0.936 / 0.514 | **0.804** | 0.663 |
| `eco` | very distant | 11,458 | 355 | 0.788 / 0.304 | **0.598** | 0.439 |

The default is precision-heavy and recall-poor, and gets more so with distance:
recall falls 0.76 → 0.60 → 0.51 → 0.30 across the ladder while precision holds
above 0.93 until the very-distant pair. Whether that is a flaw or the point
depends on the weighting — see below.

Both other parameters are near-inert, confirming arm 1 on independent data:

| Dimension (`kla` / `ani`) | Values | F1 |
|---|---|---|
| `max_evalue` | 1e-4 … 1e-50 | identical (0.855 / 0.663) |
| `max_evalue` | 1e-100 | 0.811 / 0.563 |
| `min_align_len` | 50, 100 | 0.874 / 0.682 |
| `min_align_len` | 200 *(default)* | 0.855 / 0.663 |
| `min_identity` | see below | dominant |

## Applying the criterion

The criterion, fixed before any data existed: *maximise hit-level F-score against
KO sharing, averaged over the medium and distant bands, subject to the
very-distant pairs transferring no more than the current defaults do*, and change
a default only if the improvement exceeds the spread across organisms within a
band.

### The loss function, decided on the argument

The protocol deliberately left the precision/recall weighting open, because it is
a value judgement rather than a measurement. It was settled on the asymmetry, not
by inspecting which threshold won: **a wrongly transferred reaction is worse than
a missed one.** A missing reaction can be recovered by gap-filling; a wrong one
pollutes the model and its gene associations silently and survives into
everything downstream. So **β = 0.5**, weighting precision above recall.

### Scored at β = 0.5

| ide | `kla` | `yli` | `ani` | `eco` | objective (`yli`, `ani`) |
|---|---|---|---|---|---|
| 25 | 0.861 | 0.833 | 0.771 | 0.460 | 0.802 |
| 30 | 0.896 | 0.851 | 0.796 | 0.518 | 0.823 |
| 35 | 0.917 | 0.868 | 0.814 | 0.558 | 0.841 |
| **40** *(default)* | **0.923** | 0.860 | 0.804 | **0.598** | 0.832 |
| 45 | 0.907 | 0.813 | 0.735 | 0.520 | 0.774 |

**The default identity of 40 is vindicated.** It is optimal on the close and
very-distant organisms outright, and within 0.01 of the best on the other two.
The F1-scored version of this table told a very different story — a 10-point
deficit on `ani` — and that deficit was entirely an artefact of weighting a
missed reaction as heavily as a wrong one. Under the weighting the maintainers
actually hold, the existing default is right.

### One candidate change, which the rules reject

Scoring the full grid rather than identity alone puts the optimum at
`min_align_len` **50** with identity unchanged, for every organism:

| Organism | default (40/200) | (40/50) | Δ | precision |
|---|---|---|---|---|
| `kla` | 0.923 | 0.933 | +0.009 | 0.976 → 0.976 |
| `yli` | 0.860 | 0.871 | +0.011 | 0.967 → 0.961 |
| `ani` | 0.804 | 0.815 | +0.012 | 0.936 → 0.937 |
| `eco` | 0.598 | 0.618 | +0.020 | 0.788 → 0.783 |

Recall rises 3–4 points while precision moves by at most 0.006 — nearly free, and
consistent across all four organisms. It nonetheless fails both gates:

- **Constraint**: `eco` calls rise 137 → 152, which is more than the defaults
  transfer. Violated, if narrowly.
- **Margin**: the objective improves by 0.011 against a within-band spread of
  0.056. Not met.

So the defaults stand — **`(1e-30, 200, 40)` unchanged** — and the honest note is
that the margin rule is doing questionable work here. It compares an improvement
against the spread between organisms of *different difficulty*, which is not a
noise estimate; a change that helps 4 of 4 organisms at unchanged precision is a
consistent signal, not sampling noise. Whether to adopt `min_align_len` 50 on
that basis is a maintainer's call, and one this study deliberately does not make
for itself.

## Recommendations

1. **Keep `(1e-30, 200, 40)`.** The identity default is confirmed by the agreed
   loss function; the other two parameters have no better value that clears the
   gates.
2. **Stop treating `max_evalue` as a tuning knob.** It is inert across five
   orders of magnitude in the default regime, on both ground truths. Document it
   and spend no further calibration effort there.
3. **Consider `min_align_len` 50** as a separate, explicit decision. It helps all
   four organisms at unchanged precision and fails only the margin and constraint
   rules, both of which are doing debatable work at this scale.
4. **State the weighting whenever these thresholds are discussed.** The same
   measurements recommend 25 or 40 depending on β alone.

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
