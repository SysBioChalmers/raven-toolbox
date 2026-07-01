# DeepLoc 2.1 vs Human-GEM gene localisations

A human **positive control** — and a lesson in **circularity**. Human is the core of DeepLoc 2.x's
training, so high agreement is expected and proves little on its own (pair it with the stringent
[AraCore](deeploc_aracore_benchmark.md) and [iCre1355](deeploc_icre1355_benchmark.md) tests). The
real value here is methodological: Human-GEM's own gene localisations are **partly DeepLoc-derived**,
so a naive benchmark grades DeepLoc against itself.

* Driver: [`scripts/benchmark_deeploc_humangem.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/scripts/benchmark_deeploc_humangem.py)
* Model: Human-GEM v2.0.0 (Robinson *et al.*, *Sci. Signal.* 2020,
  [doi:10.1126/scisignal.aaz1482](https://doi.org/10.1126/scisignal.aaz1482)). Gene-level: the truth
  is each gene's annotated location set in `model/genes.tsv` (`compartments`), with provenance in
  `compDataSource`.
* Scores: DeepLoc 2.1 (slow ProtT5) on 2839/2848 genes
  (`data/deeploc/humangem/Human-GEM_deeploc_00{1..6}.csv`; `ENSG…` headers). A gene's annotation is a
  set (e.g. `Nucleus;Cytosol`); DeepLoc predicts the single dominant location, so a call is correct
  when its compartment is **in** that set.

## Circularity guard (the headline)

439/2848 of Human-GEM's gene compartments were assigned **by DeepLoc 2 itself**
(`compDataSource == DeepLoc2`). Scoring against those grades DeepLoc on its own output; they are
dropped, and only SwissProt/CellAtlas-sourced genes are scored.

| gene set | n scored | agreement (DeepLoc top in annotated set) |
|---|--:|--:|
| **independent** (SwissProt/CellAtlas) | 2365 | **75.2%** |
| independent, addressable only | 2100 | **84.7%** |
| circular (DeepLoc2-sourced) | 436 | 93.8% |
| all (incl. circular) | 2801 | 78.1% |

The circular rows score **93.8%** (DeepLoc agreeing with DeepLoc) versus **75.2%** on the independent
set — a stark illustration of why provenance matters. Including them inflates the naive number by
+2.9 pp; the inflation would be larger on a model with a higher DeepLoc-sourced fraction.

**Addressability.** Dropping the circular rows also removes the entire `ce` (cell membrane) and `e`
(extracellular) truth vocabulary: in Human-GEM *every* gene annotated to those two was DeepLoc2-
sourced — the independent SwissProt/CellAtlas sources annotate no metabolic gene there. DeepLoc still
routes 265 independent genes to `ce`/`e`, which then cannot be corroborated and count as misses. On
the 7 compartments the independent truth can express, agreement is **84.7%** — the fair DeepLoc-skill
number; 75.2% is the conservative floor.

## Per-compartment (independent genes)

Precision by predicted label: among genes DeepLoc calls a compartment, how often it is in the gene's
annotated set.

| predicted | n | in annotated set | precision |
|---|--:|--:|--:|
| c (cytosol) | 811 | 667 | 82.2% |
| n (nucleus) | 246 | 228 | 92.7% |
| er (endoplasmic reticulum) | 357 | 302 | 84.6% |
| g (golgi) | 197 | 173 | 87.8% |
| m (mitochondrion) | 375 | 340 | 90.7% |
| p (peroxisome) | 45 | 38 | 84.4% |
| ly (lysosome) | 69 | 30 | 43.5% |
| ce (cell membrane) | 190 | 0 | 0.0% (off-vocabulary) |
| e (extracellular) | 75 | 0 | 0.0% (off-vocabulary) |

## Accuracy by DeepLoc confidence (independent genes)

| confidence bin | n | correct | accuracy |
|---|--:|--:|--:|
| [0.0, 0.5] | 13 | 9 | 69.2% |
| [0.5, 0.7] | 619 | 419 | 67.7% |
| [0.7, 0.9] | 1320 | 974 | 73.8% |
| [0.9, 1.0] | 413 | 376 | 91.0% |

## Conclusions

* **On its home turf DeepLoc is strong (~85% addressable, 91% at high confidence).** The major
  organelles — nucleus (93%), mitochondrion (91%), golgi (88%), ER (85%), cytosol (82%), peroxisome
  (84%) — are recovered with high precision. As expected for a training-distribution positive
  control, this is the ceiling, not a generalisation claim.
* **Lysosome is the weak organelle (43.5%)** — DeepLoc's `Lysosome/Vacuole` call corroborates the
  human lysosome annotation only ~half the time, the one clear soft spot among the addressable
  compartments.
* **The circularity correction is the transferable lesson.** 15% of Human-GEM's gene compartments are
  DeepLoc-derived and score 94% against DeepLoc — any model that bootstraps localisation from DeepLoc
  must be de-circularised before it can validate DeepLoc, and `compDataSource`-style provenance is
  what makes that possible. Confidence stays well calibrated (68% → 91%), so the `min_confidence`
  gate transfers to human.

## Reproducing

```bash
python scripts/benchmark_deeploc_humangem.py \
    --genes /path/to/Human-GEM/model/genes.tsv \
    --csv data/deeploc/humangem/Human-GEM_deeploc_*.csv --doc /tmp/deeploc_humangem.md
```
