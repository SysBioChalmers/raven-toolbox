# DeepLoc 2.1 vs iCre1355 (Chlamydomonas) compartmentalisation

A second independent non-yeast eukaryote, and the **most training-distant**: the green alga
*Chlamydomonas reinhardtii*. Like [AraCore](deeploc_aracore_benchmark.md) it exercises the
**chloroplast** yeast lacked, and adds the richest organelle set of the candidates. iCre1355 predates
DeepLoc 2 (no circularity) but is an **auto-generated** model (`iCre1355_auto.xml`), so its curation
is noisier than the hand-curated yeast-GEM / AraCore — read the numbers with that in mind.

* Driver: [`scripts/benchmark_deeploc.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/scripts/benchmark_deeploc.py) `--species icre1355`
* Model: iCre1355 (Imam *et al.*, *Plant J.* 2015,
  [doi:10.1111/tpj.13059](https://doi.org/10.1111/tpj.13059)) — 2394 reactions, compartments `c`
  Cytosol, `h` Chloroplast, `m` Mitochondria, `x` Glyoxysome, `n` Nucleus, `g` Golgi, `e`
  Extra-organism, plus `u` Thylakoid lumen, `f` Flagellum, `s` Eyespot, `i` inner-mito space.
* Scores: DeepLoc 2.1 (slow ProtT5) on the 1368 genes with an EnsemblPlants v5.5 sequence
  (`data/deeploc/icre1355/iCre1355_deeploc_00{1,2,3}.csv`; `Cre…` transcript-id headers).

DeepLoc maps to `c`, `h`, `m`, `x` (glyoxysome ← Peroxisome), `n`, `g`, `e`; the model's `er`-less,
lysosome-less scheme means DeepLoc's ER / Lysosome-Vacuole / Cell-membrane calls have no target. The
thylakoid lumen `u`, flagellum `f` and eyespot `s` are out of DeepLoc's scope (49 reactions total), so
addressable accuracy (49.3%) barely exceeds the as-is 47.9% — the low score is real, not an
addressability artefact.

## As-is (organelle call vs curation)

| compartment | n | correct | accuracy |
|---|--:|--:|--:|
| c (cytosol) | 887 | 349 | 39.3% |
| h (chloroplast) | 465 | 363 | 78.1% |
| m (mitochondria) | 243 | 79 | 32.5% |
| x (glyoxysome) | 42 | 31 | 73.8% |
| n (nucleus) | 21 | 1 | 4.8% |
| g (golgi) | 12 | 0 | 0.0% |
| u (thylakoid lumen) | 24 | 0 | 0.0% |
| f (flagellum) | 13 | 0 | 0.0% |
| s (eyespot) | 11 | 0 | 0.0% |
| **all** | **1718** | **823** | **47.9%** |

### Accuracy by DeepLoc confidence

| confidence bin | n | correct | accuracy |
|---|--:|--:|--:|
| [0.0, 0.5] | 195 | 84 | 43.1% |
| [0.5, 0.7] | 262 | 125 | 47.7% |
| [0.7, 0.9] | 541 | 247 | 45.7% |
| [0.9, 1.0] | 720 | 367 | 51.0% |

## Conclusions

* **The chloroplast generalises again (78.1%).** Across two independent photosynthetic eukaryotes —
  AraCore (89.9%) and now Chlamydomonas (78.1%) — DeepLoc recovers the plastid well, confirming the
  plastid result is not specific to land plants. The algal glyoxysome (a specialised peroxisome) also
  comes through at 73.8% via DeepLoc's Peroxisome call.
* **Cytosol (39.3%) and mitochondrion (32.5%) are poor — the hardest case in the suite.** Two causes,
  hard to separate: (i) Chlamydomonas is the furthest organism from DeepLoc's human/animal/plant
  training, and algal proteins are heavily **dual-targeted to chloroplast + mitochondrion**, so
  DeepLoc's single dominant call defaults to the chloroplast and a mito/cytosol reaction reads as a
  miss; (ii) `iCre1355_auto` is automatically compartmentalised, so its `c`/`m` truth is itself
  noisier than a hand-curated model. The flat confidence calibration (43% → 51%) — much weaker than
  yeast or AraCore — is consistent with DeepLoc being genuinely less certain on this proteome.
* **Net:** a useful stress test rather than a clean validation. The plastid/peroxisome signal
  transfers to algae; the cytosol/mito discrimination does not, and an auto-generated truth set caps
  how much can be concluded. The hand-curated [AraCore benchmark](deeploc_aracore_benchmark.md)
  remains the cleaner cross-kingdom result.

## Reproducing

```bash
python scripts/benchmark_deeploc.py --species icre1355 \
    --model /path/to/iCre1355_auto.xml \
    --csv data/deeploc/icre1355/iCre1355_deeploc_*.csv --doc /tmp/deeploc_icre1355.md
```
