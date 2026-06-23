# DeepLoc 2.1 vs yeast-GEM compartmentalisation

Model `yeastGEM_develop`, 1143 DeepLoc predictions; 2207 single-compartment GPR reactions (867 genes covered). DeepLoc 2.1 predicts 9 organelles and has no label for the organelle membranes (erm/mm/gm/vm) or the lipid particle (lp).

## 1. As-is (organelle call vs curation)

### Per-compartment accuracy

| compartment | n | correct | accuracy |
|---|--:|--:|--:|
| c | 698 | 529 | 75.8% |
| m | 219 | 147 | 67.1% |
| mm | 296 | 0 | 0.0% |
| er | 90 | 83 | 92.2% |
| erm | 348 | 0 | 0.0% |
| p | 115 | 94 | 81.7% |
| ce | 110 | 0 | 0.0% |
| lp | 146 | 0 | 0.0% |
| n | 44 | 7 | 15.9% |
| g | 13 | 3 | 23.1% |
| gm | 52 | 0 | 0.0% |
| v | 7 | 1 | 14.3% |
| vm | 60 | 0 | 0.0% |
| e | 9 | 7 | 77.8% |
| **all** | **2207** | **871** | **39.5%** |

## 2. Can membrane-type recover the lumen/membrane split?

### Membrane-type discrimination of lumen vs membrane

AUC of DeepLoc's membrane signal (1 - P(Soluble)) separating curated lumen from membrane, among reactions DeepLoc placed in the correct organelle.

| organelle | reactions | organelle recall | lumen / membrane | AUC |
|---|--:|--:|--:|--:|
| Endoplasmic reticulum (er/erm) | 438 | 89% | 83 / 307 | 0.409 |
| Mitochondrion (m/mm) | 515 | 47% | 147 / 97 | 0.923 |
| Golgi apparatus (g/gm) | 65 | 29% | 3 / 16 | 0.000 |
| Lysosome/Vacuole (v/vm) | 67 | 1% | 1 / 0 | nan |

## 3. Collapsed membranes (erm->er, mm->m, gm->g, vm->v)

Collapsing also turns lumen/membrane-spanning reactions (previously multi-compartment, so excluded) into single-compartment ones, so the truth set grows.

### Per-compartment accuracy on the organelle-collapsed model

| compartment | n | correct | accuracy |
|---|--:|--:|--:|
| c | 698 | 529 | 75.8% |
| m | 518 | 247 | 47.7% |
| er | 439 | 390 | 88.8% |
| p | 115 | 94 | 81.7% |
| ce | 110 | 0 | 0.0% |
| lp | 146 | 0 | 0.0% |
| n | 44 | 7 | 15.9% |
| g | 225 | 19 | 8.4% |
| v | 67 | 1 | 1.5% |
| e | 9 | 7 | 77.8% |
| **all** | **2371** | **1294** | **54.6%** |

## Conclusions

* **DeepLoc reproduces the major metabolic organelles well** - ER, cytoplasm, peroxisome, extracellular and mitochondrion are recovered at 50-90%. It is weak on nucleus, Golgi and vacuole, **never** predicts the cell envelope (`ce` 0/110), and has no label for the lipid particle (`lp`).
* **Membrane-type helps for the mitochondrial membrane, but not the ER.** Given the correct organelle, DeepLoc's membrane signal separates mitochondrial matrix from membrane cleanly (AUC ~0.92), so `Mitochondrion + transmembrane -> mm` is a usable rule. For the ER it does **not** discriminate (AUC ~0.41): DeepLoc flags almost all ER proteins as membrane-associated, so it cannot tell `er` lumen from `erm` membrane. A naive membrane-routing rule only *looks* good on ER because `erm` is the majority class.
* **Collapsing the organelle membranes into their lumen is the fair organelle-level target** and lifts overall accuracy from 39.5% to 54.6%, by removing the four structurally-unpredictable membrane rows. The residual gaps are `ce`, `lp` and the under-recalled nucleus/Golgi/vacuole - predictor limits, not modelling ones.

