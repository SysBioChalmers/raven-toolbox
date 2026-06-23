# DeepLoc 2.1 vs yeast-GEM compartmentalisation

Model `yeastGEM_develop`, 1143 DeepLoc predictions; 2207 single-compartment GPR reactions (867 genes
covered, 1305 in DeepLoc-addressable compartments). DeepLoc 2.1 predicts 9 organelles and has no label
for the organelle membranes (erm/mm/gm/vm) or the lipid particle (lp).

> Sections 1–3 reproduce with `scripts/benchmark_deeploc.py --species yeast`; the UniProt cross-check
> (§4) and finetuning validation (§5) with `scripts/compare_localization_sources.py` and the
> `load_deeploc` / `combine_scores` options they motivate.
>
> **Model: slow (ProtT5) DeepLoc.** §1–3 below use DeepLoc 2.1's *slow / high-quality* model. It
> lifts organelle-collapsed accuracy from the fast (ESM1b) model's **54.6% → 64.6%** and
> mitochondrial-membrane recall from **47% → 86%** (see §3), so the slow model is preferred and the
> two should never be mixed across compared models. **§4–5 were computed on the earlier fast run** and
> are pending a slow refresh; their *qualitative* conclusions hold but the exact percentages predate
> the model switch.

## 1. As-is (organelle call vs curation)

### Per-compartment accuracy

| compartment | n | correct | accuracy |
|---|--:|--:|--:|
| c | 698 | 552 | 79.1% |
| m | 219 | 152 | 69.4% |
| mm | 296 | 0 | 0.0% |
| er | 90 | 55 | 61.1% |
| erm | 348 | 0 | 0.0% |
| p | 115 | 95 | 82.6% |
| ce | 110 | 12 | 10.9% |
| lp | 146 | 0 | 0.0% |
| n | 44 | 8 | 18.2% |
| g | 13 | 3 | 23.1% |
| gm | 52 | 0 | 0.0% |
| v | 7 | 0 | 0.0% |
| vm | 60 | 0 | 0.0% |
| e | 9 | 7 | 77.8% |
| **all** | **2207** | **884** | **40.1%** |

### Accuracy by DeepLoc confidence

Reaction's top-gene organelle confidence vs whether the call is correct (vs yeast-GEM alone; the
corroboration view in §4 — vs yeast-GEM *or* UniProt — runs higher).

| confidence bin | n | correct | accuracy |
|---|--:|--:|--:|
| [0.5, 0.7] | 447 | 138 | 30.9% |
| [0.7, 0.9] | 1457 | 551 | 37.8% |
| [0.9, 1.0] | 303 | 195 | 64.4% |

## 2. Can membrane-type recover the lumen/membrane split?

### Membrane-type discrimination of lumen vs membrane

AUC of DeepLoc's membrane signal (1 - P(Soluble)) separating curated lumen from membrane, among
reactions DeepLoc placed in the correct organelle.

| organelle | reactions | organelle recall | lumen / membrane | AUC |
|---|--:|--:|--:|--:|
| Mitochondrion (m/mm) | 515 | 86% | 152 / 289 | 0.931 |
| Endoplasmic reticulum (er/erm) | 438 | 88% | 55 / 331 | 0.391 |
| Golgi apparatus (g/gm) | 65 | 5% | 3 / 0 | nan |
| Lysosome/Vacuole (v/vm) | 67 | 36% | 0 / 24 | nan |

## 3. Collapsed membranes (mm->m, erm->er, gm->g, vm->v)

Collapsing also turns lumen/membrane-spanning reactions (previously multi-compartment, so excluded)
into single-compartment ones, so the truth set grows.

### Per-compartment accuracy on the organelle-collapsed model

| compartment | n | correct | accuracy |
|---|--:|--:|--:|
| c | 698 | 552 | 79.1% |
| m | 518 | 444 | 85.7% |
| er | 439 | 386 | 87.9% |
| p | 115 | 95 | 82.6% |
| ce | 110 | 12 | 10.9% |
| lp | 146 | 0 | 0.0% |
| n | 44 | 8 | 18.2% |
| g | 225 | 3 | 1.3% |
| v | 67 | 24 | 35.8% |
| e | 9 | 7 | 77.8% |
| **all** | **2371** | **1531** | **64.6%** |

## Conclusions

* **DeepLoc reproduces the major metabolic organelles well** — cytoplasm, peroxisome, extracellular
  and mitochondrion are recovered at ~70–85%. It is weak on nucleus, Golgi and vacuole, barely
  predicts the cell envelope (`ce` 12/110), and has no label for the lipid particle (`lp`).
* **The slow model fixes the mitochondrial membrane.** Given the correct organelle, DeepLoc's
  membrane signal separates mitochondrial matrix from membrane cleanly (AUC ~0.93), and the slow
  model now *recalls* the mitochondrion in 86% of m/mm reactions (vs 47% fast), lifting the collapsed
  `m` row to 85.7%. So `Mitochondrion + transmembrane -> mm` is a usable rule. For the ER it still
  does **not** discriminate (AUC ~0.39): DeepLoc flags almost all ER proteins as membrane-associated,
  so it cannot tell `er` lumen from `erm` membrane (and the slow model's heavier membrane bias even
  lowers raw `er`-lumen accuracy while leaving collapsed ER unchanged at ~88%).
* **Collapsing the organelle membranes into their lumen is the fair organelle-level target** and
  reaches **64.6%** (slow; 54.6% fast), by removing the four structurally-unpredictable membrane rows.
  The residual gaps are `ce`, `lp` and the under-recalled nucleus/Golgi/vacuole — predictor limits,
  not modelling ones. For a cross-kingdom check that every truth compartment is addressable, see the
  [AraCore benchmark](deeploc_aracore_benchmark.md) (80.3%, plastid included).

## 4. Cross-check against UniProt — yeast-GEM is not gold truth

> Computed on the earlier fast (ESM1b) DeepLoc run; pending a slow refresh. The conclusion (two
> curated sources agree only ~90%, so the benchmark carries baked-in noise) is model-independent.

yeast-GEM's curation is a *fair indication*, not ground truth, so we triangulate with a second curated
source (UniProtKB `Subcellular location`). Gene-level, organelle-collapsed:

| comparison | agreement |
|---|---|
| DeepLoc -> yeast-GEM | 78.6% (687/874) |
| DeepLoc -> UniProt | 80.8% (678/839) |
| **UniProt <-> yeast-GEM** (two curated sources) | **89.8% (574/639)** |

The two curated sources agree only ~90%, so a yeast-GEM benchmark carries ~10% baked-in noise. On the
639 genes annotated by all three, DeepLoc's top organelle agrees with **both** 74.5% of the time;
another **9.4%** agree with UniProt but not yeast-GEM (DeepLoc likely right, yeast-GEM the outlier),
3.3% agree with yeast-GEM only, and **12.8%** agree with neither (DeepLoc's true errors). So DeepLoc's
effective error is ~13%, not the ~21% the yeast-GEM benchmark alone implies. The two curated sources
themselves agree worst on `ce` (60%) and `v` (65%) — genuinely hard compartments.

DeepLoc's own confidence is well calibrated: corroboration (top organelle in yeast-GEM **or** UniProt)
runs 54% -> 67% -> 88% -> **97%** across confidence bins `[0,0.5) [0.5,0.7) [0.7,0.9) [0.9,1]`.

## 5. Finetuning the pipeline (and the overfitting guard)

Because no single source is authoritative, finetune toward the **consensus**, not yeast-GEM alone.
Data-backed levers, all in `raven_toolbox.localization`:

| lever | how | effect |
|---|---|---|
| `combine_scores([deeploc, uniprot, …])` | weighted-sum consensus, agreement reinforced | DeepLoc+UniProt vs independent yeast-GEM: 75.3% -> **78.7%** (fast; multi-source) |
| `load_deeploc(min_confidence=0.7)` | drop low-confidence genes | corroboration 82.4% -> **88.3%**, keeps 80% of genes (slow) |
| `load_deeploc(membrane_split={"m":"mm"})` | route Mitochondrion->mm when transmembrane | mm recall ~100% among correctly-placed mito reactions (slow); threshold-optimal at 0.55 |

The slow refresh of `min_confidence` and `membrane_split`, plus the finetuned triage
`DEEPLOC_COMPARTMENT_TRUST` (mitochondrion 0.67 -> 0.86), is in the dedicated
[localisation finetuning study](localization_finetuning.md). The membrane split is **mito-only** by
design (ER's AUC ~0.39 means routing `er`/`erm` would only relabel the majority class). Validate any
tuning against the consensus or a held-out reference — never against yeast-GEM alone, or you fit its
~10% curation noise (the `combine_scores` row is the one lever still on the fast run, pending a
multi-source slow rerun).
