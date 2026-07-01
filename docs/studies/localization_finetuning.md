# Finetuning localisation hyperparameters on yeast-GEM (slow DeepLoc)

Finetunes the DeepLoc-loading hyperparameters against curated yeast-GEM, using the **slow / ProtT5**
predictions (`data/deeploc/yeast-GEM_deeploc_*.csv`). Three hyperparameters the DeepLoc-vs-curation
data determines directly: the mitochondrial membrane-routing cut, the confidence gate, and the
triage per-compartment trust table.

* Driver: [`scripts/finetune_localization_yeast.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/scripts/finetune_localization_yeast.py)
  (`--organism 559292` adds the UniProt corroboration column).
* Model `yeastGEM_develop`, 1143 slow-model predictions, 2207 single-compartment GPR reactions.
* **Overfitting guard:** yeast-GEM's curation carries ~10% noise (two curated sources agree only
  ~90% — see [§4 of the DeepLoc benchmark](deeploc_yeast_benchmark.md)), so the confidence gate is
  judged on *corroboration against the consensus* (yeast-GEM **or** UniProt), never raw agreement
  with yeast-GEM alone.

## 1. Mitochondrial membrane routing (`membrane_threshold`)

### Mitochondrion matrix vs membrane (m/mm) -- routing threshold

Among 441 reactions DeepLoc placed in `m` (152 curated `m`, 289 curated `mm`): predict `mm` when
`1 - P(Soluble) >= threshold`. Balanced accuracy = mean of the two recalls.

| threshold | mm recall | m recall | balanced acc | Youden J |
|---|--:|--:|--:|--:|
| 0.30 | 100.0% | 55.3% | 77.6% | 0.553 |
| 0.35 | 100.0% | 78.9% | 89.5% | 0.789 |
| 0.40 | 100.0% | 83.6% | 91.8% | 0.836 |
| 0.45 | 100.0% | 86.2% | 93.1% | 0.862 |
| 0.50 | 100.0% | 86.8% | 93.4% | 0.868 |
| **0.55** | 100.0% | 90.8% | 95.4% | 0.908 |
| 0.60 | 85.5% | 91.4% | 88.5% | 0.769 |
| 0.65 | 85.5% | 92.1% | 88.8% | 0.776 |
| 0.70 | 19.0% | 92.8% | 55.9% | 0.118 |
| 0.75 | 17.0% | 94.1% | 55.5% | 0.110 |
| 0.80 | 15.2% | 97.4% | 56.3% | 0.126 |
| 0.85 | 15.2% | 98.7% | 57.0% | 0.139 |
| 0.90 | 12.5% | 99.3% | 55.9% | 0.118 |

**Optimal `membrane_threshold` = 0.55** (Youden's J = 0.908). The signal is sharply bimodal: every
curated `mm` reaction clears 0.55 (100% recall) while the matrix recall keeps climbing, so the whole
plateau `0.45–0.55` is ≥93% balanced. The current **default 0.50 sits inside that plateau**; 0.55 is
~2pp better but the gap is ~6 reactions on a single species' in-sample ROC, so the default is left at
0.50 (changing it to chase an in-sample optimum would over-fit). Above ~0.58 the signal collapses —
DeepLoc gives most inner-membrane proteins a membrane score in `[0.55, 0.6)`.

## 2. Confidence gate (`min_confidence`)

### `min_confidence` gate -- coverage vs accuracy (organelle-collapsed)

Drop reactions whose top gene is below the threshold; accuracy is on what remains.

| min_confidence | reactions kept | coverage | accuracy vs yeast-GEM | corroborated (yGEM or UniProt) |
|---|--:|--:|--:|--:|
| 0.0 | 2207 | 100% | 40.1% | 82.4% |
| 0.3 | 2207 | 100% | 40.1% | 82.4% |
| 0.4 | 2207 | 100% | 40.1% | 82.4% |
| 0.5 | 2207 | 100% | 40.1% | 82.4% |
| 0.6 | 2139 | 97% | 40.3% | 82.5% |
| 0.7 | 1760 | 80% | 42.4% | 88.3% |
| 0.8 | 1003 | 45% | 45.1% | 93.1% |
| 0.9 | 303 | 14% | 64.4% | 99.0% |

Read the **corroboration** column, not raw yeast-GEM accuracy (the latter is dragged down by the
structurally-unpredictable membrane sub-compartments, not by confidence). **`min_confidence = 0.7`
remains the recommended gate**: corroboration rises 82.4% → 88.3% while keeping 80% of genes. Push to
0.9 only when precision matters more than coverage (99% corroborated, but 14% kept). DeepLoc's
probability stays well calibrated on the slow model (monotone across every bin). The default is left
at 0.0 (opt-in) so the loader never silently drops genes.

## 3. Per-compartment trust (triage)

### Refreshed `DEEPLOC_COMPARTMENT_TRUST` (slow, organelle-collapsed)

Per-compartment reliability of a DeepLoc organelle call, from the slow run. Only `mm` inherits its
lumen's trust (the mitochondrial split is the one validated routing, AUC ~0.93); the other membranes
(`erm`/`gm`/`vm`) and `ce`/`lp` stay 0 -- DeepLoc cannot reach them reliably.

| compartment | n | accuracy (trust) |
|---|--:|--:|
| er | 439 | 0.88 |
| m | 518 | 0.86 |
| p | 115 | 0.83 |
| c | 698 | 0.79 |
| e | 9 | 0.78 |
| v | 67 | 0.36 |
| n | 44 | 0.18 |
| ce | 110 | 0.11 |
| g | 225 | 0.01 |
| lp | 146 | 0.00 |

```python
DEEPLOC_COMPARTMENT_TRUST = {
    "er": 0.88, "m": 0.86, "mm": 0.86, "p": 0.83, "c": 0.79, "e": 0.78, "v": 0.36, "n": 0.18, "ce": 0.11, "g": 0.01, "lp": 0.00, "erm": 0.00, "gm": 0.00, "vm": 0.00
}
```

The slow model reshuffles the trust order versus the earlier fast values: **mitochondrion jumps
0.67 → 0.86** (the slow model's big win — see the [yeast benchmark](deeploc_yeast_benchmark.md) §3),
vacuole 0.14 → 0.36, and `mm` becomes trustworthy (0 → 0.86) now that the validated split reaches it.
Golgi drops 0.23 → 0.01 (the slow model rarely calls Golgi correctly). This refreshed table is
applied as the default in
[`raven_toolbox.localization.triage`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/localization/triage.py).

## Applied

* **`DEEPLOC_COMPARTMENT_TRUST`** — refreshed to the slow values above (triage default).
* **`min_confidence`** — recommendation refreshed to 0.7 (corroboration 88.3%, 80% kept); loader
  default unchanged (0.0, opt-in).
* **`membrane_threshold`** — default left at 0.50 (inside the optimal plateau; 0.55 marginally better
  but in-sample on one species).

## Reproducing

```bash
python scripts/finetune_localization_yeast.py \
    --model /path/to/yeast-GEM/model/yeast-GEM.xml \
    --csv data/deeploc/yeast-GEM_deeploc_*.csv --organism 559292 --doc /tmp/finetune.md
```
