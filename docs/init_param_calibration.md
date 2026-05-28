# (f)tINIT parameter calibration & input-robustness

Empirical study of ravengem's (f)tINIT parameters on a genome-scale model (Human-GEM,
Hart2015 / HCT116). Two questions:

1. **Calibration** — on clean data, which parameter values give the best speed/quality
   trade-off? (`scripts/analyze_init_params.py`)
2. **Robustness** — with the task layer always on (it is part of the pipeline, not a
   variable), how does degrading the *transcriptomics input* affect the model, and which
   parameters keep it functional and stable? (`scripts/analyze_init_robustness.py`)

Both scripts are resumable and reusable on any model/dataset; the numbers below are HCT116.
"Jaccard" is reaction-set overlap with the reference (tightest setting / clean data) — for
a model-extraction tool the reaction set is the product, and a MIP gap bounds only the
*objective*, so set-stability is tracked separately.

---

## 1. Clean-data calibration

### ftINIT MILP — `mip_gap` (single step-0 solve, big_m=100, force_on=0.1)

| mip_gap | time (s) | objective | rel.obj.gap | Jaccard vs tightest |
|--------:|---------:|----------:|------------:|--------------------:|
| 0.0002 | 48 | 49357 | ref | ref |
| 0.001 | 44 | 49357 | +0.0000 | **1.0000** |
| 0.003 | 42 | 49289 | +0.0014 | 0.9973 |
| 0.01 | 42 | 49185 | +0.0035 | 0.9935 |
| 0.03 | 52 | 49185 | +0.0035 | 0.9935 |
| 0.1 | 46 | 45615 | +0.0758 | 0.9469 |

**Solve time is essentially flat across the gap** (the model build dominates), so a tight
gap is nearly free. `mip_gap=0.001` reproduces the proven optimum exactly (Jaccard 1.0);
quality only collapses at 0.1. → **Default 0.001.** (The genome-scale staged pipeline still
needs *some* gap + a `time_limit` because the full essential-forced MILP can be much harder
than this single step — see robustness timings.)

### ftINIT MILP — `big_m` (gap=0.001, force_on=0.1)

| big_m | time (s) | rel.obj.gap | Jaccard vs big_m=100 |
|------:|---------:|------------:|---------------------:|
| 100 | 51 | ref | ref |
| 50 | 54 | +0.0006 | 0.983 |
| 25 | 53 | +0.0007 | 0.982 |
| 250 | 55 | +0.0005 | 0.984 |
| 1000 | 59 | +0.0001 | 0.986 |

At step-0 (on the *scaled* model) `big_m` barely affects objective or time, but shifts which
reactions are kept by ~2% (alternate optima). `big_m=100` is RAVEN's value and is required
for the *staged* pipeline to stay feasible (a fixed 100 is only valid with stoichiometric
rescaling — see §1.4). → **Default 100.**

### ftINIT MILP — `force_on` (gap=0.001, big_m=100)

| force_on | time (s) | rel.obj.gap | Jaccard vs 0.1 |
|---------:|---------:|------------:|---------------:|
| 0.1 | 63 | ref | ref |
| 0.02 | 69 | +0.0005 | 0.983 |
| 0.05 | 56 | +0.0000 | 0.990 |
| 0.2 | 59 | +0.0004 | 0.982 |
| 0.5 | 79 | +0.0005 | 0.985 |

`force_on` (minimum flux for a reaction to count as "on") changes the *model*, not just a
tolerance, but the reaction set is fairly insensitive (Jaccard ≥0.98) and the objective
hardly moves. → **Default 0.1** (RAVEN), no strong reason to change.

### prep scaling — `rescaleModelForINIT` `max_stoich_diff` and on/off (gap=0.001, big_m=100)

| config | time (s) | rel.obj.gap | Jaccard vs scaled msd=25 |
|--------|---------:|------------:|-------------------------:|
| scale on, msd=25 | 51 | ref | ref |
| msd=10 | 49 | +0.0075 | 0.989 |
| msd=50 | 61 | +0.0003 | 0.982 |
| msd=100 | 62 | −0.0001 | 0.986 |
| scale off | 45 | +0.0129 | 0.973 |

At step-0 even `scale=off` is feasible, but it drifts most (Jaccard 0.973, objective +1.3%);
`max_stoich_diff` 10–100 are all within ~1%. **This understates scaling's importance** — at
step-0 there is no big-M cap on the held-out transports. In the *full staged pipeline*,
`scale=off` with `big_m=100` is **infeasible** (step-1 caps transports that step-0 used
freely). → **Keep scaling on, msd=25** (RAVEN's default).

**Calibration summary (defaults are well-chosen):** `mip_gap=0.001`, `big_m=100`,
`force_on=0.1`, scaling on (`max_stoich_diff=25`). For the genome-scale staged pipeline also
set a `time_limit` (≈120–600 s/step) so a hard essential-forced step returns a near-optimal
incumbent rather than grinding.

_tINIT (`run_init`) calibration (`mip_gap`/`eps`/`prod_weight`/`big_m`): pending (sweep
in progress)._

---

## 2. Robustness to degraded transcriptomics (task layer always on)

The metabolic-task + gap-fill layer is held fixed; only the expression input is degraded.
`frac` = fraction of the 69 essential tasks the extracted model performs (`check_tasks`);
`Jaccard` = reaction-set overlap with the clean-data model.

| input | n_rxns | tasks pass | frac | Jaccard vs clean |
|-------|-------:|-----------:|-----:|-----------------:|
| **clean** | 7777 | 69/69 | 1.000 | ref |
| dropout 50% | 5968 | 67/69 | 0.971 | **0.713** |
| dropout 70% | 5113 | 68/69 | 0.986 | **0.594** |
| noise σ=1.0 | 7812 | 69/69 | 1.000 | 0.952 |
| noise σ=2.0 | 7768 | 69/69 | 1.000 | 0.919 |
| downsample 50% | 6765 | 68/69 | 0.986 | 0.815 |
| downsample 70% | 6123 | 68/69 | 0.986 | 0.728 |

(dropout = genes set to 0 → score −5; noise = ×`exp(N(0,σ))`; downsample = genes dropped →
`no_gene_score`.)

**Findings:**

* **Robust to noise, sensitive to sparsity.** Multiplicative expression noise barely changes
  the model (Jaccard 0.92–0.95, size stable, all tasks pass). Sparsity is far more damaging:
  50% dropout already drops the reaction set to **0.71 Jaccard** (and shrinks 7777→5968), 70%
  to **0.59**.
* **Sparsity shrinks the model toward the task-essential core.** Missing/zeroed genes remove
  the expression evidence for a reaction; the task layer only adds back what tasks require, so
  sparse input yields smaller, more "generic" models. Dropout (−5) is harsher than
  downsampling (−2).
* **Functionality is largely but not perfectly preserved.** With the task layer, `frac` stays
  ≥0.97, but dips to 67–68/69 under heavy sparsity — i.e. the bounded gap-fill plus the
  post-hoc low-score-gene pruning occasionally leave 1–2 essential tasks unsatisfied. (See the
  lever sweep below for whether `no_gene_score`/`force_on` recover them.)
* **Cost tracks damage.** Dropout runs are slowest (more broken tasks → more gap-fill);
  noise is cheap.

> **Tractability note (a parameter that prevents failure):** the gap-fill MILP must be bounded
> (`mip_gap`/`time_limit`). Unbounded, severe degradation (which breaks many tasks at once)
> makes it solve a hard min-cost MILP per broken task to proven optimality — observed to run
> >75 min for one 90%-dropout model. With the bound it returns a near-optimal fill quickly.

### Levers at dropout 70% — which parameter best stabilises the model?

| config | n_rxns | frac | Jaccard vs clean |
|--------|-------:|-----:|-----------------:|
| default (no_gene_score=−2, force_on=0.1) | 5113 | 0.986 | 0.594 |
| no_gene_score=−1.0 | 5110 | 0.986 | 0.593 |
| no_gene_score=−0.5 | 5128 | 0.986 | 0.593 |
| force_on=0.2 | 5159 | 0.986 | 0.600 |

**No lever recovers the drift** — Jaccard stays ~0.59 across all settings. Two reasons,
both informative:

* The information dropout destroys is simply gone; no scoring/connectivity knob reconstructs
  the missing expression evidence. You cannot tune your way out of sparse input.
* `no_gene_score` is the wrong knob *for dropout specifically*: dropout leaves genes
  *present but zero* (scored −5), whereas `no_gene_score` only governs reactions whose genes
  are **absent** from the data — i.e. the *downsampling* failure mode. So `no_gene_score` is
  a meaningful lever for missing-data sparsity (a less-negative value keeps more
  unmeasured reactions, growing the model back toward clean), but it has nothing to act on
  under dropout.

**Practical takeaway.** The robustness levers that matter are *structural*, not numeric: the
task + gap-fill layer (keeps the model functional regardless of input quality) and a bounded
gap-fill MILP (keeps it tractable). For *missing*-gene sparsity specifically, `no_gene_score`
trades model size against confidence. For noise, defaults are already robust. No parameter
restores fidelity lost to dropout — that is a property of the data, not the pipeline.

---

## Reproducing

```bash
python scripts/analyze_init_params.py    --cell HCT116 --sweeps ftinit_milp,prep_scale,tinit,ftinit_full
python scripts/analyze_init_robustness.py --cell HCT116 --algo ftinit   # then --algo tinit
```

Both reuse the cached Human-GEM preps from the validation run
([docs/humangem_validation.md](humangem_validation.md)) and are resumable.
