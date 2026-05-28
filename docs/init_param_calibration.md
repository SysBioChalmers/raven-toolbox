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

### tINIT MILP (`get_init_model`, `essential_rxns=[]`, `time_limit=400s`)

**mip_gap** (eps=1, prod_weight=0.5):

| mip_gap | time (s) | n_kept | Jaccard vs gap=0.001 |
|--------:|---------:|-------:|---------------------:|
| 0.001 | 901 | 6024 | ref |
| 0.003 | 869 | 6036 | 0.991 |
| 0.01 | 595 | 5967 | 0.968 |

Tightening the gap costs ~50% more wall time on this MILP (unlike ftINIT step-0, build
doesn't dominate); a 1% gap is ~30% faster with ~3% reaction-set drift.
→ **`mip_gap=0.001`** for stability, **0.01** for a faster looser solve.

**eps** (gap=0.005, the connectivity-flux threshold — *changes the model*):

| eps | n_kept | Jaccard vs eps=1.0 |
|----:|-------:|-------------------:|
| 0.1 | 6119 | 0.952 |
| 0.5 | 6123 | 0.952 |
| 1.0 | 6064 | ref |
| 2.0 | 6090 | 0.960 |

Each `eps` value gives a slightly different model (Jaccard ~0.95 across the range); the
reaction-set spread is ~5%. `eps=1.0` is RAVEN's default; smaller values produce *slightly*
larger models (loosen the connectivity bar). Pick by what the data justifies — see the
caveat at the top of `init.py`.

**prod_weight** (gap=0.005, the metabolite-production reward — *changes the model*):

| prod_weight | n_kept | Jaccard vs 0.5 |
|------------:|-------:|---------------:|
| 0.0 | 5973 | 0.961 |
| 0.25 | 6015 | 0.974 |
| 0.5 | 6064 | ref |
| 1.0 | 6106 | 0.955 |

A higher `prod_weight` keeps slightly more reactions (rewards more connectivity); `0.5`
(RAVEN's default) is the middle of the range. Effect is modest (~5%).

**big_m** (gap=0.005, `None` = per-reaction `ub`):

| big_m | n_kept | Jaccard vs None |
|------:|-------:|----------------:|
| None (per-rxn ub) | 6064 | ref |
| 1000 | 6064 | **1.000** |
| 250 | 6114 | 0.953 |
| 100 | 6023 | 0.930 |

`big_m=1000` is identical to `big_m=None` here because the model's `ub` is ±1000 already
(so the per-reaction cap *is* 1000). Smaller fixed caps (250, 100) shift alternate optima
by 5–7% but do not change the objective. Unlike ftINIT, tINIT has *not* been run through
`rescaleModelForINIT`, so dropping `big_m` below 1000 may invalidate the LP feasibility
region for reactions whose own bound is larger — keep the default (per-reaction `ub`).

**tINIT calibration summary:** `mip_gap=0.001` (or 0.01 for ~30% speedup at ~3% drift);
`eps`, `prod_weight`, `big_m` defaults are fine — they all change the *model*, not just
tolerance, so tune by what the data and biology call for, not by these tables.

### ftINIT full pipeline (`ftinit`, series='1+1', no-task scaled prep, `time_limit=600s`)

| config | time (s) | n_kept | Jaccard vs gap=0.001 |
|--------|---------:|-------:|---------------------:|
| **mip_gap=0.001** (default big_m=100) | 346 | 7752 | ref |
| mip_gap=0.003 | 288 | 7748 | 0.993 |
| mip_gap=0.01 | 218 | 7746 | **0.995** |
| big_m=50 (gap=0.003) | 738 | 7799 | 0.974 |
| big_m=250 (gap=0.003) | 345 | 7766 | 0.977 |

Unlike the single-step ftINIT MILP in §1.1 (where build time dominated and the gap was
free), **the full pipeline does benefit from a looser gap**: `mip_gap=0.01` is ~37 %
faster than `0.001` with Jaccard 0.995 — essentially the same model. → **For genome-scale
ftINIT, `mip_gap=0.01` (or 0.005) is the sweet spot**; keep 0.001 only if exact
reproducibility matters more than a few minutes.

`big_m=50` is actually *slower* than the default 100 (738s vs 346s) — a tighter cap makes
the LP relaxation harder for borderline reactions; `big_m=250` is the same speed as 100
but shifts the reaction set ~2 %. → **Keep `big_m=100`** (RAVEN's value, what scaling is
designed for).

### tINIT + many task-essential reactions: a structural limitation

ftINIT's task layer (gap-fill) and tINIT's task layer (forcing `essential_rxns`) are
*not equivalent*. tINIT forces every essential reaction to carry `flux ≥ eps`. With
Human-GEM's 113 task-essential reactions (the validation set), the resulting steady-state
system is infeasible regardless of `eps`:

| essentials passed to `run_init` | result |
|---|---|
| 0 (the original validation call) | ✅ ok, 6024 reactions |
| 113 (merged-survivor IDs from `prep.essential_rxns`) | ❌ `infeasible` (proven by Gurobi presolve, ~330s) |
| 260 (pre-merge IDs from `find_task_essential_reactions` cache) | ❌ `infeasible` (~480s) |

Lowering `eps` (1.0 → 0.1) does **not** fix it; the issue is that 100+ reactions cannot
simultaneously each carry a fixed positive flux in their forced direction at steady state.
ftINIT avoids this by using an *adaptive* per-reaction forcing magnitude
(`min(0.99·|previous flux|, force_on)`) so each essential is forced at a value it
*actually carried* in a prior feasible solution. tINIT's one-size-fits-all `eps`
mechanism doesn't have that escape hatch.

**Practical takeaway.** For functional context-specific models on genome-scale data, use
ftINIT — the task layer (gap-fill, adaptive essential forcing) is what makes the pipeline
robust. tINIT remains useful for the small/no-essentials case (e.g. the
expression-only baseline in the validation), but pairing it with the full task-essential
set is a known incompatibility; the tINIT robustness study below is therefore reported
with `essential_rxns=[]`.

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

### tINIT robustness — `essential_rxns=[]` (the tINIT-without-task-layer picture)

For the reasons in §1.5, tINIT cannot accept the full task-essential set as forced
reactions; this section runs `get_init_model` with `essential_rxns=[]` to show the
realistic tINIT behaviour on the same degradation gradient — i.e. the *cost of not
having ftINIT's gap-fill safety net*.

| input | n_rxns | tasks pass | frac | Jaccard vs clean |
|-------|-------:|-----------:|-----:|-----------------:|
| **clean** | 6277 | **35/69** | **0.507** | ref |
| dropout 50% | 4910 | 23/69 | 0.333 | 0.673 |
| dropout 70% | 2807 | 21/69 | 0.304 | 0.408 |
| noise σ=1.0 | 6661 | 25/69 | 0.362 | 0.878 |
| noise σ=2.0 | 6146 | 24/69 | 0.348 | 0.869 |
| downsample 50% | 5006 | 24/69 | 0.348 | 0.722 |
| downsample 70% | 3541 | 19/69 | 0.275 | 0.515 |

**The headline contrast with ftINIT:**

| | ftINIT (task layer) | tINIT (no task layer) |
|---|---|---|
| clean | 7777 rxns, **69/69 (1.000)** | 6277 rxns, **35/69 (0.507)** |
| dropout 0.7 | 5113 rxns, **68/69 (0.986)**, J 0.594 | 2807 rxns, **21/69 (0.304)**, J 0.408 |
| noise σ=2.0 | 7768 rxns, **69/69 (1.000)**, J 0.919 | 6146 rxns, **24/69 (0.348)**, J 0.869 |
| downsample 0.7 | 6123 rxns, **68/69 (0.986)**, J 0.728 | 3541 rxns, **19/69 (0.275)**, J 0.515 |

* tINIT-without-gap-fill fails roughly **half the essential tasks even on clean data**;
  ftINIT-with-gap-fill passes them all. Under degradation tINIT collapses further (down
  to 19/69 at 70 % downsample), ftINIT stays ≥67/69 throughout.
* **Reaction-set drift is comparable** under noise (Jaccard 0.87 vs 0.92) but worse for
  tINIT under sparsity (0.41 vs 0.59 at 70 % dropout) because there's no gap-fill to
  re-add structurally needed reactions.

This is *not* a critique of the tINIT algorithm — classic INIT was designed for the
no-task-layer case. It is the empirical evidence for why ftINIT's design choices (task
+ gap-fill, adaptive essential forcing) are the right ones for genome-scale tissue
model extraction, and why tINIT is mostly useful here as a baseline.

#### tINIT levers at dropout 70%

| config | n_rxns | tasks pass | frac | Jaccard vs clean |
|--------|-------:|-----------:|-----:|-----------------:|
| default (prod_weight=0.5, eps=0.1) | 2807 | 21/69 | 0.304 | 0.408 |
| prod_weight=0.0 | 2791 | 21/69 | 0.304 | 0.416 |
| prod_weight=1.0 | 3386 | 22/69 | 0.319 | 0.485 |
| prod_weight=2.0 | 3888 | 21/69 | 0.304 | 0.458 |
| eps=0.5 | 2620 | 21/69 | 0.304 | 0.391 |
| eps=1.0 | 3311 | 22/69 | 0.319 | 0.460 |

Same conclusion as the ftINIT levers: parameter tuning can nudge (`prod_weight≥1.0`
or a larger `eps` modestly grows the model and lifts Jaccard from 0.41 to ~0.48), but
**no tINIT parameter recovers anything close to ftINIT's functionality** (22/69 at best
vs ftINIT's 67–69/69 at the same dropout). The gap-fill layer, not the parameter
choice, is what bridges the gap.

---

## 3. Cross-solver portability

See [init_solver_benchmark.md](init_solver_benchmark.md) for the genome-scale
solver comparison (Gurobi/HiGHS/GLPK) and [tests/test_init_solvers.py](../tests/test_init_solvers.py)
for CI parameterised over installed MILP backends. Headline: at genome scale only Gurobi
is viable today; HiGHS fails on an upstream optlang `hybrid_interface.clone()` bug; GLPK
ignores `configuration.timeout` on MIP and ran 1 h+ without converging. Toy-scale
correctness is portable (Gurobi + GLPK give identical verdicts on the unit-test
networks), so local development works without a Gurobi licence.

---

## Reproducing

```bash
python scripts/analyze_init_params.py    --cell HCT116 --sweeps ftinit_milp,prep_scale,tinit,ftinit_full
python scripts/analyze_init_robustness.py --cell HCT116 --algo ftinit   # then --algo tinit
```

Both reuse the cached Human-GEM preps from the validation run
([docs/humangem_validation.md](humangem_validation.md)) and are resumable.
