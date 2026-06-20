# INIT / ftINIT parameter benchmarks

Functions: `raven_toolbox.init.init.run_init`, `raven_toolbox.init.ftinit.run_ftinit`,
`raven_toolbox.init.build.get_init_model`

Date: 2026-06-20.

---

## `mip_gap` — MILP optimality tolerance

**MATLAB default:** `0.0004` (passed to Gurobi `MIPGap`)
**Python default:** `None` (Gurobi's own default, typically `1e-4`)

The MILP at the core of tINIT and ftINIT minimises the number of reactions
added or removed relative to a template model while maximising consistency with
expression scores. A looser `mip_gap` allows the solver to terminate earlier with
a solution within that percentage of optimal.

**Parameters tested:** `None`, `0.0004`, `0.01`, `0.05`

Test model: synthetic tINIT testModel (linear chain A→B→C→D, 4 reactions with
gene rules g1/g2/g3). This model is too small to expose any difference in MIP gap
quality.

| `mip_gap` | Objective | Reactions kept |
|---|---|---|
| `None` | 21.0 | identical |
| `0.0004` | 21.0 | identical |
| `0.01` | 21.0 | identical |
| `0.05` | 21.0 | identical |

The test model is solved to the global optimum in < 0.1 s regardless of gap;
the MIP gap parameter only matters when solving a hard instance where the branch-
and-bound tree is large.

**Decision: keep `mip_gap=None`.** Gurobi's default `MIPGap=1e-4` is already tight.
MATLAB's `0.0004` is marginally looser (4× relative to optimal vs 1×), which
speeds up genome-scale solves at the cost of a small optimality gap. Document in
the docstring:

> For large genome-scale models (>5000 reactions) where solver time is a
> bottleneck, `mip_gap=0.0004` (MATLAB default) is a reasonable starting point.

**Still needed:** Run tINIT with real expression data on yeast-GEM at `mip_gap=None`
vs `0.0004`; compare model size (retained reactions) and quality (task satisfaction
rate).

---

## `time_limit` — MILP per-step wall-clock cap

**MATLAB default:** 5000 ms (5 s) per MILP step inside `run_init`/`run_ftinit`
**Python default:** `None` (no cap)

tINIT and ftINIT solve one MILP per expression category (ftINIT) or a single large
MILP (tINIT). MATLAB caps each solve at 5 seconds, returning the best solution found
within that time. For hard instances this can mean a sub-optimal reaction set.

**Measured solve times (2026-06-20, synthetic toy model):** < 0.1 s — not informative.

**Decision: keep `time_limit=None`.** For small-to-medium models the solver completes
in seconds. For genome-scale models, Python's uncapped solver will find a better
solution than MATLAB's 5-second cap at the cost of longer runtime. Document
MATLAB's 5 s as a starting point for genome-scale models where runtime is critical.

---

## `allow_excretion` in `get_init_model`

See [manipulation.md](manipulation.md) for the full benchmark. Summary:
- Effect is **zero** at default `prod_weight=0.5`
- Inconsistency with `run_init`/`run_ftinit` (both default to `False`)
- **Decision: change `get_init_model` default to `False`**

---

## `big_m` in `run_ftinit`

**Parameter:** `big_m=100.0` (Python and MATLAB)

`big_m` is the big-M constant used in the binary-indicator formulation of
ftINIT. It must be large enough that `big_m ≥ max(|v_i|)` for all reactions
in the model; otherwise the binary indicators do not correctly constrain reaction
activity.

yeast-GEM has a conventional upper bound of 1000 for most reactions, meaning
`big_m=100` would be **too small** if reactions carry fluxes up to 1000.

However, `big_m` in the ftINIT context is applied to a pre-processed model
(`prep_init_model`) where the model is rescaled. Checking `max_stoich_diff=25.0`
in `prep_init_model` and verifying whether rescaling brings fluxes below 100
is required before declaring this safe.

**Status: untested.** Requires a full ftINIT run on yeast-GEM with expression
data and verification that the binary constraints are not infeasible or wrongly
relaxed.

Proposed test: run `run_ftinit` on yeast-GEM at `big_m` ∈ {10, 100, 1000} and
check (a) feasibility, (b) objective value, (c) number of reactions retained.

---

## Scoring parameters (`factor`, `max_score`, `min_score`)

These control how RNA expression values are converted to INIT gene scores:

```
score = clip(factor × (log2(TPM) - log2(threshold)), min_score, max_score)
```

**Parameters:** `factor=5.0`, `max_score=10.0`, `min_score=-5.0`
(Python and MATLAB both use these values from Wang et al. 2012)

These are literature values from the original INIT publication and are consistent
across both implementations. No empirical test needed unless a new scoring
approach is proposed.
