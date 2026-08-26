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

**Genome-scale answer (2026-06-2x, Human-GEM HCT116 — see
[init_param_calibration.md](../../studies/init_param_calibration.md)):** this was
run already; it just hadn't been cross-referenced from this file. Findings:
- ftINIT single step: solve time is flat across the gap (model build dominates),
  so `mip_gap=0.001` is nearly free and reproduces the tightest-gap model exactly
  (Jaccard 1.0).
- ftINIT full genome-scale pipeline: **does** benefit from a looser gap —
  `mip_gap=0.01` is ~37% faster than `0.001` at Jaccard 0.995 (essentially the
  same model).
- tINIT (`get_init_model`): `mip_gap=0.001` for stability, `0.01` for ~30% faster
  at ~3% reaction-set drift.

**Decision: keep `mip_gap=None`** (Gurobi's default `MIPGap≈1e-4` is already at
least as tight as the measured-good `0.001`). Docstring guidance updated in
`run_init`/`run_ftinit` to cite the measured genome-scale numbers above instead
of MATLAB's never-independently-tested `0.0004`.

---

## `time_limit` — MILP per-step wall-clock cap

**MATLAB default:** 5000 ms (5 s) per MILP step inside `run_init`/`run_ftinit`
**Python default:** `None` (no cap)

tINIT and ftINIT solve one MILP per expression category (ftINIT) or a single large
MILP (tINIT). MATLAB caps each solve at 5 seconds, returning the best solution found
within that time. For hard instances this can mean a sub-optimal reaction set.

**Measured solve times (2026-06-20, synthetic toy model):** < 0.1 s — not informative.

**Genome-scale answer:** also already measured, in the same
[init_param_calibration.md](../../studies/init_param_calibration.md) robustness
study. This is the one place the two toolboxes' choices both turn out to be
wrong in the same direction (uncapped) or too aggressive (MATLAB's fixed 5 s):
genome-scale solves in that study routinely took 42–901 s depending on
configuration, and one severely-degraded-input case ran **>75 minutes**
uncapped before being killed — MATLAB's 5 s would truncate genome-scale solves
far too early for a useful incumbent; Python's `None` has a real, observed
runaway-time failure mode under hard/degraded input. The study's own working
value for genome-scale ftINIT is a `time_limit` of **≈120–600 s/step**; its
tINIT harness used **400 s**.

**Decision: keep `time_limit=None` as the code default** (still correct for
small/medium models, and for genome-scale runs on clean data where a step-0
build dominates anyway — see `init_param_calibration.md` §1.1). Docstring
guidance in `run_init`/`run_ftinit` updated to cite the measured 120–600 s
(ftINIT) / 400 s (tINIT) working values instead of MATLAB's untested 5 s, and to
flag the >75 min uncapped failure mode explicitly so users degrading input
quality know to set a cap.

---

## `allow_excretion` in `get_init_model`

See [manipulation.md](manipulation.md) for the full benchmark. Summary:
- Effect is **zero** at default `prod_weight=0.5`
- Was inconsistent with `run_init`/`run_ftinit` (both default to `False`)
- **Decision: `get_init_model` default changed to `False`** — done 2026-06-20, `6f3b57c`

---

## `big_m` in `run_ftinit`

**Parameter:** `big_m=100.0` (Python and MATLAB)

**Benchmark (2026-06-21): yeast-GEM prep_init_model flux bounds after rescaling**

`prep_init_model` calls `rescale_for_init` to normalise stoichiometric coefficients
(mean |coeff| → 1 per reaction), then **explicitly resets all bounds to ±1000**.

```
yeast-GEM (4102 reactions):
  After prep_init_model: 3078 reactions
  After rescale_for_init: max UB (finite) = 1000.00
  Reactions with UB > 100 (finite): 3061 / 3078
```

At first glance this looks like `big_m=100` is too small (3061 reactions have
UB=1000 while big_m=100). But reading the `ftinit.py` module docstring explains
why it is intentional:

> `big_m` caps a *scored* reaction's flux in its on/off (direction) constraint —
> using a fixed 100 rather than the reaction's ±1000 bound keeps the LP relaxation
> tight (what makes the genome-scale MILP tractable). Free / essential reactions
> keep their real bounds.

The key points:
1. **big_m is not a flux maximum** — it is an LP relaxation tightener. A Big-M
   constraint `v ≤ big_m × y` that is smaller than the variable bound (1000)
   makes the LP relaxation closer to the integer solution, dramatically reducing
   solve time for genome-scale MILPs.
2. **Essential and free reactions are unaffected** — they retain ±1000 bounds.
   Only *scored* reactions are capped by big_m.
3. **Stoichiometric rescaling shifts typical fluxes to O(1)** — after normalising
   stoichiometry (mean |coeff| = 1 per reaction), the biologically relevant flux
   range shifts from O(1000) to O(1). `big_m=100` >> `force_on=0.1` (the minimum
   flux to count as "on"), so the indicator binary correctly distinguishes on (flux
   ≥ 0.1) from off (flux = 0) without being artificially binding.
4. **MATLAB RAVEN also uses big_m=100** — confirming this is an intentional design
   decision, not an oversight.

**Decision: ✓ keep `big_m=100.0`.** Intentional LP-tightening parameter that
matches MATLAB RAVEN. The `rescale_for_init` normalisation means that the effective
dynamic range of scored-reaction fluxes is O(1), not O(1000), making big_m=100
a valid and appropriate relaxation tightener. Free and essential reactions are
not constrained by big_m and retain their full ±1000 bounds.

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
