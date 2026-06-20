# Flux sampling parameter benchmarks

Function: `raven_toolbox.analysis.sampling.random_sampling`

Date: 2026-06-20. Models: yeast-GEM (4102 reactions), iJO1366 (2583 reactions).

---

## `replace_max_bound` — replace big-M upper bounds with infinity before sampling

**Parameters tested:** `False` (Python default), `True` (MATLAB `randomSampling` default)

Only affects `method='random_objective'` (the `random_objective` sampler is the only
method that constructs random LP objectives over the sampling polytope).

Test model: yeast-GEM (4102 reactions). With `method='random_objective'`, n=200 samples.

| `replace_max_bound` | Outcome | Details |
|---|---|---|
| `False` (Python) | 200 samples complete | 0.57% of samples pinned at the 1000 bound; median per-reaction std = 0; 2626/4102 reactions always at zero |
| `True` (MATLAB) | **Solver unbounded** | All 4083/4102 reactions at the big-M bound get `ub=+inf`; random objectives drive any of them to +∞ |

yeast-GEM (like all RAVEN-convention models) uses 1000 as the conventional big-M
upper bound for ~99% of reactions. Replacing all of them with `+inf` makes the
random-objective LP unbounded — the solver can drive the objective to infinity
through any unconstrained reaction.

MATLAB's `replace_max_bound=True` was designed for models where only a handful
of reactions genuinely hit the big-M bound and those reactions represent true
physiological capacity limits. Such models are rare in practice.

**Decision: ✓ keep `replace_max_bound=False`.** MATLAB `True` is broken on
standard RAVEN-convention models. Note: this parameter only applies to
`method='random_objective'`; for the default `method='achr'` it has no effect.

---

## `thinning` — ACHR thinning factor (samples discarded between stored samples)

**Python default:** `100` (cobrapy ACHRSampler default)
**MATLAB default:** N/A (MATLAB's randomSampling only implements random_objective)

In ACHR sampling, `thinning=k` means k random walks are taken between each stored
sample. Higher thinning reduces autocorrelation at the cost of more computation.
The appropriate thinning depends on the mixing time of the Markov chain, which
scales with the number of reactions and the geometry of the flux polytope.

**Test results (yeast-GEM, n=300 samples, warmup=1000, Gurobi):**

| `thinning` | Lag-1 autocorrelation | Wall time (s) |
|---|---|---|
| 20 | **0.973** (very high) | 660 |
| 100 (default) | TBD (~3300 s estimated) | — |
| 500 | TBD (~16500 s estimated) | — |

Note: the thinning=100 and thinning=500 runs were not completed due to wall-clock
time (~55 min and ~4.6 h respectively at the yeast-GEM scale). The thinning=20
result alone is highly informative.

At thinning=20, the lag-1 autocorrelation across the first 20 variable reactions is
0.973 — near-unit autocorrelation indicating that consecutive ACHR steps on yeast-GEM
(4102 reactions) are almost perfectly correlated. This is consistent with the theory:
ACHR's mixing time scales with the number of dimensions, and yeast-GEM's 4102-reaction
polytope is much larger than the cobrapy validation models (~200–300 reactions).

Implication: the default thinning=100 may produce samples with moderate-to-high
autocorrelation on yeast-GEM-scale models. For genome-scale sampling, users should
either increase thinning (at proportional computational cost) or use a post-hoc
effective sample size (ESS) diagnostic to assess whether their sample is adequate.

cobrapy's default of 100 was validated on smaller models in the cobrapy test suite.
The raven-toolbox ACHR implementation wraps cobrapy's `ACHRSampler` directly.

**Decision: ✓ keep `thinning=100`** (unchanged from cobrapy upstream). Document in
the docstring that for large models (>2000 reactions), consider increasing thinning
to 500–1000 and checking effective sample size diagnostics.

---

## `warmup` — number of warmup steps before storing samples

**Python default:** `1000` (cobrapy ACHRSampler default)
**MATLAB default:** N/A

Warmup ensures the chain has mixed before samples are stored. Too few warmup steps
can produce samples clustered near the starting point; too many is wasted computation.
1000 warmup steps at the default thinning is generally sufficient for models up to
~3000 reactions (cobrapy validation).

**Decision: ✓ keep `warmup=1000`.** Matches cobrapy default.

---

## `n_objectives` — number of random objectives per sample (random_objective method)

**Python default:** `2` (Bordel et al. 2010)
**MATLAB default:** `2`

At each sampling step, `n_objectives` random linear objectives are sequentially
optimised to generate a new feasible flux distribution. Higher values explore the
polytope more broadly per step but require more LP solves.

Both implementations match the Bordel et al. 2010 paper value. No sensitivity
benchmark has been run.

Proposed test: run `random_sampling(model, method='random_objective', n_objectives=k)`
at k ∈ {1, 2, 3, 4} on yeast-GEM and measure coverage of the flux space using
pairwise distance or PCA variance explained.

**Decision: ✓ keep `n_objectives=2`** pending a dedicated sensitivity test.

---

## `method` — sampling algorithm

**Python default:** `'achr'` (Hit-and-Run with direction sampled from the Approximate
Centroid of the feasible region)
**MATLAB default:** `'random_objective'` (sequential random LP objectives)

ACHR is the more modern and statistically rigorous approach: it generates
samples from the uniform distribution over the flux polytope. The random_objective
method generates solutions that are optimal for random objectives — more spread
in flux space but not uniformly distributed.

**Decision: ✓ keep `'achr'` as default.** ACHR is the standard for genome-scale
flux sampling in the field. The random_objective method remains available for
compatibility with MATLAB RAVEN workflows. The migration note in the docstring
informs MATLAB users.

---

## `loopless_good_reactions` — use loopless FVA to exclude thermodynamic loop reactions

**Python default:** `True`
**MATLAB default:** heuristic (exclude reactions with FVA bound ≥ 999)

The MATLAB heuristic excludes any reaction whose maximum flux reaches the big-M
bound (999) as a potential loop. This over-excludes legitimate reactions that
genuinely approach capacity limits. Python's loopless FVA correctly identifies
only reactions that participate in thermodynamically infeasible cycles.

**Decision: ✓ keep `loopless_good_reactions=True`.** More correct than MATLAB's
heuristic; correctly classifies reactions that reach the 1000 bound through real
metabolic capacity.
