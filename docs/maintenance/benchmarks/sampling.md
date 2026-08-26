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

**Test results (yeast-GEM, n=300 samples, warmup=1000, Gurobi, 2026-06-20):**

| `thinning` | Lag-1 autocorrelation | Wall time (s) |
|---|---|---|
| 20 | 0.973 | 660 |
| 100 (default) | 0.926 | 841 |
| 500 | 0.849 | 927 |

**All three values show very high autocorrelation.** At thinning=500, consecutive
samples are still 85% correlated — far from independent.

The decay is slow: going from thinning=20 to thinning=100 (5×) reduces autocorrelation
by only 0.047; going from 100 to 500 (5×) reduces it by 0.077. Extrapolating, to reach
a lag-1 autocorrelation of ~0.3 (a rough threshold for near-independent samples) would
require thinning in the tens of thousands — taking weeks of compute time for 1000 samples
on yeast-GEM.

This is a known property of ACHR on large polytopes: mixing time scales with the number
of dimensions (reactions), and yeast-GEM's 4102-reaction polytope is 15–20× larger than
the cobrapy validation models (~200–300 reactions) on which thinning=100 was calibrated.

**Implications:**
1. For large genome-scale models (>2000 reactions), ACHR samples with the default
   thinning=100 are highly autocorrelated and should not be treated as independent.
2. Effective sample size (ESS) diagnostics should be run post-sampling.
3. For production analyses on yeast-GEM or Human-GEM, consider using cobrapy's
   `OptGPSampler` (which has better large-model mixing via a different hit-and-run
   direction scheme) rather than ACHR.
4. The `thinning` parameter trades compute time linearly for autocorrelation reduction —
   but the reduction is logarithmically slow, so large thinning values alone do not solve
   the problem at genome scale.

**Effective sample size (ESS) estimate:**
For an AR(1) process, ESS ≈ n × (1−ρ) / (1+ρ).
At thinning=100 (ρ=0.926): ESS ≈ 300 × 0.074 / 1.926 ≈ **12 effectively independent
samples from 300 stored samples**. To collect 100 effectively independent samples
from ACHR at thinning=100 on yeast-GEM would require ~2600 stored samples at a
wall time of ~7.5 hours.

**Decision: ✓ keep `thinning=100`** (unchanged from cobrapy upstream, which is the
correct upstream default to track). Add a docstring warning that the default is
calibrated for small models. On genome-scale models (>2000 reactions) the ACHR ESS
at thinning=100 is very low (~12 effective samples from 300 stored); users should
either dramatically increase thinning (1000+), increase n_samples to compensate,
check ESS diagnostics post-sampling, or switch to a different sampler
(`method='optgp'` via cobrapy's OptGPSampler).

**Follow-up (2026-08-26):** this ESS number is a *single-chain* diagnostic — it
doesn't check whether independent chains agree with each other. See
[sampling_convergence_calibration.md](../../studies/sampling_convergence_calibration.md)
for the between-chain (Gelman-Rubin R-hat) picture. On yeast-GEM at these exact
default settings (4 chains, `n_samples=300`, `thinning=100`, `warmup=1000`),
**the median reaction already fails the R-hat>1.1 "not converged" threshold**,
and two-thirds of all reactions fail it (96.5% fail the stricter 1.01 bar).
This is worse than the ESS number alone implies — it's not just that samples
are autocorrelated within a chain, it's that independent chains frequently
land on different distributions entirely within the default sample budget.
**Revised guidance: at genome scale, treat `random_sampling`'s default-settings
output as unconverged for most reactions, not as a caveat for a minority of
hard ones.** Whether increasing `thinning`/`n_samples` or switching to
`method='optgp'` actually fixes this (and at what wall-time cost) is not yet
measured — see the study's "Open question" section.

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
