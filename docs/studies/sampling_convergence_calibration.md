# ACHR sampling: between-chain convergence

[`sampling.md`](../maintenance/benchmarks/sampling.md)'s existing `thinning`
result is a **single-chain** diagnostic: lag-1 autocorrelation / effective
sample size, measuring how independent consecutive samples are *within* one
Markov chain. It says nothing about whether that one chain actually reached
every part of the flux polytope, or got stuck mixing well inside a sub-region.
That needs multiple independent chains from different starting points and a
check that they agree — the Gelman-Rubin R-hat diagnostic.

For each reaction, R-hat compares between-chain variance to within-chain
variance across `n_chains` independent `random_sampling` runs (different
seeds). R-hat ≈ 1.0 means the chains agree; R-hat > 1.1 (the common
convergence threshold) or > 1.01 (the stricter one used in published MCMC
work) flags a reaction whose distribution still depends on where its chain
started.

* Driver: [`scripts/analyze_sampling_convergence.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/docs/parameter-defaults/scripts/analyze_sampling_convergence.py)
* Settings matched to the existing single-chain study for direct comparison:
  `n_samples=300`, `thinning=100`, `warmup=1000` (cobrapy/RAVEN defaults).
* `n_chains=4`, run in parallel via `ProcessPoolExecutor` (one process per
  chain — cobra models aren't guaranteed thread-safe for concurrent solves).

## e_coli_core (95 reactions) — methodology validation

4 chains × 300 samples, 34.7 s wall (parallel).

| | value |
|---|---:|
| reactions scored (non-constant) | 87 / 95 |
| R-hat median | 1.0071 |
| R-hat p90 | 1.0369 |
| R-hat max | 1.3028 |
| reactions with R-hat > 1.01 | 42 (48.3%) |
| reactions with R-hat > 1.1 | 1 (1.1%) |

Worst-converged reactions:

| reaction | R-hat |
|---|---:|
| `EX_succ_e` | 1.3028 |
| `EX_co2_e` | 1.0420 |
| `CO2t` | 1.0420 |
| `FUM` | 1.0390 |
| `EX_h_e` | 1.0383 |
| `TKT2`, `RPE`, `TKT1`, `TALA`, `G6PDH2r` | 1.0369 (tied) |

**Already informative at this scale.** Even on a 95-reaction textbook model,
one reaction — `EX_succ_e`, succinate exchange, a byproduct/overflow route —
clears the "not converged" threshold (R-hat 1.30) at the default settings, and
nearly half the reactions fail the stricter 1.01 bar. This is a genuinely
different failure mode from what the single-chain ESS result showed: it's not
that samples are autocorrelated *within* a chain, it's that independent chains
land on measurably different distributions for a subset of reactions —
consistent with a byproduct-secretion pathway that's rarely favoured and only
gets explored if a chain's random walk happens to wander into that corner of
the polytope.

## yeast-GEM (4105 reactions, genome scale)

4 chains × 300 samples, 2524.4 s wall (~42 min — slower than the naive
"~same as one chain" estimate; four Gurobi processes evidently contend for
resources on a 12-core machine rather than scaling for free).

| | value |
|---|---:|
| reactions scored (non-constant) | 3364 / 4105 |
| R-hat median | **1.1671** |
| R-hat p90 | 1.6414 |
| R-hat max | 9.9416 |
| reactions with R-hat > 1.01 | 3246 (96.5%) |
| reactions with R-hat > 1.1 | 2271 (**67.5%**) |

Worst-converged reactions: `r_0318` (9.94), `r_0307` (9.65), `r_1690` (8.89),
`r_1077` (8.59), `r_2625` (5.92), `r_1072` (4.21), `r_4015` (3.88), `r_2690`
(3.43), `r_1113` (3.42), `r_1648` (3.42).

**This is a materially worse picture than e_coli_core suggested, and worse
than the existing single-chain ESS finding implied on its own.** At genome
scale, with the exact default settings (`thinning=100`, `n_samples=300`,
`warmup=1000`), the *median* reaction already exceeds the 1.1 "not converged"
threshold — meaning independent chains disagree on where a typical reaction's
flux distribution sits, not just on a tail of hard cases. Two in three
reactions fail even the loose threshold; effectively none (3.5%) pass the
strict one.

This is consistent with, and sharpens, the existing single-chain result
(ESS≈12 effective samples from 300 stored at these settings): low ESS said
samples are highly autocorrelated within a chain; R-hat now shows that beyond
being autocorrelated, at genome scale the chains often haven't reached the
same distribution *at all* within 300×100=30,000 total ACHR steps. These are
two independent lines of evidence pointing the same direction, not a
restatement of one.

**Practical reading:** the existing per-reaction flux ranges reported by
`random_sampling` at default settings on a genome-scale model should not be
trusted as converged for the majority of reactions. The existing docstring
warning (increase `thinning`/`n_samples`, check ESS, or switch to
`method='optgp'`) was correctly directioned but understated the scale of the
problem — this justifies raising it from an FYI-level note to an explicit
warning with numbers attached.

## Open question this raises but doesn't answer

Does a fix that's cheap to *describe* (bigger `thinning`, bigger `n_samples`,
`method='optgp'`) actually bring genome-scale R-hat down to a reasonable
level, and at what cost in wall time? Not measured here — each additional
4-chain genome-scale configuration costs on the order of 40 minutes (see
timing below), so this is left as a deliberate next step rather than
open-endedly sweeping configurations in the same run.

## Reproducing

```bash
python scripts/analyze_sampling_convergence.py \
    --model /path/to/yeast-GEM.xml --out work/ \
    --n-chains 4 --n-samples 300 --thinning 100 --warmup 1000
```

Results are cached per (model, n_chains, n_samples, thinning, warmup) config,
so re-running with the same settings is instant; changing `--thinning` or
`--n-samples` re-runs and caches separately, letting the R-hat vs. thinning
trade-off be explored without re-deriving already-cached chains.

**Timing caveat:** the single-chain study's 841 s was extrapolated to "about
841 s wall for 4 parallel chains too" — that estimate was wrong by ~3x
(actual: 2524 s). Four concurrent Gurobi processes on a 12-core machine
evidently contend for resources rather than scaling for free; budget for that
when planning further sweeps at this scale.
