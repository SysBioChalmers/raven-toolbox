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

## Does `method='chrr'` fix it? Yes on e_coli_core — but at a cost that may not scale

Same 4 chains × 300 samples on e_coli_core, `method='chrr'` instead of `'achr'`:

| | ACHR | CHRR |
|---|---:|---:|
| wall time | 34.7 s | **702.0 s** |
| reactions scored | 87 | 95 (0 excluded as constant) |
| R-hat median | 1.0071 | **1.0049** |
| R-hat p90 | 1.0369 | 1.0146 |
| R-hat max | 1.3028 | **1.0248** |
| R-hat > 1.01 | 48.3% | 17.9% |
| R-hat > 1.1 | 1.1% | **0.0%** |

CHRR converges properly here — the worst reaction (`ICDHyr`, 1.0248) doesn't
even reach the loose 1.1 threshold, versus ACHR's `EX_succ_e` at 1.30. This is
a real, substantial fix, not a marginal one.

The cost is the problem: **~20x slower on a 95-reaction model.** CHRR's
up-front max-volume-ellipsoid rounding step is the likely driver, and MVE
computation typically scales worse than linearly with dimension — so a naive
extrapolation of 20x to yeast-GEM's 4102 reactions (43x more reactions than
e_coli_core) could plausibly land anywhere from "worse than 20x" to much
worse, not better. At ACHR's already-measured 2524s for 4 genome-scale chains,
a proportional 20x would be ~14 hours — not attempted blind. See the bounded
probe below for what was actually measured.

## Follow-up: does reallocating the same ACHR budget help? No.

Same total step budget as the default config (`thinning × n_samples` =
30,000 either way), just distributed differently: `thinning=300,
n_samples=100` instead of `thinning=100, n_samples=300`. 4 chains, yeast-GEM.

| | default (t=100, n=300) | reallocated (t=300, n=100) |
|---|---:|---:|
| wall time | 2524.4 s | 3451.0 s |
| R-hat median | 1.1671 | 1.1636 |
| R-hat p90 | 1.6414 | 1.6542 |
| R-hat max | 9.9416 | 9.8345 |
| R-hat > 1.01 | 96.5% | 95.2% |
| R-hat > 1.1 | 67.5% | 67.2% |
| worst reaction | `r_0318` (9.94) | `r_0318` (9.83) |

Essentially no change — same worst reactions, same rough ordering, same
overall failure rate, and it took *longer* (57.5 min vs 42 min) despite equal
total steps. **This rules out "just thin more within a fixed budget" as a
fix.** If more thinning genuinely bought better mixing, spending the same
budget on longer gaps between fewer stored samples should have moved R-hat;
it didn't move it meaningfully in either direction. The non-convergence looks
structural to ACHR's mixing on this polytope, not a matter of turning an
existing dial — consistent with CHRR (a different algorithm entirely) fixing
it on e_coli_core while this reallocation, still ACHR, does not.

## Does CHRR fix it at genome scale, and is it practical?

Not run at matching scale (4 chains × 300 samples) given the ~20x
e_coli_core cost multiplier implies perhaps 14 hours. A small bounded probe
(2 chains, 20 samples, same `thinning=100`) was run instead purely to get a
real genome-scale CHRR timing number before deciding whether a full run is
worth attempting.

(cobrapy's `OptGPSampler` was not a candidate here: `random_sampling` doesn't
wire it in, only `'achr'` and `'chrr'` — see
[flux_sampling_algorithms.md](../reference/flux_sampling_algorithms.md).)

**Result: 4815.3 s (~80 min) for 2 chains × 20 samples.** This settles the
timing question on its own, independent of sample-count effects: CHRR's cost
at genome scale is dominated by a fixed per-chain cost (almost certainly the
max-volume-ellipsoid rounding step, computed once before any samples are
drawn) that doesn't shrink with a smaller sample request. A trivial 20-sample
probe already costs comparable wall time to a *full* 300-sample ACHR run
(2524 s). CHRR is not a cheap drop-in genome-scale fix with the current
implementation, regardless of what its converged quality would turn out to
be at a matching sample count.

The R-hat computed from this probe (median **5.16**, p90 86, max in the
billions) should **not** be read as "CHRR converges worse than ACHR at genome
scale." With only 20 samples per chain, within-chain variance (R-hat's
denominator) is estimated from too little data to be stable — a reaction
that happens to show near-zero variance in 20 draws by chance, combined with
any between-chain difference, produces an enormous, physically meaningless
ratio. (Contrast e_coli_core, where 300 samples/chain gave stable,
well-behaved R-hat throughout.) This run cannot distinguish "CHRR doesn't
work at genome scale" from "R-hat needs more than 20 samples to mean
anything" — telling those apart would need a genome-scale CHRR run with
enough samples for a stable R-hat, which circles back to the cost problem
above.

## Bottom line

- **Default settings are unconverged for most reactions at genome scale** —
  robust finding, confirmed by two independent lines of evidence (ESS and
  R-hat).
- **Reallocating the same ACHR budget doesn't help** — ruling out the
  cheapest possible fix.
- **CHRR fixes it on a small model, but its current genome-scale cost (a
  fixed ~80 min+ per chain before any samples are even drawn) makes it
  impractical as a drop-in fix today.** Whether CHRR *would* converge well at
  genome scale given enough samples to trust the R-hat is still open — it
  would need a run long enough to be informative, which is itself the
  problem.
- **No cheap, validated fix exists yet.** Users doing genome-scale flux
  sampling with `random_sampling`'s defaults should treat per-reaction flux
  ranges as unconverged for most reactions, not as a caveat affecting a
  minority. The most concrete unblock identified but not pursued here: CHRR's
  rounding transform is recomputed from scratch per chain/call — caching or
  reusing it across calls on the same model would remove the dominant fixed
  cost and is worth a future look, but is an engineering change, not a
  parameter default.

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
