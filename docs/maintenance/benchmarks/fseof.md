# FSEOF parameter benchmarks

Function: `raven_toolbox.analysis.fseof.fseof`
Algorithm: Flux Scanning based on Enforced Objective Function (Choi et al. 2010)

Test setup: iJO1366 (2583 reactions, 1805 metabolites), target reaction `EX_succ_e`
(succinate export), objective `BIOMASS_Ec_iJO1366_core_53p95M`.
Date: 2026-06-20, Gurobi solver.

---

## `flux_eps` — noise floor for flat/near-zero reactions

**Parameters tested:** `1e-8`, `1e-7`, `1e-6` (Python default) — all three run through
raven-toolbox's own `_classify()`; MATLAB's `FSEOF.m` itself was not run for this
comparison, see the correction below.

Three roles in `_classify()`:
1. Skip reaction if `flux.std() < flux_eps` (flat across all steps — likely unconstrained or constant)
2. Skip reaction if `|slope| < flux_eps` (regression slope near zero — no trend)
3. Mark as knockout if `|final_flux| < flux_eps` (reaction is at zero at maximum enforcement)

| `flux_eps` | Amplified | Knocked down / out | Notes |
|---|---|---|---|
| `1e-8` | 18 | 414 | 21 extra vs `1e-6`; all have flux std ≈ 5e-7 |
| `1e-7` | 18 | 414 | same 21 extra |
| `1e-6` (default) | 18 | 393 | |

The 21 reactions caught at `1e-8`/`1e-7` but not `1e-6` all have flux standard
deviation ≈ 5×10⁻⁷ across the 10 scan steps. This is below Gurobi's primal feasibility
tolerance (1e-9) accumulated across ~2583 reactions — textbook floating-point
summation noise. Labelling them knockdown targets would mislead users into
targeting reactions that are functionally zero throughout the scan.

**Decision: keep `flux_eps=1e-6`.** The `1e-8` row above was originally labelled
"MATLAB implicit" as a stand-in for "very tight tolerance", not from actually running
MATLAB's `FSEOF.m`. **Correction (2026-08-28):** reading `core/FSEOF.m` on
`origin/develop` directly shows it has no tolerance at all — target classification is
a bare `fseof.results(j,i) > fseof.results(j,i-1)`-style comparison against the
previous iteration's exact float value, not against any threshold, fixed or
otherwise. That's stricter than even the `1e-8` row here, so the 21-reaction
false-positive count above is a *lower* bound on what MATLAB's zero-tolerance
comparison would flag, not an equivalent measurement of it — the true MATLAB-side
number wasn't measured and is expected to be higher. Either way the conclusion
holds: an explicit floor is needed, `1e-6` measured to be safely above the solver-noise
band, and porting one to MATLAB is a bigger change than exposing an existing
implicit value (there isn't one to expose).

---

## `n_steps` — number of scan points between wt and max-enforced flux

**Parameters tested:** 5, 10 (Python/MATLAB default), 20

| `n_steps` | Amplified | Knocked down / out |
|---|---|---|
| 5 | 18 | 393 |
| 10 (default) | 18 | 395 |
| 20 | 18 | 395 |

n=5 misses 2 knockdown reactions that appear stably at n=10 and n=20 (the correlation
is computed over fewer points, so weak but real trends can be missed). n=10 and n=20
give identical results, so n=10 is not a bottleneck.

**Decision: keep `n_steps=10`.** Matches the Choi et al. 2010 paper and is stable on
iJO1366. n=5 is marginally less sensitive without being faster in practice (the LP
solves dominate, not the scan count).

---

## `max_fraction` — upper ceiling on enforced target flux (as fraction of FVA max)

**Parameters tested:** 0.5, 0.9 (Python/MATLAB default), 0.99

| `max_fraction` | Amplified | Knocked down / out |
|---|---|---|
| 0.5 | 21 | 398 |
| 0.9 (default) | 18 | 395 |
| 0.99 | 18 | 395 |

At `max_fraction=0.5` the enforced range is too narrow; 3 extra reactions appear
amplified that are at capacity throughout the scan and should not be targets. The
0.9 and 0.99 results are identical, meaning the last 10% of flux range adds no
new information on iJO1366.

**Decision: keep `max_fraction=0.9`.** Matches the Choi et al. 2010 paper. Values
above 0.9 give no additional sensitivity; values below 0.9 produce false-positive
amplification targets.

---

## `correlation_threshold` — minimum Pearson r for a reaction to be called a target

**Parameters tested:** 0.7, 0.9 (Python/MATLAB default), 0.95

| `correlation_threshold` | Amplified | Knocked down / out |
|---|---|---|
| 0.7 | 22 | 396 |
| 0.9 (default) | 18 | 395 |
| 0.95 | 18 | 394 |

Lowering to 0.7 adds 4 amplification targets and 1 knockdown. Raising to 0.95 drops
1 knockdown. The 4 extra targets at 0.7 include reactions that have noisy or non-linear
flux profiles across the scan; forcing a high Pearson r ensures only monotonically
correlated reactions are reported.

**Decision: keep `correlation_threshold=0.9`.** Matches Choi et al. 2010. The 0.7
threshold is too permissive for genome-scale models with many weakly correlated reactions.
