# Markov-chain flux sampling — CHRR and ACHR

This document describes the two Markov-chain Monte Carlo (MCMC) flux samplers added to
both raven-toolbox (`raven_toolbox.analysis.sample_flux_space`) and MATLAB RAVEN
(`sampleFluxSpace`): **CHRR** (Coordinate Hit-and-Run with Rounding) and **ACHR**
(Artificially Centered Hit-and-Run). It explains what each samples, why rounding matters,
and which to use for enzyme-constrained and tissue-specific models.

---

## Contents

1. [What flux sampling computes](#1-what-flux-sampling-computes)
2. [Why rounding matters: the elongated-polytope problem](#2-why-rounding-matters)
3. [ACHR](#3-achr)
4. [CHRR](#4-chrr)
5. [The maximum-volume ellipsoid (MVE) core](#5-the-maximum-volume-ellipsoid-core)
6. [CHRR vs ACHR — when it matters](#6-chrr-vs-achr)
7. [Relationship to `randomSampling` / `random_sampling`](#7-relationship-to-random-sampling)
8. [Implementation map](#8-implementation-map)
9. [Validation](#9-validation)
10. [References](#10-references)

---

## 1. What flux sampling computes

For a model with stoichiometric matrix `S` and bounds `lb <= v <= ub`, the steady-state
feasible set is the **flux polytope**

```
P = { v : S v = 0,  lb <= v <= ub }.
```

Uniform sampling of `P` answers questions that a single FBA optimum cannot: the marginal
distribution and correlation structure of fluxes consistent with the constraints, the
spread of feasible states, and how measured data narrows the feasible space. Because `S v = 0`
removes degrees of freedom, `P` is a *lower-dimensional* polytope embedded in flux space;
both samplers handle this by working in a full-dimensional parametrisation.

These MCMC samplers draw the (near-)uniform distribution of the polytope **interior**, which
is the statistically correct object for the questions above. This is distinct from the
random-objective method (`randomSampling`), which draws polytope **vertices** — see §7.

---

## 2. Why rounding matters: the elongated-polytope problem {#2-why-rounding-matters}

Hit-and-run mixing time scales with the polytope's **aspect ratio** — the ratio of its
longest to shortest axis. For a near-isotropic (ball-like) polytope, hit-and-run mixes in
`O*(d²)` steps (`d` = dimension). For a long, thin slab, an unrounded chain crawls along the
long axes and almost never moves across the thin ones, so it *appears* converged (short
within-axis autocorrelation) while having explored only a sliver of the feasible set.

This is exactly the geometry of the models this work targets:

- **Enzyme-constrained models (ecModels) with proteomics data.** Each measured enzyme
  abundance imposes a tight bound `v_j <= kcat_j · [E_j]`. Thousands of such near-binding
  constraints collapse the polytope to a thin slab concentrated near an enzyme-capacity
  ridge.
- **Models conditioned on measured fluxes.** Fixing or tightly bounding measured fluxes
  collapses the polytope onto a low-dimensional face.

Both produce a badly conditioned `P` where unrounded ACHR mixes poorly. **CHRR removes the
conditioning problem** by first applying an affine transformation that makes the polytope
near-isotropic (§5), after which the mixing time no longer depends on the original aspect
ratio.

---

## 3. ACHR

**Artificially Centered Hit-and-Run** (Kaufman & Smith 1998; the sampler in COBRA Toolbox
and cobrapy's `ACHRSampler`):

1. Generate **warmup points** — one flux vector per reaction maximised and one minimised
   (`2·nRxns` LPs), giving vertices that span `P`.
2. Maintain a running **center** = mean of all points generated so far.
3. From the current point, choose a random warmup point, take the direction *through the
   center*, and sample uniformly along the feasible chord in that direction.
4. Record every `thinning`-th point.

Every search direction is a difference of feasible points, so each step stays exactly on
`S v = 0`; no reprojection is needed for the homogeneous (`b = 0`) steady-state problem.
ACHR does **no rounding** — directions come from the warmup vertices, which on an elongated
polytope cluster near the long axes, so the thin directions are under-explored.

ACHR is lighter than CHRR (no ellipsoid solve) and is a good default on well-conditioned
models.

---

## 4. CHRR

**Coordinate Hit-and-Run with Rounding** (Haraldsdóttir et al. 2017):

1. **Reduce to full dimension.** Compute a nullspace basis `N` of the equality system and a
   particular solution `v0`, parametrising `v = v0 + N x`. Reactions whose flux range is
   ~zero (detected by FVA) are folded into the equality system first, so the reduced
   polytope `{x : A x <= b}` with `A = [N; -N]`, `b = [ub - v0; v0 - lb]` is genuinely
   full-dimensional (a non-degenerate MVE requires this).
2. **Round.** Compute the maximum-volume ellipsoid `{x : x = c + E s, ||s|| <= 1}` inscribed
   in the polytope (§5) and substitute `x = c + E y`. The rounded polytope
   `{y : (A E) y <= b - A c}` contains the unit ball and is contained in a ball of radius `d`
   (John's theorem), i.e. it is well-conditioned regardless of the original aspect ratio.
3. **Walk.** Run **coordinate** hit-and-run on the rounded polytope: pick a random coordinate
   `j`, compute the feasible chord from the maintained slack `s = b_r - A_r y` and the column
   `A_r[:,j]`, and sample uniformly. The slack updates in `O(m)` per step.
4. **Map back.** `v = v0 + N (c + E y)` for each recorded point.

The rounding (step 2) is computed **once** per model; the per-step cost afterwards equals
ACHR's. The investment buys aspect-ratio-independent mixing.

---

## 5. The maximum-volume ellipsoid (MVE) core {#5-the-maximum-volume-ellipsoid-core}

Rounding requires the **maximum-volume inscribed ellipsoid** of `{z : A z <= b}`:

```
maximize   log det E
subject to a_i^T x + || E a_i ||_2 <= b_i   for every constraint i,
```

over the centre `x` and SPD matrix `E`. This is a convex program solved by the **Zhang & Gao
(2003)** primal-dual interior-point method — the same routine COBRA's `chrrSampler` uses.
Each Newton iteration forms the ellipsoid matrix `E2 = (Aᵀ diag(y) A)⁻¹` from the dual
variables `y`, computes residuals for dual feasibility, primal slack, and complementarity,
and takes a fraction-to-boundary step. The returned lower-triangular factor `E` (with
`E Eᵀ = E2`) is the rounding transform.

A subtlety worth noting: the MVE of a triangle is its **Steiner inellipse** (centred at the
centroid, tangent at the side midpoints), *not* its incircle — the largest inscribed *ellipse*
generally has larger area than the largest inscribed *circle*. This is used as an exact
validation case (§9): the solver must recover the off-diagonal Steiner shape, confirming it
maximises volume rather than merely fitting the biggest ball.

---

## 6. CHRR vs ACHR — when it matters {#6-chrr-vs-achr}

| Model class | Polytope geometry | Recommended | Why |
|---|---|---|---|
| ecModel + proteomics | Thin slab, many near-binding enzyme bounds | **CHRR** | Rounding finds the slab geometry; ACHR warmup directions miss the thin axes |
| Flux-measured model | Collapsed onto a low-dimensional face | **CHRR** | MVE finds the effective dimensionality automatically |
| ftINIT tissue-specific model | Pruned, moderate conditioning | CHRR preferred; ACHR acceptable | Human GEM-derived (~1500–3000 rxns); rounding cost modest, mixing safer |
| Small / well-conditioned model | Near-isotropic | ACHR (lighter) or CHRR | Either mixes well; ACHR avoids the ellipsoid solve |

For the workflows that motivated this addition — ecModels constrained with proteomics and
experimentally measured fluxes, and ftINIT-generated tissue-specific models — **CHRR is the
default**. The ecModel case is precisely the elongated-polytope regime where unrounded ACHR
gives chains that look converged but are not.

---

## 7. Relationship to `randomSampling` / `random_sampling` {#7-relationship-to-random-sampling}

RAVEN's existing `randomSampling` (and raven-toolbox's `random_sampling`) implement the
**random-objective** method (Bordel et al. 2010): each sample maximises a small random linear
objective, returning a polytope **vertex**. That is a different statistical object — a robust,
fast spread of diverse *optimal* states, well-behaved on large or tightly constrained models
where MCMC mixing is hard. It does **not** approximate the uniform interior distribution.

Use `sample_flux_space` / `sampleFluxSpace` (CHRR/ACHR) when you need the (near-)uniform flux
distribution; use `random_sampling` / `randomSampling` when you want diverse vertices. They
complement each other.

---

## 8. Implementation map {#8-implementation-map}

| Component | MATLAB RAVEN | raven-toolbox (Python) |
|---|---|---|
| Dispatcher | `sampleFluxSpace.m` | `analysis.sample_flux_space` |
| CHRR | `sampleCHRR.m` | `analysis.flux_sampling._sample_chrr` |
| ACHR | `sampleACHR.m` + `sampleWarmupPoints.m` | wraps `cobra.sampling.ACHRSampler` |
| MVE solver | `sampleMaxVolEllipse.m` | `analysis.max_volume_ellipsoid` |
| Interior point | `sampleChebyshevCenter.m` | `_chebyshev_center` (scipy `linprog`) |

**Does cobrapy already have these?** cobrapy ships a mature **ACHR** (`ACHRSampler`, plus the
parallel `OptGPSampler`), so the Python side **wraps** it rather than reimplementing — both
methods are reachable through the one `sample_flux_space` entry point and return the same
`FluxSamplingResult`. cobrapy has **no CHRR** (no rounding sampler), so the Python CHRR is a
genuine addition. MATLAB RAVEN had neither, so both are implemented there from scratch.

The MATLAB filenames are prefixed `sample*` to avoid clashing with COBRA Toolbox's
`sampleCbModel` / `ACHRSampler` / `CHRRSampler` / `mve_solver` on case-insensitive
filesystems.

The Python `max_volume_ellipsoid` is the numerically validated reference (§9) that the MATLAB
`sampleMaxVolEllipse` mirrors line-for-line.

---

## 9. Validation

The MVE solver — the numerical crux of CHRR — is checked against analytic cases with known
closed-form answers:

- **Unit box** `{-1 <= x_i <= 1}` → unit ball (`E = I`, centre `0`).
- **Scaled box** `{-r_i <= x_i <= r_i}` → axis-aligned ellipsoid `E Eᵀ = diag(r_i²)`.
- **Sheared box** `{-1 <= M x <= 1}` → `E Eᵀ = M⁻¹ M⁻ᵀ` (recovers the linear map).
- **Triangle** (standard simplex) → Steiner inellipse, centre `(1/3, 1/3)`,
  `E Eᵀ = [[1/9, -1/18], [-1/18, 1/9]]`, `det = 1/108` (confirms volume maximisation, §5).

The samplers are then checked end-to-end: every sample satisfies `S v = 0` and the bounds;
chains are reproducible under a fixed seed; and CHRR's marginals on an axis-aligned box are
verified uniform (mean at the midpoint, ~half the mass on each side). See
`tests/test_analysis_flux_sampling.py` (Python) and `tSampling.m` (MATLAB).

---

## 10. References

- Kaufman, D.E. & Smith, R.L. (1998). Direction choice for accelerated convergence in
  hit-and-run sampling. *Operations Research* 46(1):84–95.
- Zhang, Y. & Gao, L. (2003). On numerical solution of the maximum volume ellipsoid problem.
  *SIAM Journal on Optimization* 14(1):53–76.
- Bordel, S., Agren, R. & Nielsen, J. (2010). Sampling the solution space of genome-scale
  metabolic networks reveals transcriptional regulation in key enzymes.
  *PLoS Computational Biology* 6(7):e1000859.
- Haraldsdóttir, H.S., Cousins, B., Thiele, I., Fleming, R.M.T. & Vempala, S. (2017). CHRR:
  coordinate hit-and-run with rounding for uniform sampling of constraint-based models.
  *Bioinformatics* 33(11):1741–1743.
- Megchelenbrink, W., Huynen, M. & Marchiori, E. (2014). optGpSampler: an improved tool for
  uniformly sampling the solution-space of genome-scale metabolic networks.
  *PLoS ONE* 9(2):e86587.
