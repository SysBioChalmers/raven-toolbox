"""Tests for CHRR / ACHR flux sampling (analysis/flux_sampling.py).

The maximum-volume ellipsoid (MVE) solver is the numerical crux of CHRR, so it is
validated directly against analytic cases (box -> unit ball, scaled/sheared box,
simplex -> incircle) before the end-to-end sampler is exercised on cobra models.
"""
import cobra
import numpy as np
import pytest

from raven_toolbox.analysis import (
    FluxSamplingResult,
    max_volume_ellipsoid,
    random_sampling,
)


# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #
@pytest.fixture
def branched_model():
    """A -> {B, C} with bounded uptake; a small bounded full-dimensional polytope.

    sup: -> A (ub 10); v_b: A -> B; v_c: A -> C; EX_B, EX_C export. No loops, so the
    flux polytope is a bounded 2-D slab (v_b + v_c = uptake), good for sampling.
    """
    m = cobra.Model("branched")
    A, B, C = (cobra.Metabolite(x, compartment="c") for x in "ABC")
    m.add_metabolites([A, B, C])

    def rxn(rid, mets, lb=0.0, ub=1000.0):
        r = cobra.Reaction(rid, lower_bound=lb, upper_bound=ub)
        r.add_metabolites(mets)
        return r

    m.add_reactions([
        rxn("sup", {A: 1}, ub=10.0),
        rxn("v_b", {A: -1, B: 1}),
        rxn("v_c", {A: -1, C: 1}),
        rxn("EX_B", {B: -1}),
        rxn("EX_C", {C: -1}),
    ])
    m.objective = "EX_B"
    return m


# --------------------------------------------------------------------------- #
# MVE solver — analytic validation
# --------------------------------------------------------------------------- #
def _box(n, half=1.0):
    """A z <= b for the box {-half <= z_i <= half}."""
    A = np.vstack([np.eye(n), -np.eye(n)])
    b = np.full(2 * n, half)
    return A, b


def test_mve_unit_box_is_unit_ball():
    A, b = _box(3, 1.0)
    x, E, converged = max_volume_ellipsoid(A, b, np.zeros(3))
    assert converged
    assert np.allclose(x, 0.0, atol=1e-5)
    assert np.allclose(E @ E.T, np.eye(3), atol=1e-4)


def test_mve_scaled_box_is_diagonal():
    # box {-r_i <= z_i <= r_i}: MVE is axis-aligned ellipsoid E2 = diag(r_i^2).
    r = np.array([2.0, 0.5, 3.0])
    A = np.vstack([np.eye(3), -np.eye(3)])
    b = np.concatenate([r, r])
    x, E, converged = max_volume_ellipsoid(A, b, np.zeros(3))
    assert converged
    assert np.allclose(x, 0.0, atol=1e-5)
    assert np.allclose(E @ E.T, np.diag(r**2), atol=1e-4)


def test_mve_sheared_box_recovers_inverse_map():
    # {-1 <= M z <= 1}: MVE matrix E2 = M^{-1} M^{-T}.
    M = np.array([[1.0, 0.5], [0.0, 1.0]])
    A = np.vstack([M, -M])
    b = np.ones(4)
    x, E, converged = max_volume_ellipsoid(A, b, np.zeros(2))
    assert converged
    Minv = np.linalg.inv(M)
    assert np.allclose(x, 0.0, atol=1e-5)
    assert np.allclose(E @ E.T, Minv @ Minv.T, atol=1e-4)


def test_mve_offcenter_box_centers_correctly():
    # box [0,4] x [1,2]: centre (2, 1.5), E2 = diag(2^2, 0.5^2).
    lb = np.array([0.0, 1.0])
    ub = np.array([4.0, 2.0])
    A = np.vstack([np.eye(2), -np.eye(2)])
    b = np.concatenate([ub, -lb])
    x, E, converged = max_volume_ellipsoid(A, b)
    assert converged
    assert np.allclose(x, [2.0, 1.5], atol=1e-4)
    assert np.allclose(E @ E.T, np.diag([4.0, 0.25]), atol=1e-4)


def test_mve_triangle_is_steiner_inellipse():
    # The maximum-area inscribed ellipse of a triangle is the Steiner inellipse,
    # centred at the centroid and tangent at the side midpoints -- NOT the incircle.
    # For the standard simplex {z >= 0, z1 + z2 <= 1} (vertices (0,0),(1,0),(0,1)):
    #   centroid = (1/3, 1/3);  E2^{-1} = [[12,6],[6,12]]  =>  E2 = [[1/9,-1/18],
    #   [-1/18,1/9]] (det E2 = 1/108, matching the Steiner-inellipse area pi/(6*sqrt3)).
    # This confirms the solver maximises volume (finds the off-diagonal Steiner
    # shape), not merely the largest inscribed ball.
    A = np.array([[-1.0, 0.0], [0.0, -1.0], [1.0, 1.0]])
    b = np.array([0.0, 0.0, 1.0])
    x, E, converged = max_volume_ellipsoid(A, b, np.array([0.25, 0.25]))
    E2_expected = np.array([[1.0 / 9, -1.0 / 18], [-1.0 / 18, 1.0 / 9]])
    assert converged
    assert np.allclose(x, [1.0 / 3, 1.0 / 3], atol=1e-4)
    assert np.allclose(E @ E.T, E2_expected, atol=1e-4)
    assert abs(np.linalg.det(E @ E.T) - 1.0 / 108) < 1e-5


def test_mve_rejects_exterior_x0():
    A, b = _box(2, 1.0)
    with pytest.raises(ValueError, match="interior"):
        max_volume_ellipsoid(A, b, np.array([5.0, 0.0]))


def test_chebyshev_center_used_when_x0_none():
    A, b = _box(2, 1.0)
    x, E, converged = max_volume_ellipsoid(A, b)  # x0=None -> Chebyshev centre
    assert converged
    assert np.allclose(x, 0.0, atol=1e-5)


# --------------------------------------------------------------------------- #
# End-to-end CHRR
# --------------------------------------------------------------------------- #
def test_chrr_returns_result_shape(branched_model):
    res = random_sampling(branched_model, method="chrr", n_samples=50, seed=1)
    assert isinstance(res, FluxSamplingResult)
    assert res.method == "chrr"
    assert res.samples.shape == (50, len(branched_model.reactions))
    assert list(res.samples.columns) == [r.id for r in branched_model.reactions]


def test_chrr_samples_are_steady_state(branched_model):
    res = random_sampling(branched_model, method="chrr", n_samples=40, seed=2)
    Sm = cobra.util.create_stoichiometric_matrix(branched_model)
    ids = [r.id for r in branched_model.reactions]
    resid = Sm @ res.samples[ids].to_numpy().T
    assert np.allclose(resid, 0.0, atol=1e-6)


def test_chrr_samples_respect_bounds(branched_model):
    res = random_sampling(branched_model, method="chrr", n_samples=60, seed=3)
    for r in branched_model.reactions:
        col = res.samples[r.id].to_numpy()
        assert (col >= r.lower_bound - 1e-6).all()
        assert (col <= r.upper_bound + 1e-6).all()


def test_chrr_reproducible_with_seed(branched_model):
    a = random_sampling(branched_model, method="chrr", n_samples=30, seed=42).samples
    b = random_sampling(branched_model, method="chrr", n_samples=30, seed=42).samples
    assert np.allclose(a.to_numpy(), b.to_numpy())


def test_chrr_different_seeds_differ(branched_model):
    a = random_sampling(branched_model, method="chrr", n_samples=30, seed=1).samples
    b = random_sampling(branched_model, method="chrr", n_samples=30, seed=2).samples
    assert not np.allclose(a.to_numpy(), b.to_numpy())


def test_chrr_explores_the_branch(branched_model):
    # 5 reactions, 3 mass balances of rank 3 -> a 2-D flux polytope (free v_b, v_c,
    # with v_b + v_c = uptake <= 10). Both branches must vary across samples.
    res = random_sampling(branched_model, method="chrr", n_samples=200, seed=7)
    assert res.n_dimensions == 2
    assert res.samples["v_b"].std() > 1e-3
    assert res.samples["v_c"].std() > 1e-3


def test_chrr_uniform_marginal_on_a_box():
    """CHRR on an axis-aligned box must give ~uniform marginals.

    Model: two independent uptake->export channels with no coupling, so the flux
    polytope is a 2-D box [0,10] x [0,5]. The recorded marginals should be uniform:
    mean near the midpoint and ~half the mass in each half-interval.
    """
    m = cobra.Model("box")
    A, B = (cobra.Metabolite(x, compartment="c") for x in "AB")
    m.add_metabolites([A, B])
    rs = [
        cobra.Reaction("inA", lower_bound=0, upper_bound=10),
        cobra.Reaction("exA", lower_bound=0, upper_bound=10),
        cobra.Reaction("inB", lower_bound=0, upper_bound=5),
        cobra.Reaction("exB", lower_bound=0, upper_bound=5),
    ]
    rs[0].add_metabolites({A: 1})
    rs[1].add_metabolites({A: -1})
    rs[2].add_metabolites({B: 1})
    rs[3].add_metabolites({B: -1})
    m.add_reactions(rs)

    res = random_sampling(m, method="chrr", n_samples=4000, thinning=40, seed=11)
    a = res.samples["inA"].to_numpy()
    b = res.samples["inB"].to_numpy()
    # Uniform on [0,10] -> mean 5; on [0,5] -> mean 2.5.
    assert abs(a.mean() - 5.0) < 0.5
    assert abs(b.mean() - 2.5) < 0.3
    # Roughly half the samples below the midpoint (uniformity, not just centring).
    assert 0.4 < (a < 5.0).mean() < 0.6
    assert 0.4 < (b < 2.5).mean() < 0.6


# --------------------------------------------------------------------------- #
# ACHR wrapper
# --------------------------------------------------------------------------- #
def test_achr_returns_result_shape(branched_model):
    res = random_sampling(branched_model, method="achr", n_samples=50, thinning=20, seed=1)
    assert isinstance(res, FluxSamplingResult)
    assert res.method == "achr"
    assert res.samples.shape == (50, len(branched_model.reactions))
    assert res.n_warmup is not None and res.n_warmup > 0


def test_achr_samples_are_steady_state(branched_model):
    res = random_sampling(branched_model, method="achr", n_samples=40, thinning=20, seed=2)
    Sm = cobra.util.create_stoichiometric_matrix(branched_model)
    ids = [r.id for r in branched_model.reactions]
    resid = Sm @ res.samples[ids].to_numpy().T
    assert np.allclose(resid, 0.0, atol=1e-6)


def test_achr_reproducible_with_seed(branched_model):
    a = random_sampling(branched_model, method="achr", n_samples=25, thinning=20, seed=5).samples
    b = random_sampling(branched_model, method="achr", n_samples=25, thinning=20, seed=5).samples
    assert np.allclose(a.to_numpy(), b.to_numpy())


# --------------------------------------------------------------------------- #
# Dispatcher
# --------------------------------------------------------------------------- #
def test_unknown_method_raises(branched_model):
    with pytest.raises(ValueError, match="Unknown method"):
        random_sampling(branched_model, method="gibbs")


def test_rejects_nonpositive_samples(branched_model):
    with pytest.raises(ValueError, match="n_samples"):
        random_sampling(branched_model, n_samples=0)
