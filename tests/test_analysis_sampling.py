"""Tests for random-objective flux sampling (analysis/sampling.py)."""
import cobra
import numpy as np
import pytest

from ravengem.analysis import (
    RandomSamplingResult,
    find_good_reactions,
    random_sampling,
)


@pytest.fixture
def model():
    """S uptake -> A -> {B export, C export}, plus a thermodynamically infeasible loop.

    sup -> A; A->B (v_b) and A->C (v_c); B,C exported. r_f/r_r form a closed cycle
    (X<->Y both directions, no in/out) that can spin arbitrarily — a loop whose
    reactions must be excluded from the random objectives.
    """
    m = cobra.Model("toy")
    A, B, C, X, Y = (cobra.Metabolite(x, compartment="c") for x in "ABCXY")
    m.add_metabolites([A, B, C, X, Y])

    def rxn(rid, mets, lb=0, ub=1000):
        r = cobra.Reaction(rid, lower_bound=lb, upper_bound=ub)
        r.add_metabolites(mets)
        return r

    rxns = [
        rxn("sup", {A: 1}, ub=10),       # substrate supply
        rxn("v_b", {A: -1, B: 1}),       # A -> B
        rxn("v_c", {A: -1, C: 1}),       # A -> C
        rxn("EX_B", {B: -1}),            # export B
        rxn("EX_C", {C: -1}),            # export C
        rxn("r_f", {X: -1, Y: 1}, lb=-1000),  # X <-> Y  ┐ closed loop
        rxn("r_r", {Y: -1, X: 1}, lb=-1000),  # Y <-> X  ┘ (no source/sink for X,Y)
    ]
    m.add_reactions(rxns)
    m.objective = "EX_B"
    return m


def test_good_reactions_excludes_loop(model):
    good = find_good_reactions(model)
    # The closed X<->Y cycle can spin to the 1000 bound -> excluded.
    assert "r_f" not in good and "r_r" not in good
    # Real flux-carrying reactions are kept.
    assert {"sup", "v_b", "EX_B"} <= set(good)


def test_returns_result_shape(model):
    res = random_sampling(model, n_samples=20, seed=1)
    assert isinstance(res, RandomSamplingResult)
    assert res.samples.shape == (20, len(model.reactions))
    assert list(res.samples.columns) == [r.id for r in model.reactions]
    assert "r_f" not in res.good_reactions


def test_samples_are_steady_state(model):
    """Every sample must satisfy S·v = 0 (mass balance)."""
    res = random_sampling(model, n_samples=15, seed=2)
    s_matrix = cobra.util.create_stoichiometric_matrix(model)
    ids = [r.id for r in model.reactions]
    for _, row in res.samples.iterrows():
        residual = s_matrix @ row[ids].to_numpy()
        assert np.allclose(residual, 0, atol=1e-6)


def test_samples_respect_bounds(model):
    res = random_sampling(model, n_samples=15, seed=3)
    for r in model.reactions:
        col = res.samples[r.id].to_numpy()
        assert (col >= r.lower_bound - 1e-6).all()
        assert (col <= r.upper_bound + 1e-6).all()


def test_seed_is_reproducible(model):
    a = random_sampling(model, n_samples=10, seed=42).samples
    b = random_sampling(model, n_samples=10, seed=42).samples
    assert np.allclose(a.to_numpy(), b.to_numpy())


def test_good_reactions_reused(model):
    """Passing good_reactions back in reproduces the FVA-derived set without recomputing."""
    good = find_good_reactions(model)
    res = random_sampling(model, n_samples=5, good_reactions=good, seed=0)
    assert res.good_reactions == good


def test_min_flux_runs(model):
    res = random_sampling(model, n_samples=8, min_flux=True, seed=5)
    assert res.samples.shape == (8, len(model.reactions))


def test_diverse_samples(model):
    """Random objectives should explore different states, not a single FBA optimum."""
    res = random_sampling(model, n_samples=40, seed=7)
    # The branch split A->B vs A->C should vary across samples.
    assert res.samples["v_b"].std() > 1e-6
    assert res.samples["v_c"].std() > 1e-6


def test_rejects_bad_n_samples(model):
    with pytest.raises(ValueError, match="n_samples"):
        random_sampling(model, n_samples=0)


def test_too_few_good_reactions(model):
    with pytest.raises(ValueError, match="usable reactions"):
        random_sampling(model, n_samples=5, good_reactions=["sup"], n_objectives=2)


def test_good_reactions_keeps_reactions_at_default_bound():
    """A legitimate reaction reaching the model's 1000 bound is not dropped as a loop.

    Regression: the old loop_bound>=1000 test wrongly excluded any reaction that
    reaches the default bound. Loopless FVA keeps it (real flux) and still drops a
    closed loop.
    """
    m = cobra.Model("b")
    a, b = (cobra.Metabolite(x, compartment="c") for x in "ab")
    m.add_metabolites([a, b])
    sup = cobra.Reaction("sup", lower_bound=0, upper_bound=1000)  # uptake to the 1000 cap
    sup.add_metabolites({a: 1})
    conv = cobra.Reaction("conv", lower_bound=0, upper_bound=1000)
    conv.add_metabolites({a: -1, b: 1})
    ex = cobra.Reaction("EX_b", lower_bound=0, upper_bound=1000)
    ex.add_metabolites({b: -1})
    m.add_reactions([sup, conv, ex])
    m.objective = "EX_b"
    good = find_good_reactions(m)
    assert {"sup", "conv", "EX_b"} <= set(good)  # all reach 1000 but are real, not loops
