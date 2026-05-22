"""Tests for set_parameters (setParam port)."""
import cobra
import pytest

from ravengem.manipulation import add_reactions_from_equations, set_parameters


@pytest.fixture
def model():
    m = cobra.Model("t")
    m.add_metabolites(
        [cobra.Metabolite("a_c", compartment="c"), cobra.Metabolite("b_c", compartment="c")]
    )
    add_reactions_from_equations(
        m,
        [
            {"id": "R1", "equation": "a_c <=> b_c"},
            {"id": "R2", "equation": "a_c <=> b_c"},
        ],
    )
    return m


def test_set_lb_ub(model):
    set_parameters(model, ["R1"], lb=-50, ub=100)
    assert model.reactions.get_by_id("R1").bounds == (-50, 100)


def test_broadcast_scalar(model):
    set_parameters(model, ["R1", "R2"], lb=0)
    assert model.reactions.get_by_id("R1").lower_bound == 0
    assert model.reactions.get_by_id("R2").lower_bound == 0


def test_per_reaction_values(model):
    set_parameters(model, ["R1", "R2"], ub=[10, 20])
    assert model.reactions.get_by_id("R1").upper_bound == 10
    assert model.reactions.get_by_id("R2").upper_bound == 20


def test_eq(model):
    set_parameters(model, "R1", eq=5)
    assert model.reactions.get_by_id("R1").bounds == (5, 5)


def test_var_band_positive(model):
    # value 100, percent 5 -> 97.5 .. 102.5
    set_parameters(model, "R1", var=(100, 5))
    lb, ub = model.reactions.get_by_id("R1").bounds
    assert lb == pytest.approx(97.5)
    assert ub == pytest.approx(102.5)


def test_var_band_negative(model):
    # negative value swaps so lb <= ub
    set_parameters(model, "R1", var=(-100, 5))
    lb, ub = model.reactions.get_by_id("R1").bounds
    assert lb == pytest.approx(-102.5)
    assert ub == pytest.approx(-97.5)


def test_reset(model):
    set_parameters(model, "R1", lb=-5, ub=5)
    set_parameters(model, "R1", reset=True)
    cfg = cobra.Configuration()
    assert model.reactions.get_by_id("R1").bounds == (cfg.lower_bound, cfg.upper_bound)


def test_objective_replaces(model):
    set_parameters(model, "R1", objective=1)
    set_parameters(model, "R2", objective=1)
    # objective replaced, so only R2 carries it
    assert model.reactions.get_by_id("R1").objective_coefficient == 0
    assert model.reactions.get_by_id("R2").objective_coefficient == 1


def test_invalid_bounds_raise(model):
    # cobra itself rejects lb > ub on assignment; we surface a ValueError either way.
    with pytest.raises(ValueError):
        set_parameters(model, "R1", lb=10, ub=1)


def test_unknown_reaction_raises(model):
    with pytest.raises(ValueError, match="not found"):
        set_parameters(model, "NOPE", lb=0)
