"""Tests for set_variance_bounds (the var mode of setParam)."""
import cobra
import pytest

from raven_toolbox.manipulation import add_reactions_from_equations, set_variance_bounds


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


def test_band_positive(model):
    set_variance_bounds(model, "R1", 100, 5)  # 97.5 .. 102.5
    lb, ub = model.reactions.get_by_id("R1").bounds
    assert lb == pytest.approx(97.5)
    assert ub == pytest.approx(102.5)


def test_band_negative_is_ordered(model):
    set_variance_bounds(model, "R1", -100, 5)
    lb, ub = model.reactions.get_by_id("R1").bounds
    assert lb == pytest.approx(-102.5)
    assert ub == pytest.approx(-97.5)
    assert lb <= ub


def test_broadcast_scalar(model):
    set_variance_bounds(model, ["R1", "R2"], 50, 10)
    for rid in ("R1", "R2"):
        lb, ub = model.reactions.get_by_id(rid).bounds
        assert lb == pytest.approx(47.5)
        assert ub == pytest.approx(52.5)


def test_per_reaction_values(model):
    set_variance_bounds(model, ["R1", "R2"], [100, 200], 0)
    assert model.reactions.get_by_id("R1").bounds == pytest.approx((100, 100))
    assert model.reactions.get_by_id("R2").bounds == pytest.approx((200, 200))


def test_length_mismatch_raises(model):
    with pytest.raises(ValueError, match="to match the reactions"):
        set_variance_bounds(model, ["R1", "R2"], [1, 2, 3], 5)


def test_unknown_reaction_raises(model):
    with pytest.raises(ValueError, match="not found"):
        set_variance_bounds(model, "NOPE", 1, 5)
