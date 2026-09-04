"""Tests for follow_changed (analysis/follow_changed.py, followChanged port)."""
import cobra
import pandas as pd
import pytest

from raven_toolbox.analysis import follow_changed, print_changed_fluxes


@pytest.fixture
def model():
    """A <-> B <-> C, plus a side branch D touched only by r3.

    r1: A -> B, r2: B -> C, r3: B -> D. metabolite_list filters exercise
    "B" (touches all three) vs "D" (touches only r3).
    """
    m = cobra.Model("t")
    a, b, c, d = (cobra.Metabolite(x, name=x, compartment="c") for x in "ABCD")
    m.add_metabolites([a, b, c, d])
    r1 = cobra.Reaction("r1", lower_bound=-1000, upper_bound=1000)
    r1.add_metabolites({a: -1, b: 1})
    r2 = cobra.Reaction("r2", lower_bound=-1000, upper_bound=1000)
    r2.add_metabolites({b: -1, c: 1})
    r3 = cobra.Reaction("r3", lower_bound=-1000, upper_bound=1000)
    r3.add_metabolites({b: -1, d: 1})
    m.add_reactions([r1, r2, r3])
    return m


def test_identical_fluxes_are_not_changed(model):
    fluxes = pd.Series({"r1": 10.0, "r2": 10.0, "r3": 10.0})
    result = follow_changed(model, fluxes, fluxes, cutoff_flux=1, cutoff_diff=0.5, cutoff_change=5)
    assert result.changed == []


def test_small_absolute_difference_excluded_by_cutoff_diff(model):
    fluxes_a = pd.Series({"r1": 10.0, "r2": 1.0, "r3": 100.0})
    fluxes_b = pd.Series({"r1": 10.0, "r2": 1.4, "r3": 102.0})
    result = follow_changed(model, fluxes_a, fluxes_b, cutoff_flux=1, cutoff_diff=0.5, cutoff_change=5)
    # r2: diff=0.4 < 0.5 -> excluded even though the ratio (1.4) is a big change.
    assert "r2" not in [c.reaction for c in result.changed]


def test_small_percent_change_excluded_by_cutoff_change(model):
    fluxes_a = pd.Series({"r1": 100.0, "r2": 1.0, "r3": 1.0})
    fluxes_b = pd.Series({"r1": 102.0, "r2": 1.0, "r3": 1.0})
    result = follow_changed(model, fluxes_a, fluxes_b, cutoff_flux=1, cutoff_diff=0.5, cutoff_change=5)
    # r1: diff=2 (>=0.5) and both clear cutoff_flux, but ratio 1.02 is within
    # [1/1.05, 1.05) -- too small a percent change, excluded.
    assert [c.reaction for c in result.changed] == []


def test_significant_change_included(model):
    fluxes_a = pd.Series({"r1": 10.0, "r2": 1.0, "r3": 1.0})
    fluxes_b = pd.Series({"r1": 12.0, "r2": 1.0, "r3": 1.0})
    result = follow_changed(model, fluxes_a, fluxes_b, cutoff_flux=1, cutoff_diff=0.5, cutoff_change=5)
    assert [c.reaction for c in result.changed] == ["r1"]
    c = result.changed[0]
    assert c.flux == 10.0
    assert c.reference_flux == 12.0
    assert c.difference == pytest.approx(-2.0)


def test_zero_flux_a_included_when_b_clears_cutoff(model):
    fluxes_a = pd.Series({"r1": 0.0, "r2": 1.0, "r3": 1.0})
    fluxes_b = pd.Series({"r1": 5.0, "r2": 1.0, "r3": 1.0})
    result = follow_changed(model, fluxes_a, fluxes_b, cutoff_flux=1, cutoff_diff=0.5, cutoff_change=5)
    assert [c.reaction for c in result.changed] == ["r1"]


def test_zero_flux_a_excluded_when_b_below_cutoff(model):
    fluxes_a = pd.Series({"r1": 0.0, "r2": 1.0, "r3": 1.0})
    fluxes_b = pd.Series({"r1": 0.5, "r2": 1.0, "r3": 1.0})
    result = follow_changed(model, fluxes_a, fluxes_b, cutoff_flux=1, cutoff_diff=0.5, cutoff_change=5)
    assert [c.reaction for c in result.changed] == []


def test_sign_flip_counts_as_changed(model):
    fluxes_a = pd.Series({"r1": 5.0, "r2": 1.0, "r3": 1.0})
    fluxes_b = pd.Series({"r1": -5.0, "r2": 1.0, "r3": 1.0})
    result = follow_changed(model, fluxes_a, fluxes_b, cutoff_flux=1, cutoff_diff=0.5, cutoff_change=5)
    assert [c.reaction for c in result.changed] == ["r1"]


def test_metabolite_list_restricts_to_touching_reactions(model):
    fluxes_a = pd.Series({"r1": 10.0, "r2": 10.0, "r3": 10.0})
    fluxes_b = pd.Series({"r1": 20.0, "r2": 20.0, "r3": 1.0})
    result = follow_changed(
        model, fluxes_a, fluxes_b, cutoff_flux=1, cutoff_diff=0.5, cutoff_change=5,
        metabolite_list=["D"],
    )
    # Only r3 touches D; r1/r2 changed too but are out of scope.
    assert [c.reaction for c in result.changed] == ["r3"]


def test_missing_metabolite_reported_not_raised(model):
    fluxes = pd.Series({"r1": 10.0, "r2": 10.0, "r3": 10.0})
    result = follow_changed(
        model, fluxes, fluxes, metabolite_list=["nonexistent metabolite"],
    )
    assert result.missing_metabolites == ["nonexistent metabolite"]
    assert result.changed == []


def test_metabolite_name_match_is_case_insensitive(model):
    fluxes_a = pd.Series({"r1": 10.0, "r2": 1.0, "r3": 1.0})
    fluxes_b = pd.Series({"r1": 20.0, "r2": 1.0, "r3": 1.0})
    result = follow_changed(
        model, fluxes_a, fluxes_b, cutoff_flux=1, cutoff_diff=0.5, cutoff_change=5,
        metabolite_list=["a"],
    )
    assert [c.reaction for c in result.changed] == ["r1"]


def test_print_changed_fluxes_reports_missing_and_returns_result(model):
    fluxes = pd.Series({"r1": 10.0, "r2": 20.0, "r3": 10.0})
    printed = []
    result = print_changed_fluxes(
        model, fluxes, pd.Series({"r1": 10.0, "r2": 10.0, "r3": 10.0}),
        cutoff_flux=1, cutoff_diff=0.5, cutoff_change=5,
        metabolite_list=["nope"], print_fn=printed.append,
    )
    assert any("Could not find" in line for line in printed)
    assert [c.reaction for c in result.changed] == []  # metabolite filter excludes everything
