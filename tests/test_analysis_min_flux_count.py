"""Tests for minimum-cardinality flux search (analysis/min_flux_count.py, getMinNrFluxes port)."""
import cobra
import pytest

from raven_toolbox.analysis import MinNrFluxesResult, get_min_nr_fluxes


@pytest.fixture
def two_source_model():
    """Two equally-capable sources feed X; a fixed demand of 5 must be met.

    Splitting flux across both sources satisfies mass balance just as well
    as using one alone, but costs an extra active reaction -- a minimal
    solution should always pick exactly one source, never both.
    """
    m = cobra.Model("t")
    x = cobra.Metabolite("X", compartment="c")
    m.add_metabolites([x])
    source1 = cobra.Reaction("source1", lower_bound=0, upper_bound=1000)
    source1.add_metabolites({x: 1})
    source2 = cobra.Reaction("source2", lower_bound=0, upper_bound=1000)
    source2.add_metabolites({x: 1})
    demand = cobra.Reaction("demand", lower_bound=5, upper_bound=5)
    demand.add_metabolites({x: -1})
    m.add_reactions([source1, source2, demand])
    return m


def test_picks_one_source_not_both(two_source_model):
    res = get_min_nr_fluxes(two_source_model, ["source1", "source2", "demand"])
    assert isinstance(res, MinNrFluxesResult)
    assert res.status == "optimal"
    assert "demand" in res.active
    used_sources = {r for r in ("source1", "source2") if r in res.active}
    assert len(used_sources) == 1


def test_fluxes_cover_every_model_reaction(two_source_model):
    res = get_min_nr_fluxes(two_source_model, ["source1", "source2", "demand"])
    assert set(res.fluxes.index) == {"source1", "source2", "demand"}
    assert res.fluxes["demand"] == pytest.approx(5.0)


def test_scores_break_the_tie_deterministically(two_source_model):
    # source2 heavily discouraged (very negative score) relative to source1's default.
    res = get_min_nr_fluxes(
        two_source_model, ["source1", "source2", "demand"], scores=[-1, -100, -1],
    )
    assert res.active == ["source1", "demand"] or set(res.active) == {"source1", "demand"}
    assert "source2" not in res.active


def test_scores_length_mismatch_raises(two_source_model):
    with pytest.raises(ValueError, match="same length"):
        get_min_nr_fluxes(two_source_model, ["source1", "source2"], scores=[-1])


def test_all_nonnegative_scores_no_crash(two_source_model):
    """MATLAB's max([]) on an all-non-negative scores input silently
    misaligns the array (a latent bug, see module docstring); this port
    clamps to 0 instead, which just makes the search unweighted rather
    than raising or corrupting anything."""
    res = get_min_nr_fluxes(
        two_source_model, ["source1", "source2", "demand"], scores=[0, 1, 2],
    )
    assert res.status == "optimal"


def test_infeasible_returns_empty_result():
    m = cobra.Model("infeasible")
    x = cobra.Metabolite("X", compartment="c")
    m.add_metabolites([x])
    source = cobra.Reaction("source", lower_bound=0, upper_bound=0)
    source.add_metabolites({x: 1})
    demand = cobra.Reaction("demand", lower_bound=5, upper_bound=5)
    demand.add_metabolites({x: -1})
    m.add_reactions([source, demand])

    res = get_min_nr_fluxes(m, ["source", "demand"])
    assert res.status != "optimal"
    assert res.active == []
    assert res.fluxes.empty


def test_default_to_minimize_is_all_reactions(two_source_model):
    res = get_min_nr_fluxes(two_source_model)
    assert res.status == "optimal"
    assert set(res.fluxes.index) == {"source1", "source2", "demand"}
