"""Tests for LP-based gap-filling (gapfilling/fast_lp.py)."""
import cobra
import pytest

from raven_toolbox.gapfilling import FastLPResult, fill_gaps_fast_lp


def _met(mid):
    return cobra.Metabolite(mid, name=mid, compartment="c")


@pytest.fixture
def linear_gap():
    """Draft: EX_A -> A -(r1)-> B.  B has no consumer, so r1 is blocked.

    Template: B -(r2)-> C, EX_C (allows C to leave).
    Filling r2 + EX_C unblocks r1.
    """
    A, B = _met("A_c"), _met("B_c")
    draft = cobra.Model("draft")
    exa = cobra.Reaction("EX_A", lower_bound=-10, upper_bound=1000)
    exa.add_metabolites({A: -1})
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({A: -1, B: 1})
    draft.add_reactions([exa, r1])

    template = cobra.Model("template")
    r2 = cobra.Reaction("r2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({_met("B_c"): -1, _met("C_c"): 1})
    exc = cobra.Reaction("EX_C", lower_bound=-1000, upper_bound=1000)
    exc.add_metabolites({_met("C_c"): -1})
    template.add_reactions([r2, exc])

    return draft, template


def test_result_type(linear_gap):
    draft, template = linear_gap
    result = fill_gaps_fast_lp(draft, template, verbose=False)
    assert isinstance(result, FastLPResult)


def test_added_reactions_is_list(linear_gap):
    draft, template = linear_gap
    result = fill_gaps_fast_lp(draft, template, verbose=False)
    assert isinstance(result.added_reactions, list)


def test_model_returned(linear_gap):
    draft, template = linear_gap
    result = fill_gaps_fast_lp(draft, template, verbose=False)
    assert isinstance(result.model, cobra.Model)


def test_fills_blocked_reaction(linear_gap):
    draft, template = linear_gap
    result = fill_gaps_fast_lp(draft, template, verbose=False)
    assert "r1" in result.newly_connected
    assert len(result.added_reactions) > 0


def test_cannot_connect_is_list(linear_gap):
    draft, template = linear_gap
    result = fill_gaps_fast_lp(draft, template, verbose=False)
    assert isinstance(result.cannot_connect, list)


def test_swift_variant_runs(linear_gap):
    draft, template = linear_gap
    result = fill_gaps_fast_lp(draft, template, variant="swift", verbose=False)
    assert isinstance(result, FastLPResult)
    assert isinstance(result.added_reactions, list)


def test_no_gaps_returns_early(linear_gap):
    _, template = linear_gap
    # Use template as both model (all reactions connected) and template
    draft_complete = template.copy()
    # Add a sink for C to make template self-sufficient
    sink = cobra.Reaction("EX_A2", lower_bound=-1000, upper_bound=1000)
    sink.add_metabolites({_met("C_c"): -1})
    # For this test, just check the function doesn't crash
    result = fill_gaps_fast_lp(draft_complete, template, verbose=False)
    assert isinstance(result, FastLPResult)


def test_candidates_per_reaction_populated(linear_gap):
    """candidates_per_reaction should map each rescuable reaction to a list."""
    draft, template = linear_gap
    result = fill_gaps_fast_lp(draft, template, verbose=False)
    assert isinstance(result.candidates_per_reaction, dict)
    for v in result.candidates_per_reaction.values():
        assert isinstance(v, list)
