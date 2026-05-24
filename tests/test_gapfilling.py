"""Tests for template-based gap-filling (gapfilling/fill.py, Phase 4b)."""
import cobra
import pytest

from ravengem.gapfilling import GapFillResult, fill_gaps, gapfill_to_objective


def _met(mid):
    return cobra.Metabolite(mid, name=mid, compartment="c")


@pytest.fixture
def draft_and_template():
    """Draft: EX_A -> A -> B (r1), but B has no consumer, so r1 is blocked.

    Template supplies B -> C (r2) and an exchange for C, which unblocks r1.
    """
    A, B = _met("A_c"), _met("B_c")
    draft = cobra.Model("draft")
    exa = cobra.Reaction("EX_A", lower_bound=-10, upper_bound=1000)
    exa.add_metabolites({A: 1})
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)  # A -> B, irreversible
    r1.add_metabolites({A: -1, B: 1})
    draft.add_reactions([exa, r1])

    template = cobra.Model("template")
    r2 = cobra.Reaction("r2", lower_bound=0, upper_bound=1000)  # B -> C
    r2.add_metabolites({_met("B_c"): -1, _met("C_c"): 1})
    exc = cobra.Reaction("EX_C", lower_bound=-1000, upper_bound=1000)
    exc.add_metabolites({_met("C_c"): -1})
    extra = cobra.Reaction("r_unneeded", lower_bound=0, upper_bound=1000)  # D -> E, irrelevant
    extra.add_metabolites({_met("D_c"): -1, _met("E_c"): 1})
    template.add_reactions([r2, exc, extra])
    return draft, template


# --------------------------------------------------------------------------- #
# Connectivity gap-fill
# --------------------------------------------------------------------------- #
def test_fill_gaps_connects_blocked_reaction(draft_and_template):
    draft, template = draft_and_template
    assert "r1" in cobra.flux_analysis.find_blocked_reactions(draft)  # precondition

    res = fill_gaps(draft, template)
    assert isinstance(res, GapFillResult)
    assert "r1" in res.newly_connected
    assert set(res.added_reactions) == {"r2", "EX_C"}  # both needed to drain B
    assert "r_unneeded" not in res.added_reactions  # irrelevant template rxn not added


def test_fill_gaps_returns_working_model_that_unblocks(draft_and_template):
    draft, template = draft_and_template
    res = fill_gaps(draft, template)
    assert {"r2", "EX_C"} <= {r.id for r in res.model.reactions}
    assert "r1" not in cobra.flux_analysis.find_blocked_reactions(res.model)
    # original draft is untouched
    assert "r2" not in {r.id for r in draft.reactions}


def test_fill_gaps_nothing_to_do_when_unblocked(draft_and_template):
    draft, template = draft_and_template
    # give the draft its own drain so r1 is not blocked
    drain = cobra.Reaction("EX_B", lower_bound=-1000, upper_bound=1000)
    drain.add_metabolites({draft.metabolites.B_c: -1})
    draft.add_reactions([drain])
    res = fill_gaps(draft, template)
    assert res.added_reactions == []
    assert res.newly_connected == []


def test_fill_gaps_scores_prefer_higher_scored_reactions():
    # Two alternative single-reaction drains for B; scores should pick the preferred one.
    A, B = _met("A_c"), _met("B_c")
    draft = cobra.Model("draft")
    exa = cobra.Reaction("EX_A", lower_bound=-10, upper_bound=1000)
    exa.add_metabolites({A: 1})
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({A: -1, B: 1})
    draft.add_reactions([exa, r1])
    template = cobra.Model("t")
    d1 = cobra.Reaction("drain1", lower_bound=-1000, upper_bound=1000)
    d1.add_metabolites({_met("B_c"): -1})
    d2 = cobra.Reaction("drain2", lower_bound=-1000, upper_bound=1000)
    d2.add_metabolites({_met("B_c"): -1})
    template.add_reactions([d1, d2])
    # Scores are penalties (higher = preferred = cheaper to include); only one drain
    # is needed, so the less-penalised drain1 is chosen.
    res = fill_gaps(draft, template, scores={"drain1": -1.0, "drain2": -5.0})
    assert res.added_reactions == ["drain1"]


# --------------------------------------------------------------------------- #
# Targeted (objective) gap-fill
# --------------------------------------------------------------------------- #
def test_gapfill_to_objective_adds_missing_reaction():
    A, B = _met("A_c"), _met("B_c")
    draft = cobra.Model("draft")
    exa = cobra.Reaction("EX_A", lower_bound=-10, upper_bound=1000)
    exa.add_metabolites({A: 1})
    bio = cobra.Reaction("BIO", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({B: -1})
    draft.add_reactions([exa, bio])
    draft.objective = "BIO"  # gap: no A -> B
    assert draft.slim_optimize(error_value=0.0) == 0.0  # infeasible objective

    template = cobra.Model("t")
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({_met("A_c"): -1, _met("B_c"): 1})
    template.add_reactions([r1])

    res = gapfill_to_objective(draft, template, lower_bound=0.1)
    assert res.added_reactions == ["r1"]
    assert res.model.slim_optimize() >= 0.1


def test_gapfill_to_objective_noop_when_already_feasible():
    A, B = _met("A_c"), _met("B_c")
    draft = cobra.Model("draft")
    exa = cobra.Reaction("EX_A", lower_bound=-10, upper_bound=1000)
    exa.add_metabolites({A: 1})
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({A: -1, B: 1})
    bio = cobra.Reaction("BIO", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({B: -1})
    draft.add_reactions([exa, r1, bio])
    draft.objective = "BIO"
    res = gapfill_to_objective(draft, cobra.Model("t"), lower_bound=0.1)
    assert res.added_reactions == []


def test_gapfill_to_objective_infeasible_raises():
    A = _met("A_c")
    draft = cobra.Model("draft")
    bio = cobra.Reaction("BIO", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({A: -1})
    draft.add_reactions([bio])
    draft.objective = "BIO"  # no way to make A
    with pytest.raises(RuntimeError, match="infeasible"):
        gapfill_to_objective(draft, cobra.Model("empty"), lower_bound=0.1)
