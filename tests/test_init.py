"""Tests for the INIT MILP (init/init.py, Phase 4c)."""
import cobra
import pytest

from ravengem.init import InitResult, run_init


def _met(mid):
    return cobra.Metabolite(mid, name=mid[:-2] if mid.endswith("_c") else mid, compartment="c")


@pytest.fixture
def model():
    """EX_A -> A -(r1)-> B -(r2)-> C -(r3)-> D, with A uptake and excretion allowed.

    r1, r2 are good (positive score); r3 is bad (negative score).
    """
    m = cobra.Model("net")
    A, B, C, D = _met("A_c"), _met("B_c"), _met("C_c"), _met("D_c")
    m.add_metabolites([A, B, C, D])
    exa = cobra.Reaction("EX_A", lower_bound=-1000, upper_bound=1000)
    exa.add_metabolites({A: -1})  # negative flux = uptake of A
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({A: -1, B: 1})
    r2 = cobra.Reaction("r2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({B: -1, C: 1})
    r3 = cobra.Reaction("r3", lower_bound=0, upper_bound=1000)
    r3.add_metabolites({C: -1, D: 1})
    m.add_reactions([exa, r1, r2, r3])
    return m


def test_keeps_positive_drops_negative(model):
    scores = {"r1": 1.0, "r2": 1.0, "r3": -1.0}
    res = run_init(model, scores, prod_weight=0.0, allow_excretion=True)
    assert isinstance(res, InitResult)
    kept = {r.id for r in res.model.reactions}
    assert {"r1", "r2"} <= kept  # positive-score, flux-consistent -> kept
    assert "r3" in res.deleted_reactions  # negative score -> removed
    assert "r3" not in kept


def test_negative_scores_emptied_when_no_reward(model):
    # All reactions negative and no production reward -> keep nothing (empty optimum).
    scores = {r.id: -1.0 for r in model.reactions}
    res = run_init(model, scores, prod_weight=0.0, allow_excretion=True)
    assert res.deleted_reactions == sorted(r.id for r in model.reactions)
    assert len(res.model.reactions) == 0


def test_essential_reaction_forced_kept(model):
    # r3 is negative-scored but essential -> must be kept despite the penalty.
    scores = {"r1": 1.0, "r2": 1.0, "r3": -1.0}
    res = run_init(model, scores, essential_rxns=["r3"], prod_weight=0.0, allow_excretion=True)
    kept = {r.id for r in res.model.reactions}
    assert "r3" in kept
    assert "r3" not in res.deleted_reactions


def test_prod_weight_pulls_in_connectivity(model):
    # With everything scored 0, no reward -> empty. With prod_weight>0, producing
    # metabolites is rewarded, so flux-carrying reactions are pulled in.
    zero = {r.id: 0.0 for r in model.reactions}
    empty = run_init(model, zero, prod_weight=0.0, allow_excretion=True)
    assert len(empty.model.reactions) == 0
    pulled = run_init(model, zero, prod_weight=0.5, allow_excretion=True)
    assert len(pulled.model.reactions) > 0


def test_present_mets_reports_producibility(model):
    scores = {"r1": 1.0, "r2": 1.0}
    res = run_init(
        model, scores, present_mets=["C", "Z"], prod_weight=0.0, allow_excretion=True
    )
    assert res.met_production["C"] is True   # A->B->C is producible
    assert res.met_production["Z"] is False  # not in the model


def test_objective_returned(model):
    res = run_init(model, {"r1": 1.0, "r2": 1.0, "r3": -1.0}, prod_weight=0.0, allow_excretion=True)
    assert res.objective == pytest.approx(2.0)  # kept r1(+1) + r2(+1), dropped r3


def test_reversible_essential_keeps_productive_path():
    """A reversible essential reaction must not be forced into a phantom fwd+rev loop.

    SRC -> a, R: a <=> b (reversible, essential), SNK: b ->. Forcing R essential
    should keep the productive path SRC->R->SNK, not delete SRC/SNK and leave R
    self-looping (the bug from forcing eps flux through both split directions).
    """
    import cobra

    m = cobra.Model("revess")
    a, b = (cobra.Metabolite(x, compartment="c") for x in "ab")
    m.add_metabolites([a, b])
    src = cobra.Reaction("SRC", lower_bound=0, upper_bound=1000)
    src.add_metabolites({a: 1})
    r = cobra.Reaction("R", lower_bound=-1000, upper_bound=1000)
    r.add_metabolites({a: -1, b: 1})
    snk = cobra.Reaction("SNK", lower_bound=0, upper_bound=1000)
    snk.add_metabolites({b: -1})
    m.add_reactions([src, r, snk])
    m.objective = "SNK"

    res = run_init(m, {"SRC": -1.0, "SNK": -1.0}, essential_rxns=["R"], prod_weight=0.0)
    kept = {rxn.id for rxn in res.model.reactions}
    assert "R" in kept
    # The productive path must be kept (SRC feeds R, SNK drains it); R can't self-loop.
    assert {"SRC", "SNK"} <= kept
    assert res.model.slim_optimize() > 1e-6  # the kept model actually carries flux
