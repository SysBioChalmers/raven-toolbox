"""Tests for FSEOF (analysis/fseof.py, Phase 5)."""
import cobra
import pytest

from ravengem.analysis import FSEOFResult, fseof


@pytest.fixture
def model():
    """S -> I, then I branches to product P (via v2) or biomass B (via v3).

    Enforcing product export (EX_P) should amplify the product branch (v1, v2) and
    suppress the biomass branch (v3), which competes for the shared intermediate I.
    """
    m = cobra.Model("cell")
    S, inter, P, B = (cobra.Metabolite(x, compartment="c") for x in ("S", "I", "P", "B"))
    m.add_metabolites([S, inter, P, B])
    sup = cobra.Reaction("sup", lower_bound=0, upper_bound=10)  # -> S (substrate supply)
    sup.add_metabolites({S: 1})
    v1 = cobra.Reaction("v1", lower_bound=0, upper_bound=1000)
    v1.add_metabolites({S: -1, inter: 1})
    v2 = cobra.Reaction("v2", lower_bound=0, upper_bound=1000)
    v2.add_metabolites({inter: -1, P: 1})
    v3 = cobra.Reaction("v3", lower_bound=0, upper_bound=1000)
    v3.add_metabolites({inter: -1, B: 1})
    ex_p = cobra.Reaction("EX_P", lower_bound=0, upper_bound=1000)  # target product export
    ex_p.add_metabolites({P: -1})
    ex_b = cobra.Reaction("EX_B", lower_bound=0, upper_bound=1000)  # biomass
    ex_b.add_metabolites({B: -1})
    m.add_reactions([sup, v1, v2, v3, ex_p, ex_b])
    v1.gene_reaction_rule = "gA"
    v2.gene_reaction_rule = "gB"
    v3.gene_reaction_rule = "gC"
    m.objective = "EX_B"
    return m


def test_returns_result_with_scan(model):
    res = fseof(model, "EX_P", n_steps=8)
    assert isinstance(res, FSEOFResult)
    assert res.scan.shape[1] == len(res.enforced) >= 2
    assert "v2" in res.scan.index  # full scan retained, indexed by reaction


def test_amplification_targets(model):
    res = fseof(model, "EX_P", n_steps=8)
    amp = set(res.amplification["reaction"])
    # the product-forming reaction is amplified as EX_P is enforced upward
    # (v1/sup run at capacity throughout, so they stay constant and aren't flagged).
    assert {"v2", "EX_P"} <= amp
    v2 = res.targets.set_index("reaction").loc["v2"]
    assert v2["slope"] > 0 and v2["correlation"] > 0.9


def test_knockdown_of_competing_branch(model):
    res = fseof(model, "EX_P", n_steps=8)
    # v3 (biomass branch) competes for I -> suppressed toward zero -> knockdown/knockout
    down = set(res.knockout["reaction"])
    assert "v3" in down
    v3 = res.targets.set_index("reaction").loc["v3"]
    assert v3["slope"] < 0
    assert v3["target_type"] in ("knockdown", "knockout")


def test_gene_targets_aggregation(model):
    res = fseof(model, "EX_P", n_steps=8)
    genes = set(res.gene_targets["gene"])
    assert {"gA", "gB", "gC"} & genes  # reaction targets mapped to their genes
    gB = res.gene_targets.set_index("gene").loc["gB"]
    assert "v2" in gB["reactions"]


def test_unproducible_target_raises(model):
    # A reaction that cannot carry positive flux is not a valid product target.
    dead = cobra.Reaction("dead", lower_bound=0, upper_bound=0)
    dead.add_metabolites({model.metabolites.P: -1})
    model.add_reactions([dead])
    with pytest.raises(ValueError, match="cannot carry positive flux"):
        fseof(model, "dead")


def test_infeasible_model_raises_clear_error(model):
    """An infeasible model (slim_optimize -> NaN) raises the clear guard, not a NaN scan."""
    model.reactions.sup.bounds = (5, 5)  # force uptake while EX_P demands more -> infeasible
    model.reactions.EX_P.bounds = (1000, 1000)
    with pytest.raises(ValueError, match="cannot carry positive flux"):
        fseof(model, "EX_P", n_steps=4)
