"""Tests for tINIT scoring + get_init_model (init/score.py, init/build.py)."""
import math

import cobra
import pytest

from raven_python.init import (
    InitModelResult,
    gene_scores_from_expression,
    get_init_model,
    score_reactions_from_genes,
)


# --------------------------------------------------------------------------- #
# score_reactions_from_genes
# --------------------------------------------------------------------------- #
@pytest.fixture
def gpr_model():
    m = cobra.Model("g")
    a = cobra.Metabolite("a_c", compartment="c")
    b = cobra.Metabolite("b_c", compartment="c")
    m.add_metabolites([a, b])
    r_complex = cobra.Reaction("r_complex")  # (g1 and g2) or g3
    r_complex.add_metabolites({a: -1, b: 1})
    m.add_reactions([r_complex])
    r_complex.gene_reaction_rule = "(g1 and g2) or g3"
    r_nogene = cobra.Reaction("r_nogene")
    r_nogene.add_metabolites({b: -1})
    m.add_reactions([r_nogene])
    return m


def test_score_isozyme_max_complex_min(gpr_model):
    # (g1 and g2) or g3 -> max(min(1, 4), 3) = max(1, 3) = 3
    scores = score_reactions_from_genes(gpr_model, {"g1": 1.0, "g2": 4.0, "g3": 3.0})
    assert scores["r_complex"] == 3.0


def test_score_no_gene_reaction_gets_default(gpr_model):
    scores = score_reactions_from_genes(gpr_model, {"g1": 1, "g2": 1, "g3": 1}, no_gene_score=-2.0)
    assert scores["r_nogene"] == -2.0


def test_score_missing_genes_omitted(gpr_model):
    # g2 missing -> complex (g1 and g2) collapses to g1=1; OR with g3=3 -> max(1,3)=3
    scores = score_reactions_from_genes(gpr_model, {"g1": 1.0, "g3": 3.0})
    assert scores["r_complex"] == 3.0
    # all genes missing -> no_gene_score
    assert score_reactions_from_genes(gpr_model, {})["r_complex"] == -2.0


def test_score_invalid_method(gpr_model):
    with pytest.raises(ValueError, match="isozyme_scoring"):
        score_reactions_from_genes(gpr_model, {}, isozyme_scoring="nonsense")


# --------------------------------------------------------------------------- #
# gene_scores_from_expression (RNA-seq path)
# --------------------------------------------------------------------------- #
def test_expression_scores_sign_and_clamp():
    expr = {"hi": 100.0, "lo": 1.0, "mid": 10.0, "zero": 0.0}
    ref = 10.0  # threshold/reference
    s = gene_scores_from_expression(expr, ref)
    assert s["hi"] == pytest.approx(min(5 * math.log(10), 10.0))  # above ref -> positive
    assert s["lo"] == pytest.approx(max(5 * math.log(0.1), -5.0))  # below ref -> negative
    assert s["mid"] == pytest.approx(0.0)  # at ref -> 0
    assert s["zero"] == -5.0  # non-positive -> floor


def test_expression_per_gene_reference():
    expr = {"g": 20.0}
    s = gene_scores_from_expression(expr, {"g": 5.0})
    assert s["g"] == pytest.approx(5 * math.log(4))


# --------------------------------------------------------------------------- #
# get_init_model pipeline
# --------------------------------------------------------------------------- #
@pytest.fixture
def model():
    m = cobra.Model("net")
    A, B, C, D = (cobra.Metabolite(x, name=x[:-2], compartment="c") for x in ("A_c", "B_c", "C_c", "D_c"))
    m.add_metabolites([A, B, C, D])
    exa = cobra.Reaction("EX_A", lower_bound=-1000, upper_bound=1000)
    exa.add_metabolites({A: -1})
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({A: -1, B: 1})
    r2 = cobra.Reaction("r2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({B: -1, C: 1})
    r3 = cobra.Reaction("r3", lower_bound=0, upper_bound=1000)
    r3.add_metabolites({C: -1, D: 1})
    m.add_reactions([exa, r1, r2, r3])
    for r, rule in (("r1", "g1"), ("r2", "g2"), ("r3", "g3")):
        m.reactions.get_by_id(r).gene_reaction_rule = rule
    return m


def test_get_init_model_from_gene_scores(model):
    # g1,g2 expressed (positive), g3 not (negative) -> keep r1,r2, drop r3.
    res = get_init_model(model, gene_scores={"g1": 5.0, "g2": 5.0, "g3": -5.0}, prod_weight=0.0)
    assert isinstance(res, InitModelResult)
    kept = {r.id for r in res.model.reactions}
    assert {"r1", "r2"} <= kept
    assert "r3" not in kept
    assert res.reaction_scores["r1"] == 5.0


def test_get_init_model_requires_one_score_source(model):
    with pytest.raises(ValueError, match="exactly one"):
        get_init_model(model)
    with pytest.raises(ValueError, match="exactly one"):
        get_init_model(model, rxn_scores={}, gene_scores={})


def test_get_init_model_essential_kept(model):
    # r3 negative-scored but essential -> kept.
    res = get_init_model(
        model, rxn_scores={"r1": 1, "r2": 1, "r3": -1}, essential_rxns=["r3"], prod_weight=0.0
    )
    assert "r3" in {r.id for r in res.model.reactions}


def test_get_init_model_removes_dead_ends(model):
    # An isolated reaction that can never carry flux is dropped as a dead end.
    X, Y = cobra.Metabolite("X_c", compartment="c"), cobra.Metabolite("Y_c", compartment="c")
    dead = cobra.Reaction("dead", lower_bound=0, upper_bound=1000)
    dead.add_metabolites({X: -1, Y: 1})  # X has no source, Y no sink (no exchange)
    model.add_reactions([dead])
    res = get_init_model(model, rxn_scores={"r1": 1, "r2": 1}, prod_weight=0.0)
    assert "dead" in res.deleted_dead_end_reactions
    assert "dead" not in {r.id for r in res.model.reactions}
