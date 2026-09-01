"""Tests for tINIT scoring + get_init_model (init/score.py, init/build.py)."""
import math

import cobra
import pytest

from raven_toolbox.init import (
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
def model(linear_chain_model_with_genes):
    # Shared linear-chain INIT model (with gene rules g1/g2/g3) — see tests/conftest.py.
    return linear_chain_model_with_genes


def test_get_init_model_from_gene_scores(model):
    # g1,g2 expressed (positive), g3 not (negative) -> keep r1,r2, drop r3.
    # allow_excretion=True: this model has no boundary reaction for B/C/D, so with
    # strict mass balance the only way r2 could carry flux would be for r3 to also
    # carry it (nothing else consumes the C that r2 produces) -- allow_excretion
    # lets the productive path run without dragging the dead-scored branch along.
    res = get_init_model(
        model, gene_scores={"g1": 5.0, "g2": 5.0, "g3": -5.0}, prod_weight=0.0,
        allow_excretion=True,
    )
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
    # allow_excretion=True: D (r3's product) has no boundary reaction, so under
    # strict mass balance forcing r3 on (essential) would be infeasible outright.
    res = get_init_model(
        model, rxn_scores={"r1": 1, "r2": 1, "r3": -1}, essential_rxns=["r3"], prod_weight=0.0,
        allow_excretion=True,
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
