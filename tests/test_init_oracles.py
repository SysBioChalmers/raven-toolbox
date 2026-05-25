"""Validate the ftINIT toy oracles and that our scoring reproduces RAVEN's.

This is Phase 4d.0: the correctness scaffold. The (ft)INIT MILP itself is not yet
ported, so the on/off-output oracles in tinit_oracles live there as constants for the
later sub-phases; here we lock down the pieces that already exist — the score→
expression inversion and scoreComplexModel-equivalent scoring (RAVEN tinitTests
T0009).
"""
import pytest
from tinit_oracles import (
    TEST_MODEL4_SCORES,
    TEST_MODEL5_SCORES,
    TEST_MODEL_SCORES,
    expr_for_rxn_score,
    make_test_model,
    make_test_model4,
    make_test_model5,
)

from ravengem.init.score import gene_scores_from_expression, score_reactions_from_genes


@pytest.mark.parametrize(
    "make_model, scores",
    [
        (make_test_model, TEST_MODEL_SCORES),
        (make_test_model4, TEST_MODEL4_SCORES),
        (make_test_model5, TEST_MODEL5_SCORES),
    ],
)
def test_scoring_reproduces_defined_scores(make_model, scores):
    """RAVEN T0009: expr_for_rxn_score → scoreComplexModel round-trips the scores."""
    model = make_model()
    expression = expr_for_rxn_score(scores)
    gene_scores = gene_scores_from_expression(expression, 1.0)
    rxn_scores = score_reactions_from_genes(model, gene_scores)
    got = [rxn_scores[r.id] for r in model.reactions]
    assert got == pytest.approx(scores, abs=1e-10)


def test_expr_for_rxn_score_inverts_scoring():
    """level = exp(score/5); 5·ln(level/1) recovers the score."""
    scores = [-5, -1, 0.5, 7, 10]
    expr = expr_for_rxn_score(scores)
    recovered = gene_scores_from_expression(expr, 1.0)
    assert [recovered[f"G{i + 1}"] for i in range(len(scores))] == pytest.approx(scores)


def test_test_model_structure():
    """Sanity: shapes, no-GPR reactions, reversibility, objective."""
    m = make_test_model()
    assert len(m.reactions) == 10 and len(m.metabolites) == 8
    no_gpr = {r.id for r in m.reactions if not r.genes}
    assert no_gpr == {"R1", "R2", "R8"}  # the reactions scored -2 (no gene)
    rev = {r.id for r in m.reactions if r.lower_bound < 0}
    assert rev == {"R2", "R3", "R4", "R9", "R10"}
    assert m.objective.expression.as_coefficients_dict()  # objective set (R8)


def test_test_model_is_feasible_for_the_task():
    """The toy model can actually make e[s] from a[s] (so the task oracle is meaningful)."""
    m = make_test_model()
    m.objective = "R8"
    assert m.slim_optimize() > 1e-6
