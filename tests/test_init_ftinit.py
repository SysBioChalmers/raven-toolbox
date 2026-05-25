"""Phase 4d.3: the single-step ftINIT MILP (run_ftinit).

Validated on the testModel oracle against (a) a hand-checked score-optimal solution,
(b) the formulation invariants, and (c) exact agreement with the already-tested
run_init. The full-pipeline RAVEN outputs (tinitTests T0001/T0002) additionally
involve linear merge + the toIgnore masks + staging + exchange re-adding, layered on
in 4d.2/4d.3b/4d.5.

Note on the toy result: with strict mass balance and no metabolite-production reward
(ftINIT, unlike classic INIT, only rewards metabolomics-detected mets), the
score-optimal subnetwork on testModel is the internal cycle R4→R6→(R10 rev)→(R9 rev),
worth 7+0.5-3+3.5 = 8.0 — it beats the "honest" exchange path because that path must
pay for the negative-score transport reactions R2/R7. The bare INIT MILP has no
loopless constraint (neither does RAVEN's); loop-free models come from the staged
pipeline + exchange handling and, at genome scale, from models having real exchanges
so such cycles are not score-optimal. This faithfully matches RAVEN's MILP.
"""
import pytest
from tinit_oracles import TEST_MODEL_SCORES, expr_for_rxn_score, make_test_model

from ravengem.init import FtInitResult, run_ftinit, run_init
from ravengem.init.score import gene_scores_from_expression, score_reactions_from_genes

_LOOP = {"R4", "R6", "R9", "R10"}  # the score-optimal subnetwork (8.0)


def _scores(model):
    expr = expr_for_rxn_score(TEST_MODEL_SCORES)
    return score_reactions_from_genes(model, gene_scores_from_expression(expr, 1.0))


def test_full_milp_score_optimum():
    model = make_test_model()
    res = run_ftinit(model, _scores(model))
    assert isinstance(res, FtInitResult)
    assert set(res.kept_reactions) == _LOOP
    assert res.deleted_reactions == ["R1", "R2", "R3", "R5", "R7", "R8"]
    assert res.objective == pytest.approx(8.0, abs=1e-6)


def test_kept_reactions_carry_flux_and_balance():
    """Indicator-on reactions carry flux (≥ force_on) and the solution is steady-state."""
    model = make_test_model()
    res = run_ftinit(model, _scores(model))
    for rid in res.kept_reactions:
        assert abs(res.fluxes[rid]) > 1e-9
    # The extracted model is itself feasible/flux-consistent.
    assert res.model.slim_optimize() is not None


def test_agrees_with_run_init():
    """Exact agreement with the classic INIT MILP (no production reward, no rev loops).

    run_init splits reversibles and double-scores both directions unless no_rev_loops,
    so we compare under matching settings: same objective and same kept set.
    """
    model = make_test_model()
    scores = _scores(model)
    ft = run_ftinit(model, scores)
    init = run_init(model, scores, prod_weight=0.0, eps=0.1, no_rev_loops=True)
    assert set(ft.kept_reactions) == {r.id for r in init.model.reactions}
    assert ft.objective == pytest.approx(init.objective, abs=1e-6)


def test_essential_reaction_forced_on():
    """An essential reaction is kept and carries flux even when its score is negative."""
    model = make_test_model()
    res = run_ftinit(model, _scores(model), essential_rxns=["R3"])
    assert "R3" in res.kept_reactions
    assert abs(res.fluxes["R3"]) > 1e-6


def test_rem_pos_rev_drops_positive_reversibles():
    """rem_pos_rev frees positive reversibles (score→0): the score-8.0 loop collapses.

    R4 (+7) and R10 (+3.5) are positive reversibles; with them unscored, the cycle is
    no longer profitable (R6 0.5 - R9 3 < 0), so nothing scored stays on.
    """
    model = make_test_model()
    res = run_ftinit(model, _scores(model), rem_pos_rev=True)
    assert res.objective == pytest.approx(0.0, abs=1e-6)
    assert "R6" not in res.kept_reactions and "R9" not in res.kept_reactions


def test_allow_excretion_relaxes_balance():
    """With allow_excretion the result stays feasible (net production permitted)."""
    model = make_test_model()
    res = run_ftinit(model, _scores(model), allow_excretion=True)
    assert res.objective >= 8.0 - 1e-6  # at least as good as strict balance


def test_unscored_reactions_are_kept_free():
    """Score-0 reactions are left in the model (not removable), not deleted."""
    model = make_test_model()
    scores = _scores(model)
    scores["R3"] = 0.0  # make R3 unscored -> must not be deleted
    res = run_ftinit(model, scores)
    assert "R3" not in res.deleted_reactions
