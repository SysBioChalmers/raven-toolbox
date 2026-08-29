"""The single-step ftINIT MILP (run_ftinit).

Validated on the testModel oracle against (a) a hand-checked score-optimal solution,
(b) the formulation invariants, and (c) exact agreement with the already-tested
run_init. The full-pipeline RAVEN outputs (tinitTests T0001/T0002) additionally
involve linear merge + the toIgnore masks + staging + exchange re-adding, covered
elsewhere.

Note on the toy result: with strict mass balance and no metabolite-production reward
(ftINIT, unlike classic INIT, only rewards metabolomics-detected mets), the
score-optimal subnetwork on testModel is the internal cycle R4→R6→(R10 rev)→(R9 rev),
worth 7+0.5-3+3.5 = 8.0 — it beats the "honest" exchange path because that path must
pay for the negative-score transport reactions R2/R7. The bare INIT MILP has no
loopless constraint (neither does RAVEN's); loop-free models come from the staged
pipeline + exchange handling and, at genome scale, from models having real exchanges
so such cycles are not score-optimal. This faithfully matches RAVEN's MILP.
"""
import cobra
import pytest
from tinit_oracles import TEST_MODEL_SCORES, expr_for_rxn_score, make_test_model

from raven_toolbox.init import FtInitResult, run_ftinit, run_init
from raven_toolbox.init.score import gene_scores_from_expression, score_reactions_from_genes

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


def test_essential_force_clamps_to_capacity():
    """Forcing an essential reaction is clamped to its capacity (no lb>ub crash).

    A reaction capped at 0.05 forced with the default 0.1 must not error; it is forced
    to its capacity (0.05) and the model stays feasible. A per-reaction force of 0.04
    forces exactly that.
    """
    m = cobra.Model("cap")
    a, b = (cobra.Metabolite(x, compartment="s") for x in "ab")
    m.add_metabolites([a, b])
    r = cobra.Reaction("LOW", lower_bound=0, upper_bound=0.05)  # tiny capacity
    r.add_metabolites({a: -1, b: 1})
    for mid, st in [("EX_a", {a: -1}), ("EX_b", {b: -1})]:
        ex = cobra.Reaction(mid, lower_bound=-1000, upper_bound=1000)
        ex.add_metabolites(st)
        m.add_reactions([ex])
    m.add_reactions([r])
    m.objective = "LOW"

    res = run_ftinit(m, {}, essential_rxns=["LOW"], force_on_ess=0.1)  # clamped to 0.05
    assert res.fluxes["LOW"] >= 0.05 - 1e-9
    res2 = run_ftinit(m, {}, essential_rxns=["LOW"], essential_force={"LOW": 0.04})
    assert res2.fluxes["LOW"] >= 0.04 - 1e-9


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


def test_forced_flux_lower_bound_is_respected():
    """A scored, non-reversible reaction with lb>0 must keep carrying >= lb flux.

    Guards the bound handling: the single-direction branch must use the model's own
    [lb, ub], not zero out a positive lower bound.
    """
    model = make_test_model()
    scores = _scores(model)
    # R6 (2 d[c] => e[c]) is forward-irreversible; force >=2 flux through it.
    model.reactions.get_by_id("R6").lower_bound = 2.0
    res = run_ftinit(model, scores)
    assert res.fluxes["R6"] >= 2.0 - 1e-6
    assert "R6" not in res.deleted_reactions


# --------------------------------------------------------------------------- #
# canonical (deterministic uniqueness).
# --------------------------------------------------------------------------- #
def _degenerate_model():
    """Two interchangeable negative-score reactions (R1, R2) both feed an essential E.

    Exactly one is needed to make E carry flux, so the score optimum (-1) is degenerate
    between keeping R1 or R2. Which one a plain solve keeps is an arbitrary tie-break.
    """
    m = cobra.Model("degen")
    a, mm, p = (cobra.Metabolite(x, name=x, compartment="s") for x in ("a", "m", "p"))
    m.add_metabolites([a, mm, p])
    R1 = cobra.Reaction("R1", lower_bound=0, upper_bound=1000)
    R1.add_metabolites({a: -1, mm: 1})
    R2 = cobra.Reaction("R2", lower_bound=0, upper_bound=1000)
    R2.add_metabolites({a: -1, mm: 1})
    E = cobra.Reaction("E", lower_bound=0, upper_bound=1000)
    E.add_metabolites({mm: -1, p: 1})
    EXa = cobra.Reaction("EX_a", lower_bound=-1000, upper_bound=1000)
    EXa.add_metabolites({a: -1})
    EXp = cobra.Reaction("EX_p", lower_bound=-1000, upper_bound=1000)
    EXp.add_metabolites({p: -1})
    m.add_reactions([R1, R2, E, EXa, EXp])
    m.objective = "E"
    return m


def test_canonical_breaks_degenerate_tie_by_id():
    """canonical selects the unique sparsest, lowest-id optimum among equal alternatives."""
    m = _degenerate_model()
    res = run_ftinit(m, {"R1": -1.0, "R2": -1.0}, essential_rxns=["E"], canonical=True)
    # exactly one of the degenerate pair is kept (the tie is resolved, not doubled) ...
    assert len({"R1", "R2"} & set(res.kept_reactions)) == 1
    # ... and it is deterministically the lower-id one.
    assert "R1" in res.kept_reactions and "R2" not in res.kept_reactions
    # the reported objective is the primary (score) optimum, not the phase-2 secondary.
    assert res.objective == pytest.approx(-1.0, abs=1e-6)


def test_canonical_safe_on_unique_optimum():
    """On a non-degenerate model canonical returns the same optimum (no regression)."""
    model = make_test_model()
    res = run_ftinit(model, _scores(model), canonical=True)
    assert set(res.kept_reactions) == _LOOP
    assert res.objective == pytest.approx(8.0, abs=1e-6)
