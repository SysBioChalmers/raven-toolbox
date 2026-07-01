"""Tests for functionality-constrained compartment assignment.

The toys are single-compartment drafts (compartment ``c``) that may be relocated into ``c``
or ``m``. Metabolites carry explicit names (cobra matches by id here, but names are set for
clarity). The decisive tests show that the *functionality* constraint — not just transport
cost — drives placement: a reaction goes to a compartment **against its own top localisation
score** when keeping its pathway able to carry flux requires it.
"""
import cobra
import pandas as pd
import pytest

from raven_toolbox.localization import (
    AssignmentProposal,
    apply_assignment,
    assign_compartments,
)
from raven_toolbox.localization.scores import LocalizationScores


def _met(mid):
    return cobra.Metabolite(mid, name=mid, compartment="c")


def _linear():
    """EX_A (uptake) -> r1: A->B -> bio: B-> (biomass). B produced only by r1."""
    m = cobra.Model("linear")
    A, B = _met("A_c"), _met("B_c")
    m.add_metabolites([A, B])
    ex = cobra.Reaction("EX_A", lower_bound=-10, upper_bound=1000)
    ex.add_metabolites({A: -1})
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({A: -1, B: 1})
    r1.gene_reaction_rule = "g1"
    bio = cobra.Reaction("bio", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({B: -1})
    m.add_reactions([ex, r1, bio])
    m.objective = "bio"
    return m


def _grows(model, tol=1e-6):
    v = model.slim_optimize(error_value=0.0)
    return v is not None and v > tol


# --------------------------------------------------------------------------- basics
def test_result_type():
    m = _linear()
    scores = LocalizationScores(pd.DataFrame({"c": [0.4], "m": [0.9]}, index=["g1"]))
    res = assign_compartments(m, scores, ["r1"], transport_cost=0.1)
    assert isinstance(res, AssignmentProposal)
    assert res.status == "optimal"


def test_score_respected_when_functional():
    # r1 prefers m (0.9 > 0.4); B is transportable and transport is cheap, so the model can
    # stay functional with r1 in m -> r1 should follow its score to m.
    m = _linear()
    scores = LocalizationScores(pd.DataFrame({"c": [0.4], "m": [0.9]}, index=["g1"]))
    res = assign_compartments(m, scores, ["r1"], transport_cost=0.1)
    assert res.placements["r1"] == ["m"]


def test_applied_model_is_functional():
    m = _linear()
    scores = LocalizationScores(pd.DataFrame({"c": [0.4], "m": [0.9]}, index=["g1"]))
    res = assign_compartments(m, scores, ["r1"], transport_cost=0.1)
    out = apply_assignment(m, res)
    assert _grows(out)  # biomass still producible after compartmentalisation


# ------------------------------------------------- functionality overrides the score
def test_functionality_overrides_score():
    # Same scores (r1 prefers m), but B is NOT transportable. Placing r1 in m would strand
    # B in m where biomass (in c) cannot reach it -> infeasible. So r1 must go to c, against
    # its higher m-score. This is the behaviour no score+transport optimiser without a
    # functionality constraint guarantees.
    m = _linear()
    scores = LocalizationScores(pd.DataFrame({"c": [0.4], "m": [0.9]}, index=["g1"]))
    res = assign_compartments(m, scores, ["r1"], transportable=["A"])  # B confined
    assert res.status == "optimal"
    assert res.placements["r1"] == ["c"]
    assert _grows(apply_assignment(m, res))


# --------------------------------------------------------------- pathway coherence
def _chain():
    """EX_A -> r1: A->X -> r2: X->B -> bio: B->. X is the pathway intermediate."""
    m = cobra.Model("chain")
    A, X, B = _met("A_c"), _met("X_c"), _met("B_c")
    m.add_metabolites([A, X, B])
    ex = cobra.Reaction("EX_A", lower_bound=-10, upper_bound=1000)
    ex.add_metabolites({A: -1})
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({A: -1, X: 1})
    r1.gene_reaction_rule = "g1"
    r2 = cobra.Reaction("r2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({X: -1, B: 1})
    r2.gene_reaction_rule = "g2"
    bio = cobra.Reaction("bio", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({B: -1})
    m.add_reactions([ex, r1, r2, bio])
    m.objective = "bio"
    return m


def test_pathway_coherence_overrides_individual_score():
    # g1 strongly prefers c; g2 prefers m. The intermediate X is non-transportable, so r1
    # and r2 *must share a compartment* (r2 can only run where r1 made X). With transport
    # free (no cost confound), the optimiser keeps the whole pathway in c — where the
    # higher-weight g1 anchors it — so r2 lands in c **against its own m-preference**.
    m = _chain()
    scores = LocalizationScores(pd.DataFrame(
        {"c": [0.9, 0.5], "m": [0.1, 0.9]}, index=["g1", "g2"]))
    res = assign_compartments(m, scores, ["r1", "r2"],
                              transportable=["A", "B"],  # X confined
                              transport_cost=0.0)
    assert res.status == "optimal"
    assert res.placements["r1"] == ["c"]
    assert res.placements["r2"] == ["c"]            # against g2's top score (m=0.9)
    assert _grows(apply_assignment(m, res))


# --------------------------------------------------------------------- gap-fill coupling
def _gap_draft():
    """EX_A -> r1: A->B; bio: B + C -> biomass. Nothing makes C, so the draft cannot grow."""
    m = cobra.Model("gap")
    A, B, C = _met("A_c"), _met("B_c"), _met("C_c")
    m.add_metabolites([A, B, C])
    ex = cobra.Reaction("EX_A", lower_bound=-10, upper_bound=1000)
    ex.add_metabolites({A: -1})
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({A: -1, B: 1})
    r1.gene_reaction_rule = "g1"
    bio = cobra.Reaction("bio", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({B: -1, C: -1})
    m.add_reactions([ex, r1, bio])
    m.objective = "bio"
    return m


def _universal():
    """Candidate database with rC: A -> C (the missing producer of C)."""
    u = cobra.Model("universal")
    A, C = _met("A_c"), _met("C_c")
    u.add_metabolites([A, C])
    rc = cobra.Reaction("rC", lower_bound=0, upper_bound=1000)
    rc.add_metabolites({A: -1, C: 1})
    u.add_reactions([rc])
    return u


def test_gapfill_restores_function():
    # The draft cannot grow (C has no producer). assign_compartments may pull rC from the
    # universal model to make biomass feasible. min_growth is passed explicitly because the
    # draft itself does not grow.
    m = _gap_draft()
    scores = LocalizationScores(pd.DataFrame({"c": [1.0], "m": [0.0]}, index=["g1"]))
    res = assign_compartments(m, scores, ["r1"], universal=_universal(), min_growth=1.0)
    assert res.status == "optimal"
    assert "rC" in res.added_reactions
    out = apply_assignment(m, res, universal=_universal())
    assert _grows(out)


def test_no_gratuitous_gapfill():
    # When the draft can already grow, a provided universal candidate costs but earns nothing,
    # so it is not added (parsimony).
    m = _linear()
    scores = LocalizationScores(pd.DataFrame({"c": [0.9], "m": [0.1]}, index=["g1"]))
    res = assign_compartments(m, scores, ["r1"], universal=_universal())
    assert res.status == "optimal"
    assert res.added_reactions == []


# --------------------------------------------------------------------- infeasible
def test_infeasible_growth_floor():
    m = _linear()
    scores = LocalizationScores(pd.DataFrame({"c": [0.4], "m": [0.9]}, index=["g1"]))
    res = assign_compartments(m, scores, ["r1"], min_growth=1e6)  # unreachable
    assert res.status != "optimal"
    assert res.placements == {}


def test_no_growth_draft_raises():
    m = _linear()
    m.reactions.EX_A.lower_bound = 0  # no uptake -> no growth
    scores = LocalizationScores(pd.DataFrame({"c": [0.4], "m": [0.9]}, index=["g1"]))
    with pytest.raises(ValueError, match="does not grow"):
        assign_compartments(m, scores, ["r1"])


# ------------------------------------------------- reaction-level multi-localization (opt-in)
def _dual_target():
    """EX_A -> r1: A->P -> r2: P->Q; bio: P + Q -> biomass. P is confined.

    g1 (on r1) prefers c; g2 (on r2) strongly prefers m. Because P cannot be transported,
    r2 running in m needs P made in m (r1 in m) while biomass needs P in c (r1 in c) — so a
    correct multi-localization places r1 in BOTH c and m.
    """
    m = cobra.Model("dual")
    A, P, Q = _met("A_c"), _met("P_c"), _met("Q_c")
    m.add_metabolites([A, P, Q])
    ex = cobra.Reaction("EX_A", lower_bound=-10, upper_bound=1000)
    ex.add_metabolites({A: -1})
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({A: -1, P: 1})
    r1.gene_reaction_rule = "g1"
    r2 = cobra.Reaction("r2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({P: -1, Q: 1})
    r2.gene_reaction_rule = "g2"
    bio = cobra.Reaction("bio", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({P: -1, Q: -1})
    m.add_reactions([ex, r1, r2, bio])
    m.objective = "bio"
    return m


_DUAL_SCORES = LocalizationScores(pd.DataFrame({"c": [0.9, 0.0], "m": [0.1, 1.0]}, index=["g1", "g2"]))


def test_mono_localization_is_the_default():
    # Default (multi_localization=False): each reaction gets exactly one compartment, and
    # functionality forces the confined pathway into c (r2 cannot reach m without P there).
    m = _dual_target()
    res = assign_compartments(m, _DUAL_SCORES, ["r1", "r2"], transportable=["A", "Q"],
                              transport_cost=0.1, multi_compartment_penalty=0.1)
    assert res.status == "optimal"
    assert len(res.placements["r1"]) == 1
    assert len(res.placements["r2"]) == 1
    assert _grows(apply_assignment(m, res))


def _multiplaced_carry_flux(model, res, eps=1e-4):
    """Every compartment a reaction is multi-placed in must carry >= eps flux (no dead duplicate).

    apply_assignment keeps the reaction id for its first compartment and names extras
    ``<rid>_<compartment>``; FVA each and require flux in any reaction placed in >1 compartment.
    """
    out = apply_assignment(model, res)
    for rid, comps in res.placements.items():
        if len(comps) < 2:
            continue
        rxn_ids = [rid] + [f"{rid}_{c}" for c in comps[1:]]
        for rxid in rxn_ids:
            with out:
                out.objective = rxid
                hi = abs(out.slim_optimize(error_value=0.0) or 0.0)
                out.reactions.get_by_id(rxid).objective_coefficient = -1
                lo = abs(out.slim_optimize(error_value=0.0) or 0.0)
            if max(hi, lo) < eps:
                return False  # a placement that can carry no flux -> dead duplicate
    return True


def test_multi_localization_is_function_driven():
    # Functionality requires r1 in BOTH c and m: P is confined, needed in c for biomass and in m for
    # r2 (which prefers m). The flux-activity coupling lets r1 occupy both — and crucially each
    # placement carries flux, so this is genuine dual-targeting, not score-harvesting. Sound on every
    # solver (it is enforced by constraints, not by tie-breaking).
    m = _dual_target()
    res = assign_compartments(m, _DUAL_SCORES, ["r1", "r2"], transportable=["A", "Q"],
                              transport_cost=0.1, multi_compartment_penalty=0.1,
                              multi_localization=True)
    assert res.status == "optimal"
    assert set(res.placements["r1"]) == {"c", "m"}     # function-driven multi-localization
    assert _grows(apply_assignment(m, res))            # still functional
    assert _multiplaced_carry_flux(m, res)             # and both r1 placements actually carry flux


def test_multi_localization_admits_no_dead_placement():
    # Regression for the dead-placement exploit. With the old capability pre-pass, raising the
    # transport cost made the solver "place" r2 in m carrying ZERO flux purely to harvest g2's
    # m-score (an outcome that survived even on Gurobi). The flux-activity coupling forbids it: any
    # extra placement must carry >= eps flux, so the solver returns genuine dual-targeting (with real
    # transports) or stays mono — never a dead duplicate.
    m = _dual_target()
    res = assign_compartments(m, _DUAL_SCORES, ["r1", "r2"], transportable=["A", "Q"],
                              transport_cost=0.3, multi_compartment_penalty=0.1,
                              multi_localization=True)
    assert res.status == "optimal"
    assert _multiplaced_carry_flux(m, res)             # no dead duplicate at any transport cost
    assert _grows(apply_assignment(m, res))


def test_multi_localization_handles_reversible_reactions():
    # Regression: a REVERSIBLE movable reaction (lb<0<ub) must still be able to carry flux under
    # multi_localization. The activity gate is an implication (aF=1 => v>=eps), not a hard bound, so
    # the reverse direction is not wrongly forbidden. A naive ``v >= eps*aF`` pins reversible flux to
    # zero and makes the model infeasible; this checks the model stays feasible and grows.
    m = _linear()
    m.reactions.r1.lower_bound = -1000.0  # make r1 reversible (it still runs forward for biomass)
    scores = LocalizationScores(pd.DataFrame({"c": [0.9], "m": [0.1]}, index=["g1"]))
    res = assign_compartments(m, scores, ["r1"], transportable=["A"],
                              multi_localization=True)
    assert res.status == "optimal"                     # was INFEASIBLE before the implication fix
    assert res.placements["r1"] == ["c"]               # B confined -> r1 stays in c, carrying flux
    assert _grows(apply_assignment(m, res))


def test_multi_localization_forbids_blocked_placement():
    # The pre-pass must not place a reaction where it is blocked. In the linear toy B is confined,
    # so r1 cannot carry flux in 'm' (no B sink there) — r1 stays mono in c even with multi on.
    m = _linear()
    scores = LocalizationScores(pd.DataFrame({"c": [0.4], "m": [0.9]}, index=["g1"]))
    res = assign_compartments(m, scores, ["r1"], transportable=["A"], multi_localization=True)
    assert res.status == "optimal"
    assert res.placements["r1"] == ["c"]               # 'm' is blocked -> not offered
    assert _grows(apply_assignment(m, res))
