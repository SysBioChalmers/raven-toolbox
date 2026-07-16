"""Tests for :func:`curation_priority` — the post-hoc "which placements/transports to review" score.

Two layers: (1) end-to-end toys run through assign_compartments so the signals fire on real
proposals, and (2) unit-level toys with hand-built AssignmentProposals that pin the exact placement /
transport structure a signal needs, so each signal (and its v1.1 hardening — dual-localisation,
sibling-compartment robustness, currency/cofactor keying, essentiality-vs-baseline) is asserted in
isolation. The combination maths (noisy-OR, the essentiality stakes multiplier) and the ranking
contract are checked on top.
"""
import cobra
import pandas as pd

from raven_toolbox.localization import (
    AssignmentProposal,
    assign_compartments,
    curation_priority,
)
from raven_toolbox.localization.curation import (
    _is_currency,
    _is_impermeant,
    _noisy_or,
    _norm_key,
)
from raven_toolbox.localization.scores import LocalizationScores


def _met(mid, name=None):
    return cobra.Metabolite(mid, name=name or mid, compartment="c")


def _proposal(placements, *, transports=(), added=(), unplaced=(), min_growth=1.0):
    return AssignmentProposal(
        placements={k: list(v) for k, v in placements.items()},
        added_transports=list(transports),
        added_reactions=list(added),
        unplaced_reactions=list(unplaced),
        min_growth=min_growth,
        status="optimal",
    )


def _linear():
    """EX_A -> r1: A->B -> bio: B->."""
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


def _gap_draft():
    """EX_A -> r1: A->B; bio: B + C ->. Nothing makes C: the draft cannot grow without gap-fill."""
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
    u = cobra.Model("universal")
    A, C = _met("A_c"), _met("C_c")
    u.add_metabolites([A, C])
    rc = cobra.Reaction("rC", lower_bound=0, upper_bound=1000)
    rc.add_metabolites({A: -1, C: 1})
    u.add_reactions([rc])
    return u


def _coa_chain():
    """EX_A -> r0(c): A->palCoA -> r1: palCoA->B -> bio: B->. r1 prefers m, so palmitoyl-CoA (an
    impermeant acyl-CoA) has to be transported into m to feed it."""
    m = cobra.Model("coa")
    A = _met("A_c")
    pal = _met("palCoA_c", name="palmitoyl-CoA")
    B = _met("B_c")
    m.add_metabolites([A, pal, B])
    ex = cobra.Reaction("EX_A", lower_bound=-10, upper_bound=1000)
    ex.add_metabolites({A: -1})
    r0 = cobra.Reaction("r0", lower_bound=0, upper_bound=1000)
    r0.add_metabolites({A: -1, pal: 1})
    r0.gene_reaction_rule = "g0"
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({pal: -1, B: 1})
    r1.gene_reaction_rule = "g1"
    bio = cobra.Reaction("bio", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({B: -1})
    m.add_reactions([ex, r0, r1, bio])
    m.objective = "bio"
    return m


# ---------------------------------------------------------------------------
# end-to-end (real proposals)
# ---------------------------------------------------------------------------

def test_output_contract_and_ranking():
    m = _linear()
    scores = LocalizationScores(pd.DataFrame({"c": [0.4], "m": [0.9]}, index=["g1"]))
    prop = assign_compartments(m, scores, ["r1"], default_compartment="c")
    cp = curation_priority(m, prop, scores)
    assert list(cp.columns) == ["target", "type", "compartment", "priority", "flags", "action",
                                "n_affected", "affected"]
    assert cp["priority"].is_monotonic_decreasing
    assert (cp["priority"] >= 0).all() and (cp["priority"] <= 1).all()
    assert cp["action"].map(lambda s: s.isascii()).all()


def test_sandwich_flags_the_relocated_reaction():
    m = _linear()
    scores = LocalizationScores(pd.DataFrame({"c": [0.4], "m": [0.9]}, index=["g1"]))
    prop = assign_compartments(m, scores, ["r1"], default_compartment="c")
    assert prop.placements["r1"] == ["m"]
    row = curation_priority(m, prop, scores).query("target == 'r1'")
    assert not row.empty
    assert "sandwich" in row.iloc[0]["flags"]


def test_override_fires_when_function_beats_evidence():
    m = _linear()
    scores = LocalizationScores(pd.DataFrame({"c": [0.4], "m": [0.9]}, index=["g1"]))
    prop = assign_compartments(m, scores, ["r1"], default_compartment="c", transportable=["A"])
    assert prop.placements["r1"] == ["c"]
    row = curation_priority(m, prop, scores).query("target == 'r1'")
    assert not row.empty
    assert "override" in row.iloc[0]["flags"]


def test_confident_aligned_placement_is_not_flagged():
    m = _linear()
    scores = LocalizationScores(pd.DataFrame({"c": [0.95], "m": [0.05]}, index=["g1"]))
    prop = assign_compartments(m, scores, ["r1"], default_compartment="c")
    assert prop.placements["r1"] == ["c"]
    cp = curation_priority(m, prop, scores)
    assert "r1" not in set(cp["target"])


def test_gapfill_reaction_is_flagged():
    m = _gap_draft()
    scores = LocalizationScores(pd.DataFrame({"c": [0.9], "m": [0.1]}, index=["g1"]))
    prop = assign_compartments(
        m, scores, ["r1"], default_compartment="c", universal=_universal(), min_growth=0.5)
    assert "rC" in prop.added_reactions
    cp = curation_priority(m, prop, scores, universal=_universal())
    rc = cp.query("target == 'rC'")
    assert not rc.empty
    assert rc.iloc[0]["type"] == "gapfill"
    assert "gapfill" in rc.iloc[0]["flags"]


def test_impermeant_transport_is_flagged():
    m = _coa_chain()
    scores = LocalizationScores(pd.DataFrame(
        {"c": [0.5, 0.1], "m": [0.4, 0.9]}, index=["g0", "g1"]))
    prop = assign_compartments(m, scores, ["r0", "r1"], default_compartment="c")
    assert prop.placements["r1"] == ["m"]
    cp = curation_priority(m, prop, scores)
    imperm = cp[cp["flags"].str.contains("impermeant")]
    assert not imperm.empty
    assert (imperm["type"] == "transport").all()


def test_essentiality_amplifies_priority():
    # Use the override scenario (r1 forced to c against m-evidence) so the base priority clears the
    # essential-candidate threshold: a strongly-supported placement is deliberately de-prioritised by
    # the evidence gate and would (correctly) not be an essentiality candidate.
    m = _linear()
    scores = LocalizationScores(pd.DataFrame({"c": [0.4], "m": [0.9]}, index=["g1"]))
    prop = assign_compartments(m, scores, ["r1"], default_compartment="c", transportable=["A"])
    assert prop.placements["r1"] == ["c"]
    with_ess = curation_priority(m, prop, scores, check_essential=True).query("target == 'r1'")
    without = curation_priority(m, prop, scores, check_essential=False).query("target == 'r1'")
    assert "essential" in with_ess.iloc[0]["flags"]
    assert with_ess.iloc[0]["priority"] > without.iloc[0]["priority"]
    assert "essential" not in without.iloc[0]["flags"]


# ---------------------------------------------------------------------------
# override hardening (dual-localisation + sibling-compartment robustness)
# ---------------------------------------------------------------------------

def test_override_not_fired_when_top_is_a_secondary_placement():
    # r1 is dual-localised in c AND m; its gene's top compartment m IS a placement -> nothing overridden.
    m = _linear()
    scores = LocalizationScores(pd.DataFrame({"c": [0.3], "m": [1.0]}, index=["g1"]))
    prop = _proposal({"r1": ["c", "m"]})
    cp = curation_priority(m, prop, scores, check_essential=False)
    assert not any("override" in f for f in cp[cp["target"] == "r1"]["flags"])


def test_override_robust_to_sibling_compartment():
    # A near-tied sibling (mito matrix vs membrane, m=1.0/mm=0.98) must NOT collapse the override score:
    # the graded value is score(top) - score(assigned) = 1.0 - 0.2, not top1 - top2 (= 0.02).
    m = _linear()
    scores = LocalizationScores(pd.DataFrame({"c": [0.2], "m": [1.0], "mm": [0.98]}, index=["g1"]))
    prop = _proposal({"r1": ["c"]})
    row = curation_priority(m, prop, scores, check_essential=False).query("target == 'r1'").iloc[0]
    assert row["flags"].startswith("override")
    assert row["priority"] > 0.6  # would be ~0.018 under the collapsing top1-top2 margin


def test_dual_localised_emits_a_row_per_compartment():
    m = _linear()
    scores = LocalizationScores(pd.DataFrame({"c": [0.5], "m": [0.5], "x": [1.0]}, index=["g1"]))
    prop = _proposal({"r1": ["c", "m"]}, unplaced=["r1"])  # no_evidence forces a flag on both rows
    cp = curation_priority(m, prop, scores, check_essential=False)
    assert set(cp[cp["target"] == "r1"]["compartment"]) == {"c", "m"}


# ---------------------------------------------------------------------------
# sandwich grading
# ---------------------------------------------------------------------------

def test_sandwich_graded_by_pass_through_fraction():
    m = cobra.Model("sw")
    A, Bx, C = _met("A_c"), _met("Bx_c"), _met("C_c")
    m.add_metabolites([A, Bx, C])
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({A: -1, Bx: -1, C: 1})  # substrates A, Bx ; product C
    r1.gene_reaction_rule = "g1"
    m.add_reactions([r1])
    scores = LocalizationScores(pd.DataFrame({"c": [0.5], "m": [0.5]}, index=["g1"]))
    partial = curation_priority(  # only A of {A,Bx} imported -> fs=0.5
        m, _proposal({"r1": ["m"]}, transports=[("A", "m"), ("C", "m")]),
        scores, check_essential=False).query("target == 'r1'").iloc[0]
    full = curation_priority(  # both substrates imported -> fs=1.0 (pure pass-through)
        m, _proposal({"r1": ["m"]}, transports=[("A", "m"), ("Bx", "m"), ("C", "m")]),
        scores, check_essential=False).query("target == 'r1'").iloc[0]
    assert "sandwich" in partial["flags"] and "sandwich" in full["flags"]
    assert full["priority"] > partial["priority"]


def test_topology_signals_gated_by_evidence_support():
    # Identical sandwich topology, opposite evidence: a placement DeepLoc strongly supports is
    # de-prioritised, an unsupported one is surfaced. This gate is what keeps the ranking selective.
    m = _linear()
    prop = _proposal({"r1": ["m"]}, transports=[("A", "m"), ("B", "m")])
    supported = curation_priority(
        m, prop, LocalizationScores(pd.DataFrame({"c": [0.05], "m": [0.95]}, index=["g1"])),
        check_essential=False).query("target == 'r1'").iloc[0]["priority"]
    unsupported = curation_priority(
        m, prop, LocalizationScores(pd.DataFrame({"c": [0.05], "m": [0.05]}, index=["g1"])),
        check_essential=False).query("target == 'r1'").iloc[0]["priority"]
    assert unsupported > 3 * supported


# ---------------------------------------------------------------------------
# essentiality vs the applied baseline (not the heuristic floor)
# ---------------------------------------------------------------------------

def test_essentiality_skipped_when_materialised_growth_below_floor():
    # A proposal whose floor (5.0) exceeds what the materialised model can reach (bio maxes at 10 with
    # uptake 10, but we set an unreachable floor higher than achievable after the essential-only path):
    # here we force the pathological case directly by claiming a floor above the applied optimum.
    m = _linear()
    scores = LocalizationScores(pd.DataFrame({"c": [0.2], "m": [1.0]}, index=["g1"]))
    # r1 placed in c against top m -> override flags it as a candidate; floor absurdly high.
    prop = _proposal({"r1": ["c"]}, min_growth=1e6)
    cp = curation_priority(m, prop, scores, check_essential=True)
    row = cp.query("target == 'r1'").iloc[0]
    # essentiality cannot discriminate (applied growth < floor), so no essential flag / no amplification.
    assert "essential" not in row["flags"]


# ---------------------------------------------------------------------------
# keying: currency normalisation and impermeant cargo
# ---------------------------------------------------------------------------

def test_norm_key_strips_prefix_and_bracket():
    assert _norm_key("M_atp") == "atp"
    assert _norm_key("atp[c]") == "atp"
    assert _norm_key("  H2O ") == "h2o"


def test_currency_recognised_across_id_conventions():
    assert _is_currency("ATP", set())
    assert _is_currency("M_atp", set())    # SBML M_ prefix (compartment already stripped by base fn)
    assert _is_currency("atp[c]", set())   # bracketed-compartment id
    assert _is_currency("h2o", set())
    assert not _is_currency("palmitoyl-CoA", set())
    assert _is_currency("myCustom", {"mycustom"})


def test_impermeant_recognises_cofactors_and_thioesters():
    assert _is_impermeant("nadh", "NADH")         # free reduced cofactor -> needs a shuttle
    assert _is_impermeant("nadph", "NADPH")
    assert _is_impermeant("fad", "FAD")
    assert _is_impermeant("palCoA", "palmitoyl-CoA")
    assert _is_impermeant("s_0001", "some acyl-ACP")
    assert not _is_impermeant("glc", "D-glucose")


# ---------------------------------------------------------------------------
# combination maths
# ---------------------------------------------------------------------------

def test_noisy_or_rewards_corroboration():
    one = _noisy_or({"sandwich": 1.0})
    two = _noisy_or({"sandwich": 1.0, "override": 1.0})
    assert two > one
    assert two < 1.0


# ---------------------------------------------------------------------------
# coupling (which reactions a curation would affect)
# ---------------------------------------------------------------------------

def test_affected_lists_coupled_reactions_at_same_compartment():
    # R, S1, S2 all placed in m and all touch the shared hub metabolite HUB -> curating R affects both.
    m = cobra.Model("hub")
    A, hub, B, D = _met("A_c"), _met("HUB_c", name="HUB"), _met("B_c"), _met("D_c")
    m.add_metabolites([A, hub, B, D])
    R = cobra.Reaction("R", lower_bound=0, upper_bound=1000)
    R.add_metabolites({A: -1, hub: 1})
    R.gene_reaction_rule = "g1"
    S1 = cobra.Reaction("S1", lower_bound=0, upper_bound=1000)
    S1.add_metabolites({hub: -1, B: 1})
    S2 = cobra.Reaction("S2", lower_bound=0, upper_bound=1000)
    S2.add_metabolites({hub: -1, D: 1})
    m.add_reactions([R, S1, S2])
    scores = LocalizationScores(pd.DataFrame({"c": [0.1], "m": [0.1]}, index=["g1"]))
    prop = _proposal({"R": ["m"], "S1": ["m"], "S2": ["m"]}, unplaced=["R"])
    row = curation_priority(m, prop, scores, check_essential=False).query("target == 'R'").iloc[0]
    assert row["n_affected"] == 2
    assert set(row["affected"].replace(" ", "").split(",")) == {"S1", "S2"}


def test_affected_ignores_reactions_in_other_compartments():
    # S2 is placed in a DIFFERENT compartment (p): it shares HUB but is not coupled to R at m.
    m = cobra.Model("hub2")
    A, hub, B, D = _met("A_c"), _met("HUB_c", name="HUB"), _met("B_c"), _met("D_c")
    m.add_metabolites([A, hub, B, D])
    R = cobra.Reaction("R", lower_bound=0, upper_bound=1000)
    R.add_metabolites({A: -1, hub: 1})
    R.gene_reaction_rule = "g1"
    S1 = cobra.Reaction("S1", lower_bound=0, upper_bound=1000)
    S1.add_metabolites({hub: -1, B: 1})
    S2 = cobra.Reaction("S2", lower_bound=0, upper_bound=1000)
    S2.add_metabolites({hub: -1, D: 1})
    m.add_reactions([R, S1, S2])
    scores = LocalizationScores(pd.DataFrame({"c": [0.1], "m": [0.1]}, index=["g1"]))
    prop = _proposal({"R": ["m"], "S1": ["m"], "S2": ["p"]}, unplaced=["R"])
    row = curation_priority(m, prop, scores, check_essential=False).query("target == 'R'").iloc[0]
    assert set(row["affected"].replace(" ", "").split(",")) == {"S1"}


def test_affected_for_transport_lists_dependent_reactions():
    m = _coa_chain()
    scores = LocalizationScores(pd.DataFrame({"c": [0.5, 0.1], "m": [0.4, 0.9]}, index=["g0", "g1"]))
    prop = assign_compartments(m, scores, ["r0", "r1"], default_compartment="c")
    t = curation_priority(m, prop, scores)
    trow = t[t["type"] == "transport"].iloc[0]
    # the palmitoyl-CoA transport is coupled to the reactions that make/consume it (r0 in c, r1 in m)
    assert trow["n_affected"] >= 1
    assert "r1" in trow["affected"]
