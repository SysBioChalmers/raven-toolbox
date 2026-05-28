"""Tests for manipulation/compartments.py — merge_compartments + copy_to_compartment."""
from __future__ import annotations

import cobra
import pytest

from ravengem.manipulation.compartments import copy_to_compartment, merge_compartments


def _two_compartment_model() -> cobra.Model:
    """A_c → B_c, A_m → B_m, and a transport A_c ↔ A_m. Multi-compartment toy."""
    m = cobra.Model("toy")
    A_c = cobra.Metabolite("A_c", name="A", compartment="c")
    A_m = cobra.Metabolite("A_m", name="A", compartment="m")
    B_c = cobra.Metabolite("B_c", name="B", compartment="c")
    B_m = cobra.Metabolite("B_m", name="B", compartment="m")
    m.add_metabolites([A_c, A_m, B_c, B_m])

    def rxn(rid, lb, ub, mets, gpr=None):
        r = cobra.Reaction(rid, lower_bound=lb, upper_bound=ub)
        r.add_metabolites(mets)
        if gpr:
            r.gene_reaction_rule = gpr
        return r
    m.add_reactions([rxn("r_c", 0, 1000, {A_c: -1, B_c: 1}, "g1"),
                     rxn("r_m", 0, 1000, {A_m: -1, B_m: 1}, "g2"),
                     rxn("tr_A", -1000, 1000, {A_c: -1, A_m: 1})])
    return m


# ----------------------------------------------------------------- merge_compartments

def test_merge_compartments_collapses_to_one():
    """A_c + A_m → A; B_c + B_m → B; transport A_c↔A_m self-cancels and is dropped."""
    m = _two_compartment_model()
    merged, deleted, dupes = merge_compartments(m)
    # Only the base ids survive.
    assert {x.id for x in merged.metabolites} == {"A", "B"}
    # The transport reaction collapsed (A → A) and was deleted.
    assert "tr_A" in deleted
    # r_c and r_m are now both A → B; one of them gets deduplicated.
    surviving = {r.id for r in merged.reactions}
    assert len(surviving & {"r_c", "r_m"}) == 1
    assert (set(dupes) | (surviving & {"r_c", "r_m"})) == {"r_c", "r_m"}


def test_merge_compartments_preserves_gpr_and_subsystem():
    m = _two_compartment_model()
    m.reactions.r_c.subsystem = "carbo"
    merged, _, _ = merge_compartments(m)
    survivor = next(r for r in merged.reactions if r.id in {"r_c", "r_m"})
    # The survivor keeps its gene rule + subsystem (cobra may sometimes lose them
    # through copy; we set them explicitly).
    assert survivor.gene_reaction_rule in {"g1", "g2"}
    if survivor.id == "r_c":
        assert survivor.subsystem == "carbo"


def test_merge_compartments_keeps_single_met_reactions_when_asked():
    """drop_single_metabolite_reactions=False keeps the collapsed transport (now A → A,
    which is empty stoichiometry after net-cancellation — still dropped, but the *one-met*
    case is the more interesting one). Use a uniport pattern to exercise it."""
    m = cobra.Model("uniport")
    A_c = cobra.Metabolite("A_c", name="A", compartment="c")
    A_m = cobra.Metabolite("A_m", name="A", compartment="m")
    H_c = cobra.Metabolite("H_c", name="H", compartment="c")
    m.add_metabolites([A_c, A_m, H_c])
    # H+ symport: A_c + H_c → A_m. After merge: A + H → A → leaves H.
    sym = cobra.Reaction("sym", lower_bound=0, upper_bound=1000)
    sym.add_metabolites({A_c: -1, H_c: -1, A_m: 1})
    m.add_reactions([sym])
    merged_drop, deleted_drop, _ = merge_compartments(m, drop_single_metabolite_reactions=True)
    assert "sym" in deleted_drop
    merged_keep, deleted_keep, _ = merge_compartments(m, drop_single_metabolite_reactions=False)
    # With keep, sym survives as a one-met reaction (consumes H).
    assert "sym" not in deleted_keep
    assert "sym" in {r.id for r in merged_keep.reactions}


def test_merge_compartments_deduplicate_off_keeps_both():
    m = _two_compartment_model()
    merged, _, dupes = merge_compartments(m, deduplicate_reactions=False)
    assert dupes == []
    assert {"r_c", "r_m"} <= {r.id for r in merged.reactions}


# ----------------------------------------------------------------- copy_to_compartment

def test_copy_to_compartment_basic():
    """Copy r_c into 'p' (peroxisome): a new reaction r_c_p with metabolites in p."""
    m = _two_compartment_model()
    out, new_rxns, new_mets = copy_to_compartment(m, ["r_c"], "p",
                                                    target_compartment_name="peroxisome")
    assert "r_c_p" in [r.id for r in out.reactions]
    new_r = out.reactions.r_c_p
    assert {x.compartment for x in new_r.metabolites} == {"p"}
    assert "A_p" in [x.id for x in out.metabolites]
    assert "B_p" in [x.id for x in out.metabolites]
    assert new_rxns == ["r_c_p"]
    assert set(new_mets) == {"A_p", "B_p"}
    # Original still there.
    assert "r_c" in [r.id for r in out.reactions]


def test_copy_to_compartment_preserves_gpr_and_bounds():
    m = _two_compartment_model()
    out, _, _ = copy_to_compartment(m, ["r_c"], "p")
    new_r = out.reactions.r_c_p
    assert new_r.gene_reaction_rule == "g1"
    assert new_r.lower_bound == 0 and new_r.upper_bound == 1000


def test_copy_to_compartment_delete_original_is_a_move():
    m = _two_compartment_model()
    out, _, _ = copy_to_compartment(m, ["r_c"], "p", delete_original=True)
    assert "r_c" not in [r.id for r in out.reactions]
    assert "r_c_p" in [r.id for r in out.reactions]


def test_copy_to_compartment_idempotent():
    """Calling twice doesn't add the reaction twice."""
    m = _two_compartment_model()
    out, _, _ = copy_to_compartment(m, ["r_c"], "p")
    out2, new_rxns, _ = copy_to_compartment(out, ["r_c"], "p")
    assert new_rxns == []  # nothing added on second call
    assert len([r for r in out2.reactions if r.id == "r_c_p"]) == 1


def test_copy_to_compartment_unknown_reaction_raises():
    m = _two_compartment_model()
    with pytest.raises(ValueError, match="not in model"):
        copy_to_compartment(m, ["does_not_exist"], "p")


def test_copy_to_compartment_custom_suffix():
    m = _two_compartment_model()
    out, new_rxns, _ = copy_to_compartment(m, ["r_c"], "p", id_suffix="copy1")
    assert new_rxns == ["r_c_copy1"]
    assert "A_copy1" in [x.id for x in out.metabolites]
