"""Tests for simplifyModel reduction modes."""
import cobra
import pytest

from ravengem.manipulation import (
    add_reactions_from_equations,
    constrain_reversible_reactions,
    group_linear_reactions,
    remove_dead_end_reactions,
    remove_duplicate_reactions,
)

# --- remove_dead_end_reactions --------------------------------------------

def test_dead_end_removed():
    m = cobra.Model("t")
    m.add_metabolites([cobra.Metabolite(x, compartment="c") for x in ("a", "b", "dead")])
    add_reactions_from_equations(
        m,
        [
            {"id": "R_in", "equation": " --> a"},
            {"id": "R1", "equation": "a --> b"},
            {"id": "R_out", "equation": "b --> "},
            {"id": "R_dead", "equation": "a --> dead"},  # 'dead' only produced
        ],
    )
    removed_rxns, removed_mets = remove_dead_end_reactions(m)
    assert "R_dead" in removed_rxns
    assert "dead" in removed_mets
    # the productive path survives
    assert {"R_in", "R1", "R_out"} <= {r.id for r in m.reactions}


def test_dead_end_respects_reserved():
    m = cobra.Model("t")
    m.add_metabolites([cobra.Metabolite(x, compartment="c") for x in ("a", "dead")])
    add_reactions_from_equations(
        m, [{"id": "R_in", "equation": " --> a"}, {"id": "R_dead", "equation": "a --> dead"}]
    )
    removed_rxns, _ = remove_dead_end_reactions(m, reserved=["R_dead"])
    assert "R_dead" not in removed_rxns
    assert "R_dead" in {r.id for r in m.reactions}


# --- remove_duplicate_reactions -------------------------------------------

def test_duplicates_removed():
    m = cobra.Model("t")
    m.add_metabolites([cobra.Metabolite(x, compartment="c") for x in ("a", "b")])
    add_reactions_from_equations(
        m,
        [
            {"id": "R1", "equation": "a --> b", "bounds": (0, 1000)},
            {"id": "R2", "equation": "a --> b", "bounds": (0, 1000)},  # duplicate of R1
            {"id": "R3", "equation": "a --> b", "bounds": (0, 500)},   # different bounds
        ],
    )
    removed = remove_duplicate_reactions(m)
    assert len(removed) == 1  # one of R1/R2 removed
    assert {"R3"} <= {r.id for r in m.reactions}
    assert sum(r.id in ("R1", "R2") for r in m.reactions) == 1


def test_duplicates_keep_reserved():
    m = cobra.Model("t")
    m.add_metabolites([cobra.Metabolite(x, compartment="c") for x in ("a", "b")])
    add_reactions_from_equations(
        m,
        [
            {"id": "R1", "equation": "a --> b", "bounds": (0, 1000)},
            {"id": "R2", "equation": "a --> b", "bounds": (0, 1000)},
        ],
    )
    remove_duplicate_reactions(m, reserved=["R1"])
    assert "R1" in {r.id for r in m.reactions}  # reserved one kept


# --- constrain_reversible_reactions ---------------------------------------

def test_forward_only_reversible_constrained():
    m = cobra.Model("t")
    m.add_metabolites([cobra.Metabolite(x, compartment="c") for x in ("a", "b")])
    add_reactions_from_equations(
        m,
        [
            {"id": "R_in", "equation": " --> a", "bounds": (0, 1000)},
            {"id": "R1", "equation": "a <=> b", "bounds": (-1000, 1000)},  # can only go fwd
            {"id": "R_out", "equation": "b --> ", "bounds": (0, 1000)},
        ],
    )
    changed = constrain_reversible_reactions(m)
    assert "R1" in changed
    assert m.reactions.get_by_id("R1").lower_bound == 0  # constrained to forward


def test_truly_reversible_unchanged():
    m = cobra.Model("t")
    m.add_metabolites([cobra.Metabolite(x, compartment="c") for x in ("a", "b")])
    add_reactions_from_equations(
        m,
        [
            {"id": "R_in", "equation": " <=> a", "bounds": (-1000, 1000)},
            {"id": "R1", "equation": "a <=> b", "bounds": (-1000, 1000)},
            {"id": "R_out", "equation": "b <=> ", "bounds": (-1000, 1000)},
        ],
    )
    changed = constrain_reversible_reactions(m)
    assert "R1" not in changed  # can go both ways


# --- group_linear_reactions -----------------------------------------------

def test_linear_chain_merged():
    m = cobra.Model("t")
    m.add_metabolites([cobra.Metabolite(x, compartment="c") for x in ("a", "b", "c")])
    add_reactions_from_equations(
        m,
        [
            {"id": "R1", "equation": "a --> b"},  # b: single producer
            {"id": "R2", "equation": "b --> c"},  # b: single consumer
        ],
    )
    n_before = len(m.reactions)
    group_linear_reactions(m)
    # b is eliminated; R1+R2 merged into one reaction a --> c
    assert "b" not in m.metabolites
    assert len(m.reactions) < n_before


def test_group_linear_discards_genes():
    m = cobra.Model("t")
    m.add_metabolites([cobra.Metabolite(x, compartment="c") for x in ("a", "b", "c")])
    add_reactions_from_equations(
        m,
        [
            {"id": "R1", "equation": "a --> b", "gene_reaction_rule": "G1"},
            {"id": "R2", "equation": "b --> c", "gene_reaction_rule": "G2"},
        ],
    )
    group_linear_reactions(m)
    assert len(m.genes) == 0


# --- regression: incremental merge collapses a long chain (known_issues.md D1) ---

def test_group_linear_merges_long_chain_in_one_pass():
    """The incremental scan still flattens a 5-reaction linear chain — the
    correctness property the original O(n²·m) restart-after-merge loop had."""
    m = cobra.Model("t")
    m.add_metabolites([cobra.Metabolite(x, compartment="c") for x in "abcdef"])
    add_reactions_from_equations(
        m,
        [
            {"id": "R_in", "equation": " --> a"},
            {"id": "R1", "equation": "a --> b"},
            {"id": "R2", "equation": "b --> c"},
            {"id": "R3", "equation": "c --> d"},
            {"id": "R4", "equation": "d --> e"},
            {"id": "R5", "equation": "e --> f"},
            {"id": "R_out", "equation": "f --> "},
        ],
    )
    group_linear_reactions(m)
    # All the chain's internal metabolites are gone.
    assert {x for x in m.metabolites if x.id in {"b", "c", "d", "e"}} == set()


# --- regression: NaN FVA on infeasible model (known_issues.md C1) ----------

def test_constrain_reversible_raises_on_infeasible():
    """An infeasible model produces NaN FVA ranges; the old abs(NaN) < eps
    check silently treated those as 'truly reversible'. Now raises."""
    m = cobra.Model("t")
    a, b = (cobra.Metabolite(x, compartment="c") for x in ("a", "b"))
    m.add_metabolites([a, b])
    # Force a contradiction: r requires production AND consumption of a, but
    # nothing else produces a.
    r = cobra.Reaction("r", lower_bound=-1, upper_bound=1)
    r.add_metabolites({a: -1, b: 1})
    forced = cobra.Reaction("forced", lower_bound=5, upper_bound=10)  # infeasible
    forced.add_metabolites({a: -1})
    m.add_reactions([r, forced])
    with pytest.raises(RuntimeError, match="infeasible"):
        constrain_reversible_reactions(m)
