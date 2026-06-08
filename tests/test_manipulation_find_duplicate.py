"""Tests for raven_python.manipulation.find_duplicate_reactions."""
from __future__ import annotations

import cobra

from raven_python.manipulation import find_duplicate_reactions


def _mk_model() -> cobra.Model:
    m = cobra.Model("m")
    a = cobra.Metabolite("a_c", compartment="c")
    b = cobra.Metabolite("b_c", compartment="c")
    c = cobra.Metabolite("c_c", compartment="c")
    m.add_metabolites([a, b, c])
    return m


def test_no_duplicates_returns_empty():
    m = _mk_model()
    r1 = cobra.Reaction("r1")
    r1.add_metabolites({m.metabolites.a_c: -1, m.metabolites.b_c: 1})
    r2 = cobra.Reaction("r2")
    r2.add_metabolites({m.metabolites.b_c: -1, m.metabolites.c_c: 1})
    m.add_reactions([r1, r2])
    assert find_duplicate_reactions(m) == []


def test_same_stoichiometry_grouped():
    m = _mk_model()
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r2 = cobra.Reaction("r2", lower_bound=-100, upper_bound=100)
    for r in (r1, r2):
        r.add_metabolites({m.metabolites.a_c: -1, m.metabolites.b_c: 1})
    m.add_reactions([r1, r2])
    groups = find_duplicate_reactions(m)
    assert len(groups) == 1
    assert {r.id for r in groups[0]} == {"r1", "r2"}


def test_ignore_direction_default_groups_reverse_pair():
    """yeast-GEM's findDuplicatedRxns matches A→B with B→A. That's the default."""
    m = _mk_model()
    r1 = cobra.Reaction("r1")
    r1.add_metabolites({m.metabolites.a_c: -1, m.metabolites.b_c: 1})
    r2 = cobra.Reaction("r2")
    r2.add_metabolites({m.metabolites.a_c: 1, m.metabolites.b_c: -1})  # reversed
    m.add_reactions([r1, r2])
    groups = find_duplicate_reactions(m)
    assert len(groups) == 1


def test_ignore_direction_false_keeps_them_separate():
    m = _mk_model()
    r1 = cobra.Reaction("r1")
    r1.add_metabolites({m.metabolites.a_c: -1, m.metabolites.b_c: 1})
    r2 = cobra.Reaction("r2")
    r2.add_metabolites({m.metabolites.a_c: 1, m.metabolites.b_c: -1})
    m.add_reactions([r1, r2])
    assert find_duplicate_reactions(m, ignore_direction=False) == []


def test_three_duplicates_in_one_group():
    m = _mk_model()
    rxns = []
    for i in range(3):
        r = cobra.Reaction(f"r{i}")
        r.add_metabolites({m.metabolites.a_c: -1, m.metabolites.b_c: 1})
        rxns.append(r)
    m.add_reactions(rxns)
    groups = find_duplicate_reactions(m)
    assert len(groups) == 1
    assert len(groups[0]) == 3


def test_ignores_bounds_and_gpr_differences():
    """Bounds and GPRs are intentionally ignored — only stoichiometry."""
    m = _mk_model()
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r2 = cobra.Reaction("r2", lower_bound=-50, upper_bound=50)
    r1.add_metabolites({m.metabolites.a_c: -1, m.metabolites.b_c: 1})
    r2.add_metabolites({m.metabolites.a_c: -1, m.metabolites.b_c: 1})
    r1.gene_reaction_rule = "g1"
    r2.gene_reaction_rule = "g2"
    m.add_reactions([r1, r2])
    assert len(find_duplicate_reactions(m)) == 1
