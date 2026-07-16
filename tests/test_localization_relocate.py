"""Tests for :func:`relocate_reactions` — apply a curator's firm placement and repair the model.

The toys pin a reaction (or gene-sibling group) into a new compartment and assert the *consequential*
changes the function must make: transports bridged into the new compartment, orphaned transports at the
old one removed, gene-linked reactions co-moved (but isozyme-shared ones left), impermeant cargo warned
about, and the materialised model re-certified for growth.
"""
import cobra

from raven_toolbox.localization import (
    AssignmentProposal,
    relocate_reactions,
)


def _met(mid, name=None):
    return cobra.Metabolite(mid, name=name or mid, compartment="c")


def _proposal(placements, *, transports=(), min_growth=0.1):
    return AssignmentProposal(
        placements={k: list(v) for k, v in placements.items()},
        added_transports=list(transports),
        min_growth=min_growth,
        status="optimal",
    )


def _chain3(shared_gene=False, isozyme=False):
    """EX_A -> r1: A->B -> r2: B->C -> bio: C->. Genes: r1=g1; r2=g1 (shared) / g1+g2 (isozyme) / g2."""
    m = cobra.Model("chain3")
    A, B, C = _met("A_c"), _met("B_c"), _met("C_c")
    m.add_metabolites([A, B, C])
    ex = cobra.Reaction("EX_A", lower_bound=-10, upper_bound=1000)
    ex.add_metabolites({A: -1})
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({A: -1, B: 1})
    r1.gene_reaction_rule = "g1"
    r2 = cobra.Reaction("r2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({B: -1, C: 1})
    r2.gene_reaction_rule = "g1" if shared_gene else ("g1 or g2" if isozyme else "g2")
    bio = cobra.Reaction("bio", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({C: -1})
    m.add_reactions([ex, r1, r2, bio])
    m.objective = "bio"
    return m


def _grows(model):
    v = model.slim_optimize(error_value=0.0)
    return v is not None and v > 1e-9


def test_relocation_adds_bridging_transports_and_stays_functional():
    m = _chain3()
    prop = _proposal({"r1": ["c"], "r2": ["c"]})
    res = relocate_reactions(m, prop, {"r2": "m"}, biomass_reaction="bio")
    assert res.moved["r2"] == ("c", "m")
    # B (made by r1 in c, consumed by r2 in m) and C (made by r2 in m, consumed by bio in c) must bridge
    assert ("B", "m") in res.transports_added
    assert ("C", "m") in res.transports_added
    assert res.proposal.placements["r2"] == ["m"]
    assert res.certified and res.growth_after > 1e-9
    assert _grows(res.apply(m))


def test_orphaned_transport_at_old_compartment_is_removed():
    m = _chain3()
    prop = _proposal({"r1": ["c"], "r2": ["m"]}, transports=[("B", "m"), ("C", "m")])
    res = relocate_reactions(m, prop, {"r2": "c"}, biomass_reaction="bio")
    # r2 back in c -> nothing touches B or C in m any more -> both transports drop
    assert ("B", "m") in res.transports_removed
    assert ("C", "m") in res.transports_removed
    assert res.proposal.added_transports == []


def test_gene_sibling_is_co_moved():
    m = _chain3(shared_gene=True)  # r1 and r2 both catalysed by g1
    prop = _proposal({"r1": ["c"], "r2": ["c"]})
    res = relocate_reactions(m, prop, {"r1": "m"}, biomass_reaction="bio")
    assert "r2" in res.co_moved
    assert res.moved["r2"] == ("c", "m")
    assert res.proposal.placements["r2"] == ["m"]


def test_isozyme_shared_reaction_is_not_co_moved():
    m = _chain3(isozyme=True)  # r2 = g1 or g2; only g1 is moved, so r2 stays (ambiguous)
    prop = _proposal({"r1": ["c"], "r2": ["c"]})
    res = relocate_reactions(m, prop, {"r1": "m"}, biomass_reaction="bio")
    assert "r2" not in res.co_moved
    assert res.proposal.placements["r2"] == ["c"]


def test_gene_siblings_can_be_disabled():
    m = _chain3(shared_gene=True)
    prop = _proposal({"r1": ["c"], "r2": ["c"]})
    res = relocate_reactions(m, prop, {"r1": "m"}, biomass_reaction="bio", move_gene_siblings=False)
    assert res.co_moved == []
    assert res.proposal.placements["r2"] == ["c"]


def test_impermeant_transport_is_warned():
    m = cobra.Model("coa")
    A = _met("A_c")
    pal = _met("palCoA_c", name="palmitoyl-CoA")
    C = _met("C_c")
    m.add_metabolites([A, pal, C])
    ex = cobra.Reaction("EX_A", lower_bound=-10, upper_bound=1000)
    ex.add_metabolites({A: -1})
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({A: -1, pal: 1})
    r1.gene_reaction_rule = "g1"
    r2 = cobra.Reaction("r2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({pal: -1, C: 1})
    r2.gene_reaction_rule = "g2"
    bio = cobra.Reaction("bio", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({C: -1})
    m.add_reactions([ex, r1, r2, bio])
    m.objective = "bio"
    prop = _proposal({"r1": ["c"], "r2": ["c"]})
    res = relocate_reactions(m, prop, {"r2": "m"}, biomass_reaction="bio")
    assert ("palCoA", "m") in res.transports_added
    assert any("impermeant" in w for w in res.warnings)


def test_noop_when_already_in_target_compartment():
    m = _chain3()
    prop = _proposal({"r1": ["c"], "r2": ["c"]})
    res = relocate_reactions(m, prop, {"r1": "c"}, biomass_reaction="bio")
    assert res.moved == {}
    assert res.proposal is not prop  # a copy, never the aliased original
    assert res.proposal.placements == prop.placements
    assert res.growth_after > 1e-9  # growth is genuinely measured, not a false 0.0/not-certified


def test_unrelated_dual_localisation_is_preserved():
    # a reaction dual-localised in the input and untouched by the curator must keep BOTH compartments.
    m = _chain3()
    prop = _proposal({"r1": ["c", "m"], "r2": ["c"]})
    res = relocate_reactions(m, prop, {"r2": "m"}, biomass_reaction="bio")
    assert res.proposal.placements["r1"] == ["c", "m"]  # not collapsed to ['c']


def test_boundary_and_missing_reactions_are_skipped():
    m = _chain3()
    prop = _proposal({"r1": ["c"], "r2": ["c"]})
    res = relocate_reactions(m, prop, {"EX_A": "m", "does_not_exist": "m"}, biomass_reaction="bio")
    assert res.moved == {}
    assert any("boundary" in w for w in res.warnings)
    assert any("not in model" in w for w in res.warnings)


def test_summary_is_readable_ascii():
    m = _chain3()
    prop = _proposal({"r1": ["c"], "r2": ["c"]})
    res = relocate_reactions(m, prop, {"r2": "m"}, biomass_reaction="bio")
    s = res.summary()
    assert s.isascii()
    assert "r2" in s and "->" in s
