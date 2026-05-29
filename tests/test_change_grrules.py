"""Tests for change_gene_reaction_rules (changeGrRules port)."""
import cobra
import pytest

from raven_python.manipulation import add_reactions_from_equations, change_gene_reaction_rules


@pytest.fixture
def model():
    m = cobra.Model("t")
    m.add_metabolites(
        [cobra.Metabolite("a_c", compartment="c"), cobra.Metabolite("b_c", compartment="c")]
    )
    add_reactions_from_equations(
        m,
        [
            {"id": "R1", "equation": "a_c --> b_c", "gene_reaction_rule": "G1"},
            {"id": "R2", "equation": "a_c --> b_c"},
        ],
    )
    return m


def test_replace_rule_and_create_genes(model):
    (rxn,) = change_gene_reaction_rules(model, {"R1": "G2 and G3"})
    assert rxn.gene_reaction_rule == "G2 and G3"
    assert {g.id for g in rxn.genes} == {"G2", "G3"}
    assert {"G2", "G3"} <= {g.id for g in model.genes}


def test_append_rule(model):
    change_gene_reaction_rules(model, {"R1": "G4"}, replace=False)
    # (G1) or (G4), normalised by cobra
    assert model.reactions.get_by_id("R1").gene_reaction_rule == "G1 or G4"


def test_append_when_empty_is_just_new(model):
    change_gene_reaction_rules(model, {"R2": "G5"}, replace=False)
    assert model.reactions.get_by_id("R2").gene_reaction_rule == "G5"


def test_batch(model):
    changed = change_gene_reaction_rules(model, {"R1": "GA", "R2": "GB"})
    assert [r.id for r in changed] == ["R1", "R2"]


def test_unknown_reaction_errors(model):
    with pytest.raises(ValueError, match="not found"):
        change_gene_reaction_rules(model, {"NOPE": "G1"})
