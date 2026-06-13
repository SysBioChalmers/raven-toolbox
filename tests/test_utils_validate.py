"""Tests for check_model (the surviving checks of checkModelStruct)."""
import cobra
import pytest

from raven_toolbox.manipulation import add_reactions_from_equations
from raven_toolbox.utils import ModelIssue, check_model


def _categories(issues, category):
    return [i.object_id for i in issues if i.category == category]


@pytest.fixture
def model():
    m = cobra.Model("t")
    m.add_metabolites(
        [
            cobra.Metabolite("a_c", name="A", compartment="c"),
            cobra.Metabolite("b_c", name="B", compartment="c"),
        ]
    )
    add_reactions_from_equations(
        m, [{"id": "R1", "equation": "a_c --> b_c", "gene_reaction_rule": "G1"}]
    )
    m.reactions.get_by_id("R1").objective_coefficient = 1
    return m


def test_clean_model_has_no_issues(model):
    assert check_model(model) == []


def test_orphan_metabolite(model):
    model.add_metabolites([cobra.Metabolite("orphan_c", name="Orphan", compartment="c")])
    assert "orphan_c" in _categories(check_model(model), "orphan_metabolite")


def test_orphan_gene(model):
    model.genes.append(cobra.core.gene.Gene("G_lonely"))
    assert "G_lonely" in _categories(check_model(model), "orphan_gene")


def test_empty_reaction(model):
    model.add_reactions([cobra.Reaction("R_empty")])
    assert "R_empty" in _categories(check_model(model), "empty_reaction")


def test_empty_metabolite_name(model):
    model.add_metabolites([cobra.Metabolite("noname_c", compartment="c")])
    # also an orphan, but we check the name category specifically
    assert "noname_c" in _categories(check_model(model), "empty_metabolite_name")


def test_duplicate_name_compartment(model):
    # second metabolite named "A" in compartment c
    dup = cobra.Metabolite("a2_c", name="A", compartment="c")
    model.add_metabolites([dup])
    model.reactions.get_by_id("R1").add_metabolites({dup: -1})  # keep it used
    issues = [i for i in check_model(model) if i.category == "duplicate_name_compartment"]
    assert len(issues) == 1
    assert "a_c" in issues[0].message and "a2_c" in issues[0].message


def test_no_objective(model):
    model.reactions.get_by_id("R1").objective_coefficient = 0
    cats = [i.category for i in check_model(model)]
    assert "objective" in cats


def test_multiple_objectives(model):
    add_reactions_from_equations(model, [{"id": "R2", "equation": "b_c --> a_c"}])
    model.reactions.get_by_id("R2").objective_coefficient = 1
    obj_issues = [i for i in check_model(model) if i.category == "objective"]
    assert len(obj_issues) == 1
    assert "Multiple" in obj_issues[0].message


def test_returns_model_issue_instances(model):
    model.add_reactions([cobra.Reaction("R_empty")])
    assert all(isinstance(i, ModelIssue) for i in check_model(model))
