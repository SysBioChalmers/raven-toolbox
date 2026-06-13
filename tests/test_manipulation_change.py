"""Tests for raven_toolbox.manipulation.change (changeRxns port)."""
import cobra
import pytest

from raven_toolbox.manipulation import add_reactions_from_equations, change_reaction_equations


@pytest.fixture
def model():
    m = cobra.Model("t")
    m.add_metabolites(
        [
            cobra.Metabolite("a_c", name="A", compartment="c"),
            cobra.Metabolite("b_c", name="B", compartment="c"),
            cobra.Metabolite("c_c", name="C", compartment="c"),
        ]
    )
    add_reactions_from_equations(
        m,
        [
            {
                "id": "R1",
                "equation": "a_c <=> b_c",
                "name": "first",
                "bounds": (-30, 70),
                "gene_reaction_rule": "G1 or G2",
                "subsystem": "sub",
            },
            {"id": "R2", "equation": "a_c --> c_c"},
        ],
    )
    return m


def test_changes_stoichiometry(model):
    (rxn,) = change_reaction_equations(model, {"R1": "a_c --> 2 c_c"})
    assert rxn.id == "R1"
    assert {m.id: rxn.get_coefficient(m.id) for m in rxn.metabolites} == {
        "a_c": -1.0,
        "c_c": 2.0,
    }


def test_preserves_other_fields(model):
    before = model.reactions.get_by_id("R1")
    name, bounds, subsystem = before.name, before.bounds, before.subsystem
    genes = {g.id for g in before.genes}

    change_reaction_equations(model, {"R1": "a_c --> c_c"})

    after = model.reactions.get_by_id("R1")
    assert after.name == name
    assert after.bounds == bounds  # bounds untouched, per RAVEN
    assert after.subsystem == subsystem
    assert {g.id for g in after.genes} == genes


def test_preserves_reaction_order(model):
    order_before = [r.id for r in model.reactions]
    change_reaction_equations(model, {"R1": "b_c --> c_c"})
    assert [r.id for r in model.reactions] == order_before


def test_bounds_not_changed_by_arrow(model):
    # R1 starts reversible (-30, 70); a --> arrow must NOT make it irreversible.
    change_reaction_equations(model, {"R1": "a_c --> b_c"})
    assert model.reactions.get_by_id("R1").bounds == (-30, 70)


def test_name_mode(model):
    (rxn,) = change_reaction_equations(
        model, {"R2": "A --> C"}, mets_by="name", compartment="c"
    )
    assert {m.id for m in rxn.metabolites} == {"a_c", "c_c"}


def test_can_introduce_new_met(model):
    change_reaction_equations(
        model, {"R2": "a_c --> d_c"}, compartment="c"
    )
    assert "d_c" in model.metabolites
    assert model.reactions.get_by_id("R2").get_coefficient("d_c") == 1.0


def test_unknown_reaction_errors(model):
    with pytest.raises(ValueError, match="not found"):
        change_reaction_equations(model, {"NOPE": "a_c --> b_c"})


def test_multiple_reactions(model):
    changed = change_reaction_equations(model, {"R1": "a_c --> c_c", "R2": "b_c --> c_c"})
    assert [r.id for r in changed] == ["R1", "R2"]
    assert model.reactions.get_by_id("R2").get_coefficient("b_c") == -1.0
