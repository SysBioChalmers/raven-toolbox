"""Tests for compartment selectors (getMetsInComp / getRxnsInComp ports)."""
import cobra
import pytest

from ravengem.manipulation import add_reactions_from_equations
from ravengem.utils import get_metabolites_in_compartment, get_reactions_in_compartment


@pytest.fixture
def model():
    m = cobra.Model("t")
    m.add_metabolites(
        [
            cobra.Metabolite("a_c", compartment="c"),
            cobra.Metabolite("b_c", compartment="c"),
            cobra.Metabolite("a_m", compartment="m"),
        ]
    )
    add_reactions_from_equations(
        m,
        [
            {"id": "R_c", "equation": "a_c --> b_c"},   # fully in c
            {"id": "R_t", "equation": "a_c --> a_m"},   # transport c<->m
        ],
    )
    return m


def test_metabolites_in_compartment(model):
    assert {m.id for m in get_metabolites_in_compartment(model, "c")} == {"a_c", "b_c"}
    assert {m.id for m in get_metabolites_in_compartment(model, "m")} == {"a_m"}


def test_reactions_in_compartment_fully_contained(model):
    # default: exclude transport reactions touching other compartments
    assert {r.id for r in get_reactions_in_compartment(model, "c")} == {"R_c"}


def test_reactions_in_compartment_include_partial(model):
    assert {r.id for r in get_reactions_in_compartment(model, "c", include_partial=True)} == {
        "R_c",
        "R_t",
    }
    # m only appears in the transport reaction
    assert get_reactions_in_compartment(model, "m", include_partial=True)[0].id == "R_t"
    assert get_reactions_in_compartment(model, "m") == []  # none fully in m


def test_unknown_compartment_errors(model):
    with pytest.raises(ValueError, match="not in the model"):
        get_metabolites_in_compartment(model, "x")
    with pytest.raises(ValueError, match="not in the model"):
        get_reactions_in_compartment(model, "x")
