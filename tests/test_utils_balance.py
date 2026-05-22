"""Tests for get_elemental_balance (getElementalBalance port)."""
import cobra
import pytest

from ravengem.utils import ElementalBalance, get_elemental_balance


@pytest.fixture
def model():
    m = cobra.Model("t")
    m.add_metabolites(
        [
            cobra.Metabolite("a_c", formula="C6H12O6", charge=0, compartment="c"),
            cobra.Metabolite("b_c", formula="C6H12O6", charge=0, compartment="c"),
            cobra.Metabolite("c_c", formula="C3H6O3", charge=0, compartment="c"),
            cobra.Metabolite("n_c", compartment="c"),  # no formula
        ]
    )
    r_bal = cobra.Reaction("R_bal"); m.add_reactions([r_bal])
    r_bal.build_reaction_from_string("a_c --> b_c")        # C6H12O6 -> C6H12O6
    r_unbal = cobra.Reaction("R_unbal"); m.add_reactions([r_unbal])
    r_unbal.build_reaction_from_string("a_c --> c_c")      # C6H12O6 -> C3H6O3
    r_unknown = cobra.Reaction("R_unknown"); m.add_reactions([r_unknown])
    r_unknown.build_reaction_from_string("a_c --> n_c")    # n_c has no formula
    return m


def test_balanced(model):
    (res,) = get_elemental_balance(model, ["R_bal"])
    assert res == ElementalBalance("R_bal", "balanced", {})


def test_unbalanced_reports_imbalance(model):
    (res,) = get_elemental_balance(model, ["R_unbal"])
    assert res.status == "unbalanced"
    # products - reactants: C3H6O3 - C6H12O6 = -C3H6O3
    assert res.imbalance == {"C": -3.0, "H": -6.0, "O": -3.0}


def test_missing_formula_is_unknown_not_silently_wrong(model):
    # cobra's check_mass_balance alone would silently report an imbalance here;
    # we flag it as unknown instead.
    (res,) = get_elemental_balance(model, ["R_unknown"])
    assert res.status == "unknown"
    assert res.imbalance == {}


def test_all_reactions_default(model):
    results = get_elemental_balance(model)
    assert {r.reaction_id: r.status for r in results} == {
        "R_bal": "balanced",
        "R_unbal": "unbalanced",
        "R_unknown": "unknown",
    }


def test_charge_excluded(model):
    # give a charge imbalance but keep elements balanced -> still "balanced"
    model.metabolites.get_by_id("b_c").charge = 1
    (res,) = get_elemental_balance(model, ["R_bal"])
    assert res.status == "balanced"
