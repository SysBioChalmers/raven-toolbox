"""Tests for merge_models (mergeModels port)."""
import cobra
import pytest

from raven_toolbox.manipulation import add_reactions_from_equations, merge_models


def _model(mid, mets, reactions):
    m = cobra.Model(mid)
    m.add_metabolites(mets)
    add_reactions_from_equations(m, reactions)
    return m


@pytest.fixture
def model_a():
    return _model(
        "A",
        [
            cobra.Metabolite("glc_c", name="Glucose", compartment="c"),
            cobra.Metabolite("g6p_c", name="G6P", compartment="c"),
        ],
        [{"id": "HEX", "equation": "glc_c --> g6p_c", "gene_reaction_rule": "GA"}],
    )


@pytest.fixture
def model_b():
    # same Glucose[c] compound but a DIFFERENT id
    return _model(
        "B",
        [
            cobra.Metabolite("glucose_c", name="Glucose", compartment="c"),
            cobra.Metabolite("lac_c", name="Lactate", compartment="c"),
        ],
        [{"id": "LDH", "equation": "glucose_c --> lac_c", "gene_reaction_rule": "GB"}],
    )


def test_unifies_metabolites_by_name_comp(model_a, model_b):
    merged = merge_models([model_a, model_b])
    glucoses = [m for m in merged.metabolites if m.name == "Glucose" and m.compartment == "c"]
    assert len(glucoses) == 1  # glc_c and glucose_c unified
    # both reactions reference the same merged Glucose object
    hex_glc = [m for m in merged.reactions.get_by_id("HEX").metabolites if m.name == "Glucose"][0]
    ldh_glc = [m for m in merged.reactions.get_by_id("LDH").metabolites if m.name == "Glucose"][0]
    assert hex_glc is ldh_glc


def test_match_by_id_keeps_distinct(model_a, model_b):
    merged = merge_models([model_a, model_b], match_by="id")
    glucoses = [m for m in merged.metabolites if m.name == "Glucose"]
    assert len(glucoses) == 2  # glc_c and glucose_c are distinct by id


def test_all_reactions_kept(model_a, model_b):
    merged = merge_models([model_a, model_b])
    assert {"HEX", "LDH"} <= {r.id for r in merged.reactions}


def test_reaction_id_collision_renamed(model_a):
    # two models with the same reaction id but different chemistry
    other = _model(
        "B",
        [cobra.Metabolite("glc_c", name="Glucose", compartment="c"),
         cobra.Metabolite("x_c", name="X", compartment="c")],
        [{"id": "HEX", "equation": "glc_c --> x_c"}],
    )
    merged = merge_models([model_a, other])
    assert "HEX" in {r.id for r in merged.reactions}
    assert "HEX_B" in {r.id for r in merged.reactions}  # renamed with source id


def test_genes_merged(model_a, model_b):
    merged = merge_models([model_a, model_b])
    assert {"GA", "GB"} <= {g.id for g in merged.genes}


def test_provenance_recorded(model_a, model_b):
    merged = merge_models([model_a, model_b])
    assert merged.reactions.get_by_id("HEX").notes["origin"] == "A"
    assert merged.reactions.get_by_id("LDH").notes["origin"] == "B"
    assert merged.genes.get_by_id("GA").notes["origin"] == "A"


def test_compartments_preserved(model_a):
    model_a.compartments = {"c": "cytoplasm"}
    merged = merge_models([model_a, model_a.copy()])
    assert merged.compartments.get("c") == "cytoplasm"


def test_single_model_returns_copy(model_a):
    merged = merge_models([model_a])
    assert merged is not model_a
    assert {r.id for r in merged.reactions} == {r.id for r in model_a.reactions}


def test_three_models(model_a, model_b):
    c = _model("C", [cobra.Metabolite("co2_c", name="CO2", compartment="c")],
               [{"id": "SINK", "equation": "co2_c -->"}])
    merged = merge_models([model_a, model_b, c])
    assert {"HEX", "LDH", "SINK"} <= {r.id for r in merged.reactions}


def test_bad_match_by(model_a, model_b):
    with pytest.raises(ValueError, match="match_by"):
        merge_models([model_a, model_b], match_by="oops")


# --- regression: formula/charge conflict ------------------------------------

def test_formula_conflict_warns():
    """Two models sharing a name[comp] but with different formulas warn instead
    of silently keeping the first."""
    a = _model("A",
        [cobra.Metabolite("g1", name="Glucose", formula="C6H12O6", compartment="c")],
        [{"id": "EX_A", "equation": "g1 -->"}])
    b = _model("B",
        [cobra.Metabolite("g2", name="Glucose", formula="C6H12O7", compartment="c")],
        [{"id": "EX_B", "equation": "g2 -->"}])
    with pytest.warns(UserWarning, match="different formulas"):
        merged = merge_models([a, b])
    # The merge still picks the first-seen — the test asserts the warning fired
    # and the model survives.
    assert "EX_A" in merged.reactions and "EX_B" in merged.reactions


def test_charge_conflict_warns():
    a = _model("A",
        [cobra.Metabolite("g1", name="Glucose", formula="C6H12O6", charge=0, compartment="c")],
        [{"id": "EX_A", "equation": "g1 -->"}])
    b = _model("B",
        [cobra.Metabolite("g2", name="Glucose", formula="C6H12O6", charge=-1, compartment="c")],
        [{"id": "EX_B", "equation": "g2 -->"}])
    with pytest.warns(UserWarning, match="different charges"):
        merge_models([a, b])
