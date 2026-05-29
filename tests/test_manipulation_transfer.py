"""Tests for add_reactions_from_model (addRxnsGenesMets port)."""
import cobra
import pytest

from ravengem.manipulation import add_reactions_from_equations, add_reactions_from_model


@pytest.fixture
def draft():
    m = cobra.Model("draft")
    m.add_metabolites(
        [cobra.Metabolite("glc_c", name="Glucose", formula="C6H12O6", compartment="c")]
    )
    # an existing reaction so glc_c is in use and we have an id to test skipping
    add_reactions_from_equations(m, [{"id": "R_existing", "equation": "glc_c <=>"}])
    return m


@pytest.fixture
def source():
    m = cobra.Model("source")
    m.add_metabolites(
        [
            # same name[comp] as draft's glc_c but a DIFFERENT id
            cobra.Metabolite("glucose_c", name="Glucose", formula="C6H12O6", compartment="c"),
            cobra.Metabolite("atp_c", name="ATP", formula="C10H16N5O13P3", charge=-4, compartment="c"),
            cobra.Metabolite("g6p_c", name="G6P", formula="C6H13O9P", compartment="c"),
        ]
    )
    add_reactions_from_equations(
        m,
        [
            {
                "id": "HEX",
                "equation": "glucose_c + atp_c --> g6p_c",
                "name": "hexokinase",
                "bounds": (0, 1000),
                "gene_reaction_rule": "G1",
                "subsystem": "glycolysis",
            },
            {"id": "R_existing", "equation": "glucose_c <=>"},  # id already in draft
        ],
    )
    return m


def test_metabolite_matched_by_name_comp_not_id(draft, source):
    add_reactions_from_model(draft, source, "HEX")
    hex_rxn = draft.reactions.get_by_id("HEX")
    # Glucose reused from the draft (id glc_c), NOT the source's glucose_c
    assert "glc_c" in {m.id for m in hex_rxn.metabolites}
    assert "glucose_c" not in draft.metabolites


def test_new_metabolites_added_with_metadata(draft, source):
    add_reactions_from_model(draft, source, "HEX")
    assert "atp_c" in draft.metabolites and "g6p_c" in draft.metabolites
    assert draft.metabolites.get_by_id("g6p_c").formula == "C6H13O9P"
    assert draft.metabolites.get_by_id("atp_c").charge == -4


def test_reaction_copied_with_bounds_and_name(draft, source):
    (rxn,) = add_reactions_from_model(draft, source, "HEX")
    assert rxn.id == "HEX"
    assert rxn.name == "hexokinase"
    assert rxn.bounds == (0, 1000)
    assert rxn.subsystem == "glycolysis"
    assert {m.id: rxn.get_coefficient(m.id) for m in rxn.metabolites} == {
        "glc_c": -1.0,
        "atp_c": -1.0,
        "g6p_c": 1.0,
    }


def test_genes_true_copies_gpr_and_creates_genes(draft, source):
    add_reactions_from_model(draft, source, "HEX", genes=True)
    assert draft.reactions.get_by_id("HEX").gene_reaction_rule == "G1"
    assert "G1" in draft.genes


def test_genes_false_no_gpr(draft, source):
    add_reactions_from_model(draft, source, "HEX", genes=False)
    assert draft.reactions.get_by_id("HEX").gene_reaction_rule == ""


def test_genes_string_override(draft, source):
    add_reactions_from_model(draft, source, "HEX", genes="G9 or G10")
    assert draft.reactions.get_by_id("HEX").gene_reaction_rule == "G9 or G10"


def test_skips_already_present(draft, source):
    added = add_reactions_from_model(draft, source, ["HEX", "R_existing"])
    assert [r.id for r in added] == ["HEX"]


def test_all_present_raises(draft, source):
    with pytest.raises(ValueError, match="already in the model"):
        add_reactions_from_model(draft, source, "R_existing")


def test_unknown_source_reaction_raises(draft, source):
    with pytest.raises(ValueError, match="not found in the source model"):
        add_reactions_from_model(draft, source, "NOPE")


def test_note_and_confidence_stored(draft, source):
    (rxn,) = add_reactions_from_model(draft, source, "HEX", note="from KEGG", confidence=2)
    assert rxn.notes["note"] == "from KEGG"
    assert rxn.notes["confidence_score"] == 2


# --- regression: intra-batch met-id minting collision (known_issues.md A3) ---

def test_intra_batch_id_minting_unique():
    """Two source mets whose ids both collide with the draft and whose name[comp]
    differs both get routed through new-id minting. The fix tracks ids minted in
    the current batch so the two don't collapse to the same generated id."""
    draft = cobra.Model("draft")
    draft.add_metabolites([
        cobra.Metabolite("atp_c", name="ATP-draft", compartment="c"),
        cobra.Metabolite("adp_c", name="ADP-draft", compartment="c"),
    ])
    source = cobra.Model("source")
    source.add_metabolites([
        cobra.Metabolite("atp_c", name="ATP-source", compartment="c"),
        cobra.Metabolite("adp_c", name="ADP-source", compartment="c"),
    ])
    rxn = cobra.Reaction("R1", lower_bound=0, upper_bound=1000)
    source.add_reactions([rxn])
    rxn.add_metabolites({
        source.metabolites.get_by_id("atp_c"): -1,
        source.metabolites.get_by_id("adp_c"): 1,
    })
    add_reactions_from_model(draft, source, "R1")
    # Both source mets minted distinct ids (m1 and m2) — not a collision.
    new_ids = sorted(m.id for m in draft.metabolites if m.id not in ("atp_c", "adp_c"))
    assert len(new_ids) == 2 and len(set(new_ids)) == 2
