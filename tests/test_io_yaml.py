"""Tests for ravengem.io.yaml (readYAMLmodel / writeYAMLmodel port)."""
import textwrap
from pathlib import Path

import cobra
import pytest

from ravengem.io import read_yaml_model, write_yaml_model

# A small RAVEN/Metabolic-Atlas-style YAML: metaData block, RAVEN-only per-entry
# fields (smiles/deltaG, confidence_score/references, protein) and a foreign
# GECKO ec-rxns section.
RAVEN_YAML = textwrap.dedent(
    """
    metaData:
      id: testModel
      name: Test Model
      taxonomy: taxonomy/559292
      defaultLB: "-1000"
    compartments:
      c: cytoplasm
    metabolites:
      - id: a_c
        name: A
        compartment: c
        formula: C6H12O6
        charge: 0
        smiles: C(C)O
        deltaG: 12.5
      - id: b_c
        name: B
        compartment: c
    reactions:
      - id: R1
        name: rxn one
        metabolites:
          a_c: -1
          b_c: 1
        lower_bound: -1000
        upper_bound: 1000
        gene_reaction_rule: G1
        subsystem: glyco
        confidence_score: 2
        references: "PMID:123"
    genes:
      - id: G1
        name: gene one
        protein: P12345
    ec-rxns:
      - id: R1
        kcat: 100
    """
)


@pytest.fixture
def yaml_file(tmp_path) -> Path:
    p = tmp_path / "model.yml"
    p.write_text(RAVEN_YAML, encoding="utf-8")
    return p


def test_reads_standard_model(yaml_file):
    model = read_yaml_model(yaml_file)
    assert len(model.metabolites) == 2
    assert len(model.reactions) == 1
    assert len(model.genes) == 1
    r = model.reactions.get_by_id("R1")
    assert r.bounds == (-1000, 1000)
    assert {m.id: r.get_coefficient(m.id) for m in r.metabolites} == {"a_c": -1, "b_c": 1}


def test_metadata_sets_identity_and_is_preserved(yaml_file):
    model = read_yaml_model(yaml_file)
    # cobra alone would leave model.id is None; metaData restores it
    assert model.id == "testModel"
    assert model.name == "Test Model"
    assert model.notes["metaData"]["taxonomy"] == "taxonomy/559292"
    assert model.notes["metaData"]["defaultLB"] == "-1000"


def test_raven_only_fields_routed_by_meaning(yaml_file):
    model = read_yaml_model(yaml_file)
    a = model.metabolites.get_by_id("a_c")
    # chemical identifiers go to annotation, not notes
    assert a.annotation["smiles"] == "C(C)O"
    assert "smiles" not in a.notes
    # non-standard numeric/provenance data goes to notes
    assert a.notes["deltaG"] == 12.5
    r = model.reactions.get_by_id("R1")
    assert r.notes["confidence_score"] == 2
    assert r.notes["references"] == "PMID:123"
    assert model.genes.get_by_id("G1").notes["protein"] == "P12345"


def test_inchis_routed_to_annotation_as_inchi(tmp_path):
    p = tmp_path / "m.yml"
    p.write_text(
        "compartments: {c: cyt}\n"
        "metabolites:\n"
        "  - id: x_c\n    name: X\n    compartment: c\n    inchis: 'InChI=1S/CH4/h1H4'\n"
        "reactions: []\ngenes: []\n",
        encoding="utf-8",
    )
    model = read_yaml_model(p)
    assert model.metabolites.get_by_id("x_c").annotation["inchi"] == "InChI=1S/CH4/h1H4"


def test_foreign_sections_preserved(yaml_file):
    model = read_yaml_model(yaml_file)
    assert "ec-rxns" in model.notes["_yaml_sections"]
    assert model.notes["_yaml_sections"]["ec-rxns"][0]["kcat"] == 100


def test_round_trip(yaml_file, tmp_path):
    model = read_yaml_model(yaml_file)
    out = tmp_path / "out.yml"
    write_yaml_model(model, out)
    reloaded = read_yaml_model(out)

    assert reloaded.id == "testModel"
    assert reloaded.name == "Test Model"
    assert reloaded.notes["metaData"]["taxonomy"] == "taxonomy/559292"
    a = reloaded.metabolites.get_by_id("a_c")
    assert a.notes["deltaG"] == 12.5
    assert a.annotation["smiles"] == "C(C)O"
    assert reloaded.reactions.get_by_id("R1").notes["confidence_score"] == 2
    assert reloaded.genes.get_by_id("G1").notes["protein"] == "P12345"
    assert reloaded.notes["_yaml_sections"]["ec-rxns"][0]["id"] == "R1"


def test_write_lifts_extras_to_top_level(yaml_file, tmp_path):
    # Confirm RAVEN-only fields are emitted as per-entry top-level keys, not buried in notes.
    model = read_yaml_model(yaml_file)
    out = tmp_path / "out.yml"
    write_yaml_model(model, out)
    text = out.read_text()
    assert "deltaG:" in text
    assert "confidence_score:" in text
    assert "metaData:" in text


# Optional smoke test against a real yeast-GEM file if present locally.
_YEAST = Path("/home/eduardk/github/GECKO/tutorials/full_ecModel/models/yeast-GEM.yml")


@pytest.mark.skipif(not _YEAST.exists(), reason="real yeast-GEM.yml not available")
def test_real_yeast_gem_preserves_identity_and_deltaG():
    model = read_yaml_model(_YEAST)
    # cobra.io.load_yaml_model gives model.id is None here; we restore it from metaData
    assert model.id == "yeastGEM_develop"
    assert model.notes["metaData"]["taxonomy"] == "taxonomy/559292"
    # a RAVEN-only field cobra would have dropped
    assert any("deltaG" in m.notes for m in model.metabolites)
