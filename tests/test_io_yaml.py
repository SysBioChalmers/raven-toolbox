"""Tests for raven_python.io.yaml against the RAVEN fa281a1 (cobra-native !!omap) schema."""
from pathlib import Path

import cobra
import pytest
from cobra.io.yaml import yaml as cobra_yaml

from raven_python.io import read_yaml_model, write_yaml_model

# A model laid out exactly as RAVEN writeYAMLmodel (fa281a1) emits: cobra-native
# structure, RAVEN-only fields as top-level per-entry keys, smiles/ec-code inside
# the annotation block, metaData provenance-only, id/name/version top-level.
RAVEN_DOC = {
    "metabolites": [
        {
            "id": "s_0001",
            "name": "ATP",
            "compartment": "c",
            "formula": "C10H16N5O13P3",
            "charge": -4,
            "inchis": "InChI=1S/CH4",
            "deltaG": 12.5,
            "notes": "a metabolite note",
            "metFrom": "KEGG",
            "annotation": {"kegg.compound": ["C00002"], "smiles": ["C1=NC2"]},
        },
        {"id": "s_0002", "name": "ADP", "compartment": "c"},
    ],
    "reactions": [
        {
            "id": "R1",
            "name": "rxn one",
            "metabolites": {"s_0001": -1, "s_0002": 1},
            "lower_bound": -1000.0,
            "upper_bound": 1000.0,
            "gene_reaction_rule": "G1",
            "subsystem": "glycolysis",
            "confidence_score": 2,
            "references": "PMID:123",
            "rxnFrom": "manual",
            "notes": "a reaction note",
            "deltaG": -5.0,
            "annotation": {"ec-code": ["1.1.1.1"]},
        }
    ],
    "genes": [
        {"id": "G1", "name": "gene one", "protein": "P12345", "annotation": {"uniprot": ["P12345"]}}
    ],
    "id": "testModel",
    "name": "Test Model",
    "compartments": {"c": "cytoplasm"},
    "version": "1.0",
    "metaData": {"date": "2026-05-23", "taxonomy": "taxonomy/559292", "defaultLB": "-1000"},
    "ec-rxns": [{"id": "R1", "kcat": 100.0}],
}


@pytest.fixture
def yaml_file(tmp_path) -> Path:
    p = tmp_path / "model.yml"
    with open(p, "w", encoding="utf-8") as fh:
        cobra_yaml.dump(RAVEN_DOC, fh)
    return p


def test_standard_content(yaml_file):
    model = read_yaml_model(yaml_file)
    assert model.id == "testModel"
    assert model.name == "Test Model"
    assert {m.id for m in model.metabolites} == {"s_0001", "s_0002"}
    r = model.reactions.get_by_id("R1")
    assert r.bounds == (-1000.0, 1000.0)
    assert r.subsystem == "glycolysis"
    assert r.gene_reaction_rule == "G1"


def test_annotation_owned_by_cobra(yaml_file):
    # smiles / ec-code / miriam live in the annotation block (cobra reads them)
    model = read_yaml_model(yaml_file)
    assert model.metabolites.get_by_id("s_0001").annotation["smiles"] == ["C1=NC2"]
    assert model.metabolites.get_by_id("s_0001").annotation["kegg.compound"] == ["C00002"]
    assert model.reactions.get_by_id("R1").annotation["ec-code"] == ["1.1.1.1"]
    assert model.genes.get_by_id("G1").annotation["uniprot"] == ["P12345"]


def test_raven_only_fields_captured(yaml_file):
    model = read_yaml_model(yaml_file)
    a = model.metabolites.get_by_id("s_0001")
    assert a.notes["inchis"] == "InChI=1S/CH4"
    assert a.notes["deltaG"] == 12.5
    assert a.notes["note"] == "a metabolite note"  # RAVEN metNotes string, no crash
    assert a.notes["metFrom"] == "KEGG"
    assert "smiles" not in a.notes  # smiles stays in annotation
    r = model.reactions.get_by_id("R1")
    assert r.notes["confidence_score"] == 2
    assert r.notes["references"] == "PMID:123"
    assert r.notes["rxnFrom"] == "manual"
    assert r.notes["note"] == "a reaction note"
    assert r.notes["deltaG"] == -5.0
    assert model.genes.get_by_id("G1").notes["protein"] == "P12345"


def test_model_level_extras(yaml_file):
    model = read_yaml_model(yaml_file)
    assert model.notes["metaData"]["taxonomy"] == "taxonomy/559292"
    assert model.notes["version"] == "1.0"
    assert model.notes["_yaml_sections"]["ec-rxns"][0]["kcat"] == 100.0


def test_round_trip(yaml_file, tmp_path):
    model = read_yaml_model(yaml_file)
    out = tmp_path / "out.yml"
    write_yaml_model(model, out)
    reloaded = read_yaml_model(out)

    assert reloaded.id == "testModel"
    assert reloaded.notes["version"] == "1.0"
    assert reloaded.notes["metaData"]["taxonomy"] == "taxonomy/559292"
    a = reloaded.metabolites.get_by_id("s_0001")
    assert a.notes["deltaG"] == 12.5
    assert a.notes["note"] == "a metabolite note"
    assert a.annotation["smiles"] == ["C1=NC2"]
    r = reloaded.reactions.get_by_id("R1")
    assert r.notes["confidence_score"] == 2
    assert reloaded.genes.get_by_id("G1").notes["protein"] == "P12345"
    assert reloaded.notes["_yaml_sections"]["ec-rxns"][0]["id"] == "R1"


def test_extra_notes_not_dropped_when_free_text_note_present(yaml_file, tmp_path):
    """An entry with both a RAVEN free-text note and an extra note keeps both on write."""
    model = read_yaml_model(yaml_file)
    a = model.metabolites.get_by_id("s_0001")
    a.notes["note"] = "free text"
    a.notes["custom"] = "extra value"  # a non-RAVEN note that must not be silently lost
    out = tmp_path / "out.yml"
    write_yaml_model(model, out)
    text = out.read_text()
    assert "extra value" in text  # the leftover note survives serialization


def test_gzipped_round_trip(yaml_file, tmp_path):
    # A .yml.gz path is transparently gzipped on write and read.
    model = read_yaml_model(yaml_file)
    out = tmp_path / "out.yml.gz"
    write_yaml_model(model, out)
    assert out.read_bytes()[:2] == b"\x1f\x8b"  # gzip magic
    reloaded = read_yaml_model(out)
    assert reloaded.id == "testModel"
    assert {m.id for m in reloaded.metabolites} == {"s_0001", "s_0002"}


def test_output_is_cobra_readable(yaml_file, tmp_path):
    # The written file must load with stock cobra (it's cobra's native format).
    model = read_yaml_model(yaml_file)
    out = tmp_path / "out.yml"
    write_yaml_model(model, out)
    cobra_model = cobra.io.load_yaml_model(str(out))
    assert cobra_model.id == "testModel"
    assert {m.id for m in cobra_model.metabolites} == {"s_0001", "s_0002"}
    # RAVEN-only fields land in cobra notes; smiles in annotation
    assert cobra_model.metabolites.get_by_id("s_0001").annotation["smiles"] == ["C1=NC2"]


def test_write_emits_raven_top_level_keys(yaml_file, tmp_path):
    model = read_yaml_model(yaml_file)
    out = tmp_path / "out.yml"
    write_yaml_model(model, out)
    text = out.read_text()
    # RAVEN-only fields are lifted back to top-level entry keys, not buried in notes
    assert "inchis:" in text
    assert "deltaG:" in text
    assert "confidence_score:" in text
    assert "metaData:" in text


def test_legacy_id_in_metadata(tmp_path):
    # Older RAVEN files nest id/name under metaData and have no top-level id.
    legacy = {
        "metabolites": [{"id": "a_c", "name": "A", "compartment": "c"}],
        "reactions": [],
        "genes": [],
        "compartments": {"c": "cyt"},
        "metaData": {"id": "legacyModel", "name": "Legacy"},
    }
    p = tmp_path / "legacy.yml"
    with open(p, "w", encoding="utf-8") as fh:
        cobra_yaml.dump(legacy, fh)
    model = read_yaml_model(p)
    assert model.id == "legacyModel"
    assert model.name == "Legacy"


# Optional smoke test against a real model file if present.
_YEAST = Path("/home/eduardk/github/GECKO/tutorials/full_ecModel/models/yeast-GEM.yml")


@pytest.mark.skipif(not _YEAST.exists(), reason="real yeast-GEM.yml not available")
def test_real_yeast_gem_loads():
    model = read_yaml_model(_YEAST)
    assert len(model.reactions) > 1000
    # legacy file: identity comes from metaData
    assert model.id
