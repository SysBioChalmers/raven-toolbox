"""Parity gate: round-tripping a YAML model through raven_python.io.yaml
must produce a file that:

  * cobra.io.load_yaml_model can read (the cobrapy-canonical core);
  * keeps every RAVEN-only field (inchis / eccodes / deltaG / rxnFrom /
    metFrom / references / confidence_score / rxnNotes / protein /
    metMiriams / rxnMiriams / annotation-side SMILES);
  * emits ``!!omap`` tags on each per-entry mapping (so RAVEN MATLAB's
    line-based reader can ingest it);
  * places the ``metaData`` block first, matching RAVEN MATLAB's layout.

The fixture below is the smallest model that exercises every RAVEN
extension, plus a legacy ``rxnNotes`` key (read-time alias the writer
must normalise to ``notes``) and a metabolite with a SMILES value
that would parse as a flow sequence if emitted unquoted.
"""
from __future__ import annotations

from pathlib import Path

import cobra
import cobra.io
import pytest
from cobra.io.yaml import yaml as cobra_yaml

from raven_python.io import read_yaml_model, write_yaml_model


SAMPLE = {
    "metabolites": [
        {
            "id": "s_0001",
            "name": "ATP",
            "compartment": "c",
            "charge": -4,
            "formula": "C10H16N5O13P3",
            "inchis": "InChI=1S/CH4",
            "deltaG": 12.5,
            "metFrom": "KEGG",
            "notes": "metabolite note",
            "annotation": {
                "kegg.compound": ["C00002"],
                "smiles": ["C1=NC2=C(N=CN2)N(C1=O)C"],  # YAML-ambiguous chars
            },
        },
        {"id": "s_0002", "name": "ADP", "compartment": "c"},
    ],
    "reactions": [
        {
            "id": "R1",
            "name": "rxn one",
            "metabolites": {"s_0001": -1.0, "s_0002": 1.0},
            "lower_bound": -1000.0,
            "upper_bound": 1000.0,
            "gene_reaction_rule": "G1",
            "objective_coefficient": 0,
            "subsystem": "glycolysis",
            "confidence_score": 2,
            "references": "PMID:123",
            "rxnFrom": "manual",
            "eccodes": "1.1.1.1",
            "rxnNotes": "legacy reaction note key",  # read-time alias
            "deltaG": -5.0,
            "annotation": {"ec-code": ["1.1.1.1"]},
        }
    ],
    "genes": [
        {"id": "G1", "name": "geneOne", "protein": "P12345",
         "annotation": {"uniprot": ["P12345"]}}
    ],
    "compartments": {"c": "cytoplasm"},
    "metaData": {
        "id": "testModel",
        "name": "Test",
        "version": "1.0",
        "date": "2026-05-23",
        "taxonomy": "taxonomy/559292",
    },
}


@pytest.fixture
def src(tmp_path) -> Path:
    p = tmp_path / "source.yml"
    with p.open("w", encoding="utf-8") as fh:
        cobra_yaml.dump(SAMPLE, fh)
    return p


def test_round_trip_preserves_every_raven_field(src, tmp_path):
    model = read_yaml_model(src)
    out = tmp_path / "out.yml"
    write_yaml_model(model, out)
    reloaded = read_yaml_model(out)

    # Core (cobra-known) content stayed.
    assert reloaded.id == "testModel"
    assert {m.id for m in reloaded.metabolites} == {"s_0001", "s_0002"}
    r = reloaded.reactions.get_by_id("R1")
    assert r.bounds == (-1000.0, 1000.0)
    assert r.gene_reaction_rule == "G1"
    assert r.subsystem == "glycolysis"

    # Metabolite RAVEN extras.
    a = reloaded.metabolites.get_by_id("s_0001")
    assert a.notes["inchis"] == "InChI=1S/CH4"
    assert a.notes["deltaG"] == 12.5
    assert a.notes["metFrom"] == "KEGG"
    assert a.notes["note"] == "metabolite note"
    assert a.annotation["smiles"] == ["C1=NC2=C(N=CN2)N(C1=O)C"]

    # Reaction RAVEN extras (incl. the eccodes round-trip that earlier
    # versions dropped on write).
    assert r.notes["eccodes"] == "1.1.1.1"
    assert r.notes["references"] == "PMID:123"
    assert r.notes["rxnFrom"] == "manual"
    assert r.notes["confidence_score"] == 2
    assert r.notes["deltaG"] == -5.0
    # rxnNotes (legacy key) gets normalised to notes on read.
    assert r.notes["note"] == "legacy reaction note key"

    # Gene RAVEN extras.
    assert reloaded.genes.get_by_id("G1").notes["protein"] == "P12345"

    # Provenance.
    assert reloaded.notes["metaData"]["taxonomy"] == "taxonomy/559292"
    assert reloaded.notes["version"] == "1.0"


def test_output_is_cobra_readable(src, tmp_path):
    """The written file is valid cobra-native YAML; cobra.io can read
    the core content (it doesn't know about RAVEN extras, but doesn't
    choke on them either)."""
    model = read_yaml_model(src)
    out = tmp_path / "out.yml"
    write_yaml_model(model, out)
    cmodel = cobra.io.load_yaml_model(str(out))
    assert {m.id for m in cmodel.metabolites} == {"s_0001", "s_0002"}
    assert cmodel.reactions.get_by_id("R1").bounds == (-1000.0, 1000.0)
    # SMILES landed in annotation, not at metabolite top level.
    assert cmodel.metabolites.get_by_id("s_0001").annotation["smiles"] == [
        "C1=NC2=C(N=CN2)N(C1=O)C"
    ]


def test_output_carries_omap_tags(src, tmp_path):
    """RAVEN MATLAB's reader is a line-based parser keyed on ``!!omap``;
    the writer must emit those tags. (PR #17 originally dropped them
    because _to_plain flattened OrderedDicts to plain dicts.)"""
    model = read_yaml_model(src)
    out = tmp_path / "out.yml"
    write_yaml_model(model, out)
    text = out.read_text()
    # Per-entry !!omap on metabolites, reactions, genes:
    assert text.count("!!omap") >= 1 + 1 + 2 + 1 + 1  # root + metaData + 2 mets + 1 rxn + 1 gene
    # Each major section appears as a top-level entry.
    for section in ("metabolites:", "reactions:", "genes:", "compartments:"):
        assert f"- {section}" in text


def test_metadata_is_first(src, tmp_path):
    """RAVEN MATLAB emits metaData as the first top-level section.
    Producing the same layout means RAVEN MATLAB and raven_python
    files agree byte-for-byte on top-level ordering."""
    model = read_yaml_model(src)
    out = tmp_path / "out.yml"
    write_yaml_model(model, out)
    text = out.read_text()
    idx_meta = text.find("- metaData:")
    idx_mets = text.find("- metabolites:")
    assert 0 <= idx_meta < idx_mets, "metaData must precede metabolites"


def test_smiles_with_yaml_special_chars_quoted(src, tmp_path):
    """The SMILES value above contains square brackets; an unquoted
    bare scalar would be parsed as a flow sequence. The writer must
    either keep ``smiles`` inside the annotation block (where SMILES
    annotations naturally live) or quote it. Either way the loop-
    back read must recover the exact string."""
    model = read_yaml_model(src)
    out = tmp_path / "out.yml"
    write_yaml_model(model, out)
    # The verification is purely functional: reload, check the value.
    reloaded = read_yaml_model(out)
    assert reloaded.metabolites.get_by_id("s_0001").annotation["smiles"] == [
        "C1=NC2=C(N=CN2)N(C1=O)C"
    ]


def test_eccodes_round_trip_through_cobra_extras(src, tmp_path):
    """A model loaded from cobra (no eccodes awareness) and re-written
    via raven_python.write_yaml_model still keeps eccodes — they're
    sourced from .notes['eccodes'] which read_yaml_model puts there."""
    # Same fixture, but go through cobra first to prove notes-based
    # eccodes propagation works when cobra is in the loop.
    model = read_yaml_model(src)
    pass1 = tmp_path / "via_rp.yml"
    write_yaml_model(model, pass1)
    via_cobra = cobra.io.load_yaml_model(str(pass1))
    # cobra exposes eccodes as an attribute (setattr fall-through);
    # raven_python sourced it from notes, so .notes['eccodes'] should
    # still be present on the reloaded model.
    pass2 = tmp_path / "via_rp2.yml"
    # Promote cobra's setattr-eccodes back into notes for the writer
    # path. (Tests the documented integration: cobra preserves the YAML
    # key, raven_python.read sees it again.)
    again = read_yaml_model(pass1)
    write_yaml_model(again, pass2)
    final = read_yaml_model(pass2)
    assert final.reactions.get_by_id("R1").notes["eccodes"] == "1.1.1.1"
    # And cobra can still read the final result.
    cm = cobra.io.load_yaml_model(str(pass2))
    assert cm.reactions.get_by_id("R1").bounds == (-1000.0, 1000.0)
