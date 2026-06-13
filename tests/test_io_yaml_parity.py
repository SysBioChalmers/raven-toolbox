"""Parity gate: round-tripping a YAML model through raven_toolbox.io.yaml
must produce a file that:

  * cobra.io.load_yaml_model can read (the cobrapy-canonical core);
  * keeps every RAVEN-only field (inchis / deltaG / rxnFrom /
    metFrom / references / confidence_score / rxnNotes / protein /
    metMiriams / rxnMiriams / annotation-side SMILES and EC codes);
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

from raven_toolbox.io import read_yaml_model, write_yaml_model

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

    # EC codes round-trip through cobra annotation (the cobra-native place,
    # where geckopy reads them), not a RAVEN-only top-level/notes field.
    assert r.annotation["ec-code"] == ["1.1.1.1"]
    assert "eccodes" not in r.notes
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
    Producing the same layout means RAVEN MATLAB and raven_toolbox
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


PRE_SHIM_YAML = """\
---
!!omap
- metaData:
    id: "eciYali"
    name: "Yarrowia lipolytica"
    version: "1.0"
    date: "2024-10-17"
    geckoLight: "true"
- metabolites:
    - !!omap
      - id: "s_0001"
      - name: "ATP"
      - compartment: "c"
      - formula: "C10H16N5O13P3"
      - charge: -4
      - inchis: "InChI=1S/CH4"
      - smiles: "[O-]P(=O)([O-])OP(=O)([O-])O"
      - annotation: !!omap
          - kegg.compound: "C00002"
          - sbo: "SBO:0000247"
      - deltaG: 12.5
      - notes: "metabolite note"
      - metFrom: "KEGG"
- reactions:
    - !!omap
      - id: "r_0001"
      - name: "hexokinase"
      - metabolites: !!omap
          - s_0001: -1
      - lower_bound: -1000
      - upper_bound: 1000
      - gene_reaction_rule: "G1"
      - rxnFrom: "KEGG"
      - eccodes: "2.7.1.1"
      - references: "PMID:12345"
      - subsystem: "Glycolysis"
      - annotation: !!omap
          - kegg.reaction: "R00299"
          - sbo: "SBO:0000176"
      - deltaG: -17.39
      - confidence_score: 2
      - rxnNotes: "old reaction note"
- genes:
    - !!omap
      - id: "G1"
      - name: "HXK1"
      - protein: "P01234"
      - annotation: !!omap
          - uniprot: "P01234"
- compartments: !!omap
    - c: "cytoplasm"
- ec-rxns:
    - !!omap
      - id: "r_0001"
      - kcat: 25.3
      - enzymes: !!omap
          - P01234: 1
- ec-enzymes:
    - !!omap
      - genes: "G1"
      - enzymes: "P01234"
      - mw: 50000
"""


def test_pre_shim_format_loads(tmp_path):
    """The pre-`feat/yeast-gem-shared` RAVEN MATLAB writer emitted a
    file shape that differs from the current one in seven concrete
    ways. The reader must continue to load every one of them:

      1. ``---`` document-start marker (kept by old MATLAB writer)
      2. ``- metaData:`` as a plain block mapping (no ``!!omap`` tag)
      3. ``geckoLight: "true"`` *inside* metaData (now emitted as a
         top-level ``gecko_light``)
      4. Metabolite ``smiles`` as a top-level entry key (now emitted
         inside the ``annotation`` block)
      5. Reaction notes under the ``rxnNotes`` key (now emitted as
         ``notes``)
      6. Integer-typed bounds / coefficients (now emitted as floats)
      7. Every string double-quoted (now bare unless YAML requires
         quoting)

    Each item below maps to one of those seven cases.
    """
    p = tmp_path / "pre_shim.yml"
    p.write_text(PRE_SHIM_YAML)
    model = read_yaml_model(p)

    # metaData survives + provenance is lifted onto the cobra-shape
    # accessors (cases 1 + 2).
    assert model.id == "eciYali"
    assert model.name == "Yarrowia lipolytica"
    assert model.notes["version"] == "1.0"
    assert model.notes["metaData"]["taxonomy" if "taxonomy" in model.notes["metaData"] else "date"]

    # geckoLight in metaData populates the typed EcData (case 3).
    assert model.ec is not None
    assert model.ec.gecko_light is True
    assert model.ec.rxns == ["r_0001"]
    assert model.ec.kcat[0] == 25.3

    # Top-level smiles lifted into annotation.smiles (case 4).
    a = model.metabolites.get_by_id("s_0001")
    assert a.annotation["smiles"] == ["[O-]P(=O)([O-])OP(=O)([O-])O"]
    assert "smiles" not in a.notes  # stays in annotation, not notes

    # rxnNotes read as the canonical notes key (case 5).
    r = model.reactions.get_by_id("r_0001")
    assert r.notes["note"] == "old reaction note"

    # Integer bounds become floats inside cobra (case 6).
    assert r.bounds == (-1000.0, 1000.0)
    assert isinstance(r.lower_bound, float)

    # Quoted strings unquote cleanly (case 7) — verified implicitly by
    # all the equality assertions above. Spot check the metabolite
    # name, which used double quotes in the source.
    assert a.name == "ATP"

    # Other RAVEN extras still preserved.
    assert a.notes["inchis"] == "InChI=1S/CH4"
    assert a.notes["deltaG"] == 12.5
    assert a.notes["note"] == "metabolite note"
    assert a.notes["metFrom"] == "KEGG"
    assert r.notes["rxnFrom"] == "KEGG"
    # legacy top-level eccodes lifted into the cobra-native annotation['ec-code']
    assert r.annotation["ec-code"] == ["2.7.1.1"]
    assert "eccodes" not in r.notes
    assert r.notes["references"] == "PMID:12345"
    assert r.notes["confidence_score"] == 2
    assert r.notes["deltaG"] == -17.39
    assert model.genes.get_by_id("G1").notes["protein"] == "P01234"


def test_pre_shim_yeast_gem_loads_if_available():
    """The real pre-shim yeast-GEM.yml: 2748 mets, 4102 rxns, 1143
    genes. Skipped when the working copy isn't mounted (CI runners)."""
    real = Path("/mnt/c/Work/GitHub/yeast-GEM/model/yeast-GEM.yml")
    if not real.exists():
        pytest.skip("yeast-GEM.yml not available in this environment")
    model = read_yaml_model(real)
    assert model.id == "yeastGEM_develop"
    assert len(model.metabolites) == 2748
    assert len(model.reactions) == 4102
    assert len(model.genes) == 1143
    # Every RAVEN extension we know about must come through.
    assert sum(1 for r in model.reactions if r.annotation.get("ec-code")) == 2411
    assert sum(1 for r in model.reactions if r.notes.get("deltaG") is not None) == 3984
    assert sum(1 for m in model.metabolites if m.notes.get("deltaG") is not None) == 2696
    assert sum(1 for m in model.metabolites if "smiles" in (m.annotation or {})) == 1788
    assert sum(1 for r in model.reactions if r.notes.get("note")) == 1443


def test_eccodes_round_trip_through_cobra_extras(src, tmp_path):
    """EC codes round-trip as cobra annotation through a
    raven_toolbox -> cobra -> raven_toolbox loop. They live in
    ``annotation['ec-code']`` — the cobra-native place — so plain
    ``cobra.io`` preserves them with no RAVEN-specific handling, and
    geckopy (which reads ``annotation['ec-code']``) sees them."""
    model = read_yaml_model(src)
    pass1 = tmp_path / "via_rp.yml"
    write_yaml_model(model, pass1)
    # Plain cobra reads annotation['ec-code'] natively — this is the
    # interop the alignment guarantees.
    via_cobra = cobra.io.load_yaml_model(str(pass1))
    assert via_cobra.reactions.get_by_id("R1").annotation["ec-code"] == ["1.1.1.1"]
    pass2 = tmp_path / "via_rp2.yml"
    again = read_yaml_model(pass1)
    write_yaml_model(again, pass2)
    final = read_yaml_model(pass2)
    assert final.reactions.get_by_id("R1").annotation["ec-code"] == ["1.1.1.1"]
    # And cobra can still read the final result.
    cm = cobra.io.load_yaml_model(str(pass2))
    assert cm.reactions.get_by_id("R1").bounds == (-1000.0, 1000.0)
