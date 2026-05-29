"""Tests for raven_python.io.excel (exportToExcelFormat port, export only)."""
import cobra
import pytest

openpyxl = pytest.importorskip("openpyxl")

from raven_python.io import export_to_excel
from raven_python.manipulation import add_reactions_from_equations


@pytest.fixture
def model():
    m = cobra.Model("yeastGEM")
    m.name = "Yeast"
    m.compartments = {"c": "cytoplasm"}
    m.notes["metaData"] = {"taxonomy": "taxonomy/559292", "defaultLB": "-1000"}
    m.add_metabolites(
        [
            cobra.Metabolite("atp_c", name="ATP", formula="C10H16N5O13P3", charge=-4, compartment="c"),
            cobra.Metabolite("adp_c", name="ADP", compartment="c"),
        ]
    )
    m.metabolites.atp_c.annotation = {"kegg.compound": ["C00002"], "smiles": ["C1=NC"]}
    m.metabolites.atp_c.notes = {"inchis": "InChI=1S/X"}
    add_reactions_from_equations(
        m,
        [{"id": "R1", "equation": "atp_c <=> adp_c", "name": "rxn one",
          "gene_reaction_rule": "G1", "subsystem": "glycolysis"}],
    )
    r = m.reactions.R1
    r.annotation = {"ec-code": ["1.1.1.1"], "kegg.reaction": ["R00001"]}
    r.notes = {"confidence_score": 2, "note": "a note", "references": "PMID:1"}
    r.objective_coefficient = 1
    return m


def _wb(path):
    return openpyxl.load_workbook(path)


def test_sheets_present(model, tmp_path):
    out = tmp_path / "m.xlsx"
    export_to_excel(model, out)
    wb = _wb(out)
    assert set(wb.sheetnames) == {"RXNS", "METS", "COMPS", "GENES", "MODEL"}


def test_rxns_sheet(model, tmp_path):
    out = tmp_path / "m.xlsx"
    export_to_excel(model, out)
    ws = _wb(out)["RXNS"]
    header = [c.value for c in ws[1]]
    row = {header[i]: c.value for i, c in enumerate(ws[2])}
    assert row["ID"] == "R1"
    assert row["NAME"] == "rxn one"
    assert "ATP[c]" in row["EQUATION"] and "<=>" in row["EQUATION"]
    assert row["EC-NUMBER"] == "1.1.1.1"
    assert row["GENE ASSOCIATION"] == "G1"
    assert row["SUBSYSTEM"] == "glycolysis"
    assert row["OBJECTIVE"] == 1
    assert row["CONFIDENCE SCORE"] == 2
    assert row["NOTE"] == "a note"
    assert row["MIRIAM"] == "kegg.reaction/R00001"  # ec-code excluded (own column)


def test_mets_sheet(model, tmp_path):
    out = tmp_path / "m.xlsx"
    export_to_excel(model, out)
    ws = _wb(out)["METS"]
    header = [c.value for c in ws[1]]
    rows = {
        r[header.index("REPLACEMENT ID")].value: {header[i]: c.value for i, c in enumerate(r)}
        for r in ws.iter_rows(min_row=2)
    }
    atp = rows["atp_c"]
    assert atp["ID"] == "ATP[c]"
    assert atp["NAME"] == "ATP"
    assert atp["InChI"] == "InChI=1S/X"
    assert atp["COMPOSITION"] is None  # suppressed when InChI present
    assert atp["CHARGE"] == -4
    assert atp["MIRIAM"] == "kegg.compound/C00002"  # smiles excluded


def test_model_sheet(model, tmp_path):
    out = tmp_path / "m.xlsx"
    export_to_excel(model, out)
    ws = _wb(out)["MODEL"]
    header = [c.value for c in ws[1]]
    row = {header[i]: c.value for i, c in enumerate(ws[2])}
    assert row["ID"] == "yeastGEM"
    assert row["NAME"] == "Yeast"
    assert row["TAXONOMY"] == "taxonomy/559292"
    assert row["DEFAULT LOWER"] == "-1000"


def test_genes_sheet(model, tmp_path):
    out = tmp_path / "m.xlsx"
    export_to_excel(model, out)
    ws = _wb(out)["GENES"]
    header = [c.value for c in ws[1]]
    row = {header[i]: c.value for i, c in enumerate(ws[2])}
    assert row["NAME"] == "G1"


def test_no_genes_skips_sheet(tmp_path):
    m = cobra.Model("t")
    m.add_metabolites([cobra.Metabolite("a_c", compartment="c")])
    add_reactions_from_equations(m, [{"id": "R1", "equation": "a_c -->"}])
    out = tmp_path / "m.xlsx"
    export_to_excel(m, out)
    assert "GENES" not in _wb(out).sheetnames
