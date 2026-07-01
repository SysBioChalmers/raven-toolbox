"""Tests for evidence-aware transport scoring: annotation parsing + the per-metabolite cost."""
import pathlib

import cobra
import pandas as pd

from raven_toolbox.localization import (
    TransporterAnnotation,
    annotate_transporters,
    evidence_aware_transport_cost,
)


def _model():
    """Malate present in cytosol + mito (transportable); glucose in cytosol only."""
    m = cobra.Model("t")
    mets = [
        cobra.Metabolite("mal_c", name="malate", compartment="c"),
        cobra.Metabolite("mal_m", name="malate", compartment="m"),
        cobra.Metabolite("glc_c", name="glucose", compartment="c"),
    ]
    m.add_metabolites(mets)
    return m


SUBSTRATE = {"mal": {"organic_acid"}, "glc": {"sugar"}}


def _substrate_of(met):
    return SUBSTRATE.get(met.id.rsplit("_", 1)[0], set())


# --------------------------------------------------------------------------- annotate_transporters
def test_annotate_parses_string_lists_and_dedups():
    table = pd.DataFrame({
        "gene": ["CARR1", "CARR1", "CARR2"],
        "confidence": [0.4, 0.9, 0.7],
        "families": ["PF00153", "PF00153", "PF07690"],
        "substrate_classes": ["organic_acid", "organic_acid;nucleotide", "sugar"],
        "mechanism": [None, "antiport", "symport"],
    })
    ann = annotate_transporters(table)
    assert set(ann) == {"CARR1", "CARR2"}
    # CARR1 keeps the higher-confidence row and the union of substrate classes
    assert ann["CARR1"].confidence == 0.9
    assert ann["CARR1"].substrate_classes == frozenset({"organic_acid", "nucleotide"})
    assert ann["CARR1"].mechanism == "antiport"
    assert ann["CARR2"].substrate_classes == frozenset({"sugar"})


# --------------------------------------------------------------------- evidence_aware_transport_cost
def test_supported_transport_is_cheaper_unsupported_is_full():
    ann = {"CARR1": TransporterAnnotation("CARR1", 0.9,
                                          substrate_classes=frozenset({"organic_acid"}))}
    gene_comps = {"CARR1": {"m"}}  # a mito carrier
    costs = evidence_aware_transport_cost(
        _model(), ann, gene_comps, substrate_of=_substrate_of, base_cost=0.5)
    # malate: substrate matches + carrier sits at mito (a compartment malate occupies) -> supported
    assert costs["mal"] == 0.5 * (1 - 0.9)
    # glucose: no sugar carrier -> full prior
    assert costs["glc"] == 0.5


def test_substrate_mismatch_gets_no_discount():
    ann = {"CARR1": TransporterAnnotation("CARR1", 0.9,
                                          substrate_classes=frozenset({"amino_acid"}))}
    costs = evidence_aware_transport_cost(
        _model(), ann, {"CARR1": {"m"}}, substrate_of=_substrate_of, base_cost=0.5)
    assert costs["mal"] == 0.5  # organic_acid != amino_acid


def test_membrane_mismatch_gets_no_discount():
    # sugar carrier, but localised to the nucleus where glucose does not occur -> no support
    ann = {"CARR2": TransporterAnnotation("CARR2", 0.8, substrate_classes=frozenset({"sugar"}))}
    costs = evidence_aware_transport_cost(
        _model(), ann, {"CARR2": {"n"}}, substrate_of=_substrate_of, base_cost=0.5)
    assert costs["glc"] == 0.5


def test_output_is_a_valid_transport_cost_mapping():
    # every metabolite base is present, so the result is a self-contained transport_cost mapping
    ann = {"CARR1": TransporterAnnotation("CARR1", 1.0,
                                          substrate_classes=frozenset({"organic_acid"}))}
    costs = evidence_aware_transport_cost(
        _model(), ann, {"CARR1": {"m"}}, substrate_of=_substrate_of)
    assert set(costs) == {"mal", "glc"}
    assert costs["mal"] == 0.0 and costs["glc"] == 0.5


# --------------------------------------------------------------------------- coarse tables
def test_family_tables_use_valid_classes():
    from raven_toolbox.localization.transporter_tables import (
        COARSE_CLASSES,
        PFAM_TRANSPORTERS,
        TC_FAMILY_CLASS,
    )
    for _name, classes in PFAM_TRANSPORTERS.values():
        assert classes <= COARSE_CLASSES
    for classes in TC_FAMILY_CLASS.values():
        assert classes <= COARSE_CLASSES


# ----------------------------------------------------------- annotation back-end (integration)
def test_annotate_proteome_finds_yeast_transporters():
    """hmmsearch + diamond backend on a real yeast proteome. Skips offline / without the binaries."""
    import pytest

    from raven_toolbox.localization import annotate_proteome

    fasta = pathlib.Path("data/deeploc/yeast-GEM_proteins_001.fasta")
    if not fasta.exists():
        pytest.skip("yeast proteome FASTA not present")
    try:
        ann = annotate_proteome(fasta, threads=2)
    except Exception as exc:  # binaries / db download unavailable (e.g. offline CI)  # noqa: BLE001
        pytest.skip(f"transporter backend unavailable: {exc}")
    assert len(ann) > 20  # a genome-scale proteome has many transporters
    assert any("sugar" in a.substrate_classes for a in ann.values())        # HXT sugar transporters
    mcf = [a for a in ann.values() if "PF00153" in a.families]              # mitochondrial carriers
    assert mcf and any("carboxylate" in a.substrate_classes for a in mcf)
    assert any(a.substrate_chebi for a in ann.values())                     # curated TCDB substrate ChEBIs


# --------------------------------------------------------- metabolite -> coarse class (model side)
def test_default_substrate_of_classifies_key_metabolites():
    from raven_toolbox.localization import default_substrate_of

    def m(name):
        return cobra.Metabolite("x", name=name)

    assert "carboxylate" in default_substrate_of(m("(S)-malate"))
    assert "carboxylate" in default_substrate_of(m("citrate(3-)"))
    assert "sugar" in default_substrate_of(m("D-glucose"))
    assert "amino_acid" in default_substrate_of(m("L-glutamate"))
    assert {"nucleotide", "cofactor_vitamin"} <= default_substrate_of(m("NADPH"))
    assert default_substrate_of(m("some unidentifiable compound")) == frozenset()  # safe default


# ------------------------------------------------------- ChEBI substrate layer (specific matching)
def test_default_metabolite_chebi_normalises():
    from raven_toolbox.localization import default_metabolite_chebi

    m = cobra.Metabolite("x")
    m.annotation["chebi"] = "CHEBI:15589"
    assert default_metabolite_chebi(m) == frozenset({"CHEBI:15589"})
    m.annotation["chebi"] = ["12345", "CHEBI:67890"]  # bare + prefixed, list-valued
    assert default_metabolite_chebi(m) == frozenset({"CHEBI:12345", "CHEBI:67890"})
    assert default_metabolite_chebi(cobra.Metabolite("y")) == frozenset()  # un-annotated -> empty


def test_substrate_ontology_graded_match(tmp_path):
    import gzip

    from raven_toolbox.localization import SubstrateOntology

    rel = tmp_path / "rel.tsv.gz"
    with gzip.open(rel, "wt", encoding="utf-8") as fh:
        fh.write("CHEBI:2\tis_a\tCHEBI:1\n")                    # 2 is_a 1
        fh.write("CHEBI:3\tis_a\tCHEBI:2\n")                    # 3 is_a 2 is_a 1
        fh.write("CHEBI:5\tis_conjugate_base_of\tCHEBI:4\n")    # 5 <-> 4 (symmetric bridge)
        fh.write("CHEBI:99\talt_id\tCHEBI:3\n")                 # 99 is a secondary/deprecated id of 3
    sub = tmp_path / "sub.tsv"
    sub.write_text("2.A.1.1.1\tCHEBI:1;CHEBI:4\n", encoding="utf-8")
    onto = SubstrateOntology.load(relations_path=rel, substrates_path=sub)

    assert onto.substrates_of(["2.A.1.1.1"]) == frozenset({"CHEBI:1", "CHEBI:4"})
    assert onto.match(["CHEBI:1"], ["CHEBI:1"]) == 1.0          # exact
    assert onto.match(["CHEBI:2"], ["CHEBI:1"]) == 0.9          # 1 hop up
    assert onto.match(["CHEBI:3"], ["CHEBI:1"]) == 0.8          # 2 hops up (nearest wins)
    assert onto.match(["CHEBI:5"], ["CHEBI:4"]) == 0.9          # conjugate bridge is bidirectional
    assert onto.match(["CHEBI:99"], ["CHEBI:1"]) == 0.8         # secondary id normalised to 3 first
    assert onto.match(["CHEBI:3"], ["CHEBI:99"]) == 1.0         # normalised on the substrate side too
    assert onto.match(["CHEBI:9"], ["CHEBI:1"]) == 0.0          # unrelated
    assert onto.match(["CHEBI:1"], []) == 0.0                   # gene has no curated substrate


def test_evidence_uses_chebi_when_coarse_class_missing(tmp_path):
    import gzip

    from raven_toolbox.localization import SubstrateOntology

    rel = tmp_path / "rel.tsv.gz"
    with gzip.open(rel, "wt", encoding="utf-8") as fh:
        fh.write("CHEBI:3\tis_a\tCHEBI:2\n")
        fh.write("CHEBI:2\tis_a\tCHEBI:1\n")
    sub = tmp_path / "sub.tsv"
    sub.write_text("x\tCHEBI:1\n", encoding="utf-8")
    onto = SubstrateOntology.load(relations_path=rel, substrates_path=sub)

    model = cobra.Model("t")
    mets = []
    for comp in ("c", "m"):
        met = cobra.Metabolite(f"glc_{comp}", compartment=comp)
        met.annotation["chebi"] = "CHEBI:3"  # rolls up to the carrier substrate CHEBI:1 in 2 hops
        mets.append(met)
    rxn = cobra.Reaction("t_glc")
    rxn.add_metabolites({mets[0]: -1, mets[1]: 1})
    model.add_reactions([rxn])
    # carrier gene has NO coarse class, only a curated substrate ChEBI -> coarse misses, ChEBI catches it
    ann = {"G": TransporterAnnotation("G", 1.0, substrate_chebi=frozenset({"CHEBI:1"}))}
    costs = evidence_aware_transport_cost(
        model, ann, {"G": {"c"}}, substrate_of=lambda _m: [], ontology=onto, base_cost=0.5)
    assert costs["glc"] == 0.5 * (1.0 - 0.8)  # 2-hop ChEBI match -> weight 0.8
