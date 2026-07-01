"""Tests for evidence-aware transport scoring: annotation parsing + the per-metabolite cost."""
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
