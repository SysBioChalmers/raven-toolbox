"""Tests for get_kegg_model_for_organism (KEGG organism-ID mode, step 3b.4)."""
from pathlib import Path

import cobra
import pandas as pd
import pytest

from ravengem.reconstruction.kegg import (
    build_kegg_tables,
    build_reference_model,
    get_kegg_model_for_organism,
    get_kegg_model_for_organism_from_artefacts,
    parse_kegg_compounds,
    parse_kegg_dump,
    parse_kegg_reactions,
)

DUMP = Path(__file__).parent / "data" / "kegg_dump"


@pytest.fixture(scope="module")
def artefacts():
    reactions = parse_kegg_reactions(DUMP)
    compounds = parse_kegg_compounds(DUMP)
    linked = {ko for r in reactions for ko in r.kos}
    from ravengem.reconstruction.kegg import parse_kegg_kos

    kos = parse_kegg_kos(DUMP, keep=linked)
    model = build_reference_model(reactions, compounds)
    tables = build_kegg_tables(reactions, kos)
    return model, tables


def _build(artefacts, organism_id, **kw):
    model, tables = artefacts
    return get_kegg_model_for_organism(
        organism_id,
        model,
        tables["ko_reaction"],
        tables["organism_gene_ko"],
        rxn_flags=tables["rxn_flags"],
        **kw,
    )


# --------------------------------------------------------------------------- #
# Core behaviour
# --------------------------------------------------------------------------- #
def test_eco_keeps_only_its_reactions(artefacts):
    # eco has b0001 -> K00002 -> R00100 only.
    model = _build(artefacts, "eco")
    assert {r.id for r in model.reactions} == {"R00100"}
    assert model.id == "eco"


def test_eco_gpr_and_gene_annotation(artefacts):
    model = _build(artefacts, "eco")
    r = model.reactions.get_by_id("R00100")
    assert r.gene_reaction_rule == "b0001"
    assert model.genes.get_by_id("b0001").annotation["kegg.genes"] == "eco:b0001"
    assert r.notes["note"].startswith("Included by get_kegg_model_for_organism")


def test_bsu_or_joins_multiple_genes(artefacts):
    # bsu has BSU31050 + BSU31060, both -> K01194 -> R00010.
    model = _build(artefacts, "bsu")
    r = model.reactions.get_by_id("R00010")
    assert set(r.genes) == {model.genes.get_by_id("BSU31050"), model.genes.get_by_id("BSU31060")}
    assert r.gene_reaction_rule == "BSU31050 or BSU31060"


def test_case_insensitive_organism(artefacts):
    assert "R00010" in _build(artefacts, "BSU").reactions


def test_orphan_metabolites_pruned(artefacts):
    # eco keeps only R00100 (C00002, C00003); trehalose/glucose mets should go.
    model = _build(artefacts, "eco")
    assert {m.id for m in model.metabolites} == {"C00002", "C00003"}


def test_reference_model_unmodified(artefacts):
    reference, _ = artefacts
    before = len(reference.reactions)
    _build(artefacts, "eco")
    assert len(reference.reactions) == before  # worked on a copy
    assert len(reference.genes) == 0


# --------------------------------------------------------------------------- #
# Spontaneous handling
# --------------------------------------------------------------------------- #
def test_spontaneous_reaction_kept_without_genes(artefacts):
    # R00100 is spontaneous; for bsu it has no genes but is kept (no GPR).
    model = _build(artefacts, "bsu", keep_spontaneous=True)
    assert "R00100" in model.reactions
    assert model.reactions.get_by_id("R00100").gene_reaction_rule == ""


def test_spontaneous_dropped_when_disabled(artefacts):
    model = _build(artefacts, "bsu", keep_spontaneous=False)
    assert "R00100" not in model.reactions
    assert "R00010" in model.reactions  # the gene-backed reaction stays


# --------------------------------------------------------------------------- #
# Quality filters take precedence over having genes
# --------------------------------------------------------------------------- #
def _tiny_general_case():
    ref = cobra.Model("KEGG")
    a = cobra.Metabolite("C1", compartment="s")
    b = cobra.Metabolite("C2", compartment="s")
    ref.add_metabolites([a, b])
    rxn = cobra.Reaction("R1")
    ref.add_reactions([rxn])
    rxn.add_metabolites({a: -1, b: 1})
    ko_reaction = pd.DataFrame([("K1", "R1")], columns=["ko", "reaction"])
    ogk = pd.DataFrame([("xyz", "g1", "K1")], columns=["organism", "gene", "ko"])
    flags = pd.DataFrame(
        [("R1", False, False, False, True)],
        columns=["reaction", "spontaneous", "undefined_stoich", "incomplete", "general"],
    )
    return ref, ko_reaction, ogk, flags


def test_general_filter_drops_reaction_with_genes():
    ref, ko_reaction, ogk, flags = _tiny_general_case()
    model = get_kegg_model_for_organism("xyz", ref, ko_reaction, ogk, rxn_flags=flags)
    assert "R1" not in model.reactions  # general + keep_general=False (default)


def test_general_kept_when_enabled():
    ref, ko_reaction, ogk, flags = _tiny_general_case()
    model = get_kegg_model_for_organism(
        "xyz", ref, ko_reaction, ogk, rxn_flags=flags, keep_general=True
    )
    assert model.reactions.get_by_id("R1").gene_reaction_rule == "g1"


# --------------------------------------------------------------------------- #
# Validation + artefact loading
# --------------------------------------------------------------------------- #
def test_unknown_organism_raises(artefacts):
    with pytest.raises(ValueError, match="no genes"):
        _build(artefacts, "zzz")


def test_domain_mode_not_implemented(artefacts):
    with pytest.raises(NotImplementedError, match="getPhylDist"):
        _build(artefacts, "eukaryotes")


def test_from_artefacts_roundtrip(tmp_path):
    parse_kegg_dump(DUMP, tmp_path)
    model = get_kegg_model_for_organism_from_artefacts("eco", tmp_path)
    assert {r.id for r in model.reactions} == {"R00100"}
    assert model.reactions.get_by_id("R00100").gene_reaction_rule == "b0001"
