"""Tests for get_kegg_model_for_organism (KEGG organism-ID mode, step 3b.4).

The ``kegg_dump`` fixture (tests/conftest.py) is a small, fully fictional dump —
no real KEGG content is committed.
"""
import cobra
import pandas as pd
import pytest

from raven_toolbox.reconstruction.kegg import (
    build_kegg_tables,
    build_reference_model,
    get_kegg_model_for_organism,
    get_kegg_model_for_organism_from_artefacts,
    parse_kegg_compounds,
    parse_kegg_dump,
    parse_kegg_reactions,
)


@pytest.fixture(scope="module")
def artefacts(kegg_dump):
    reactions = parse_kegg_reactions(kegg_dump)
    compounds = parse_kegg_compounds(kegg_dump)
    linked = {ko for r in reactions for ko in r.kos}
    from raven_toolbox.reconstruction.kegg import parse_kegg_kos

    kos = parse_kegg_kos(kegg_dump, keep=linked)
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
def test_bbb_keeps_only_its_reactions(artefacts):
    # bbb has GENE03 -> K90002 -> R90100 only.
    model = _build(artefacts, "bbb")
    assert {r.id for r in model.reactions} == {"R90100"}
    assert model.id == "bbb"


def test_bbb_gpr_and_gene_annotation(artefacts):
    model = _build(artefacts, "bbb")
    r = model.reactions.get_by_id("R90100")
    assert r.gene_reaction_rule == "GENE03"
    assert model.genes.get_by_id("GENE03").annotation["kegg.genes"] == "bbb:GENE03"
    assert r.notes["note"].startswith("Included by get_kegg_model_for_organism")


def test_aaa_or_joins_multiple_genes(artefacts):
    # aaa has GENE01 + GENE02, both -> K90001 -> R90010.
    model = _build(artefacts, "aaa")
    r = model.reactions.get_by_id("R90010")
    assert set(r.genes) == {model.genes.get_by_id("GENE01"), model.genes.get_by_id("GENE02")}
    assert r.gene_reaction_rule == "GENE01 or GENE02"


def test_case_insensitive_organism(artefacts):
    assert "R90010" in _build(artefacts, "AAA").reactions


def test_orphan_metabolites_pruned(artefacts):
    # bbb keeps only R90100 (C90002, C90003); other mets should go.
    model = _build(artefacts, "bbb")
    assert {m.id for m in model.metabolites} == {"C90002", "C90003"}


def test_reference_model_unmodified(artefacts):
    reference, _ = artefacts
    before = len(reference.reactions)
    _build(artefacts, "bbb")
    assert len(reference.reactions) == before  # worked on a copy
    assert len(reference.genes) == 0


# --------------------------------------------------------------------------- #
# Spontaneous handling
# --------------------------------------------------------------------------- #
def test_spontaneous_reaction_kept_without_genes(artefacts):
    # R90100 is spontaneous; for aaa it has no genes but is kept (no GPR).
    model = _build(artefacts, "aaa", keep_spontaneous=True)
    assert "R90100" in model.reactions
    assert model.reactions.get_by_id("R90100").gene_reaction_rule == ""


def test_spontaneous_dropped_when_disabled(artefacts):
    model = _build(artefacts, "aaa", keep_spontaneous=False)
    assert "R90100" not in model.reactions
    assert "R90010" in model.reactions  # the gene-backed reaction stays


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
        _build(artefacts, "qqq")


def test_domain_mode_needs_taxonomy(artefacts):
    with pytest.raises(ValueError, match="taxonomy"):
        _build(artefacts, "eukaryotes")


def test_domain_mode_keeps_all_domain_organisms(artefacts, kegg_dump):
    # Prokaryotes (aaa + bbb) -> R90010 (aaa genes) and R90100 (bbb genes).
    model = _build(artefacts, "prokaryotes", taxonomy=kegg_dump / "taxonomy")
    assert "R90010" in model.reactions
    assert "R90100" in model.reactions
    # Genes are organism-qualified in domain mode to stay distinct.
    assert {g.id for g in model.reactions.get_by_id("R90010").genes} == {
        "aaa:GENE01",
        "aaa:GENE02",
    }


def test_domain_mode_eukaryotes(artefacts, kegg_dump):
    # Eukaryotes (ccc) -> R90010 via ccc:GENE04/GENE05; bbb-only R90100 absent of
    # genes but it is spontaneous, so kept without GPR.
    model = _build(artefacts, "eukaryotes", taxonomy=kegg_dump / "taxonomy")
    assert {g.id for g in model.reactions.get_by_id("R90010").genes} == {
        "ccc:GENE04",
        "ccc:GENE05",
    }


def test_from_artefacts_roundtrip(kegg_dump, tmp_path):
    parse_kegg_dump(kegg_dump, tmp_path)
    model = get_kegg_model_for_organism_from_artefacts("bbb", tmp_path)
    assert {r.id for r in model.reactions} == {"R90100"}
    assert model.reactions.get_by_id("R90100").gene_reaction_rule == "GENE03"


def test_from_artefacts_domain_mode_auto_resolves_taxonomy(kegg_dump, tmp_path):
    # Domain mode must auto-resolve the taxonomy artefact from artefact_dir, without the
    # caller passing taxonomy=.
    import gzip
    import shutil

    parse_kegg_dump(kegg_dump, tmp_path)
    with open(kegg_dump / "taxonomy", "rb") as src, gzip.open(tmp_path / "taxonomy.gz", "wb") as out:
        shutil.copyfileobj(src, out)
    model = get_kegg_model_for_organism_from_artefacts("prokaryotes", tmp_path)
    # Same prokaryote-domain content (aaa + bbb) as the explicit-taxonomy build.
    assert "R90010" in model.reactions and "R90100" in model.reactions
    assert {g.id for g in model.reactions.get_by_id("R90010").genes} == {
        "aaa:GENE01",
        "aaa:GENE02",
    }
