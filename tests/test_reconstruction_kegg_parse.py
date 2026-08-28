"""Tests for the KEGG dump parser (reconstruction/kegg/parse.py, step 3b.2).

The ``kegg_dump`` fixture (tests/conftest.py) is a small, fully fictional dump —
no real KEGG content is committed.
"""
import pytest

from raven_toolbox.reconstruction.kegg import (
    build_kegg_tables,
    build_reference_model,
    parse_kegg_compounds,
    parse_kegg_dump,
    parse_kegg_kos,
    parse_kegg_reactions,
    read_kegg_table,
    write_kegg_tables,
)


@pytest.fixture(scope="module")
def reactions(kegg_dump):
    return parse_kegg_reactions(kegg_dump)


@pytest.fixture(scope="module")
def compounds(kegg_dump):
    return parse_kegg_compounds(kegg_dump)


@pytest.fixture(scope="module")
def kos(kegg_dump):
    linked = {ko for r in parse_kegg_reactions(kegg_dump) for ko in r.kos}
    return parse_kegg_kos(kegg_dump, keep=linked)


# --------------------------------------------------------------------------- #
# Reactions
# --------------------------------------------------------------------------- #
def test_reactions_parsed(reactions):
    assert {r.id for r in reactions} == {"R90010", "R90100", "R90200", "R90300", "R90400"}


def test_reaction_fields(reactions):
    r = next(r for r in reactions if r.id == "R90010")
    assert r.name == "fictional glucohydrolase analogue"
    assert r.eccodes == ["9.9.9.99"]
    assert r.kos == ["K90001"]
    # rn01199 is an overview map (rn011..) and must be skipped.
    assert r.pathways == ["rn09500"]


def test_stoichiometry_cached(reactions):
    """parse_kegg_reactions populates the cached stoichiometry so
    build_reference_model doesn't have to re-parse."""
    r = next(r for r in reactions if r.id == "R90010")
    assert r.stoichiometry  # non-empty
    # Reactants negative, products positive.
    assert all(c != 0 for c in r.stoichiometry.values())
    assert any(c < 0 for c in r.stoichiometry.values())
    assert any(c > 0 for c in r.stoichiometry.values())


def test_spontaneous_flag(reactions):
    assert next(r for r in reactions if r.id == "R90100").spontaneous
    assert not next(r for r in reactions if r.id == "R90010").spontaneous


def test_general_flag(reactions):
    assert next(r for r in reactions if r.id == "R90300").general


def test_undefined_stoich_flag(reactions):
    assert next(r for r in reactions if r.id == "R90200").undefined_stoich
    assert not next(r for r in reactions if r.id == "R90010").undefined_stoich


def test_mapformula_makes_irreversible(reactions):
    # R90100 is drawn one direction in its only map -> irreversible.
    assert not next(r for r in reactions if r.id == "R90100").reversible
    # R90010 is drawn in conflicting directions across maps -> stays reversible.
    assert next(r for r in reactions if r.id == "R90010").reversible


# --------------------------------------------------------------------------- #
# Compounds
# --------------------------------------------------------------------------- #
def test_compound_first_name_only(compounds):
    water = next(c for c in compounds if c.id == "C90001")
    assert water.name == "Aqualike"
    assert water.chebi == ["CHEBI:95377"]
    assert water.pubchem == ["9303"]


def test_inchi_overrides_formula(compounds):
    glucose = next(c for c in compounds if c.id == "C90031")
    assert glucose.inchi.startswith("InChI=")
    assert glucose.formula == ""  # cleared when an InChI is available
    assert glucose.chebi == ["CHEBI:94167", "CHEBI:97634"]


# --------------------------------------------------------------------------- #
# KOs / genes
# --------------------------------------------------------------------------- #
def test_kos_limited_to_keep(kos):
    # K90099 is unlinked (excluded by keep); K90009 is referenced but absent.
    assert {ko.id for ko in kos} == {"K90001", "K90002"}


def test_ko_genes_lowercased_and_stripped(kos):
    k = next(ko for ko in kos if ko.id == "K90001")
    assert k.name == "fictional trehalase analogue [EC:9.9.9.99]"
    assert ("aaa", "GENE01") in k.genes  # '(alias1)' suffix stripped, org lowercased
    assert ("ccc", "GENE05") in k.genes


# --------------------------------------------------------------------------- #
# Reference model
# --------------------------------------------------------------------------- #
def test_reference_model_is_gene_free(reactions, compounds):
    model = build_reference_model(reactions, compounds)
    assert len(model.genes) == 0
    for rxn in model.reactions:
        assert rxn.gene_reaction_rule == ""


def test_empty_reaction_dropped(reactions, compounds):
    model = build_reference_model(reactions, compounds)
    assert "R90400" not in model.reactions  # C90007 <=> C90007 cancels out
    assert "C90007" not in model.metabolites  # and its only metabolite is unused


def test_reaction_bounds_follow_reversibility(reactions, compounds):
    model = build_reference_model(reactions, compounds)
    assert model.reactions.get_by_id("R90010").bounds == (-1000.0, 1000.0)
    assert model.reactions.get_by_id("R90100").bounds == (0.0, 1000.0)


def test_reaction_stoichiometry_and_annotation(reactions, compounds):
    model = build_reference_model(reactions, compounds)
    r = model.reactions.get_by_id("R90010")
    coefs = {m.id: c for m, c in r.metabolites.items()}
    assert coefs == {"C91083": -1.0, "C90001": -1.0, "C90031": 2.0}
    assert r.annotation["kegg.orthology"] == ["K90001"]
    assert r.annotation["ec-code"] == ["9.9.9.99"]


def test_metabolite_annotation(reactions, compounds):
    model = build_reference_model(reactions, compounds)
    glucose = model.metabolites.get_by_id("C90031")
    assert glucose.name == "Glucolike"
    assert glucose.annotation["inchi"].startswith("InChI=")


# --------------------------------------------------------------------------- #
# Tables
# --------------------------------------------------------------------------- #
def test_ko_reaction_table(reactions, kos):
    tables = build_kegg_tables(reactions, kos)
    pairs = set(map(tuple, tables["ko_reaction"].to_numpy()))
    assert ("K90001", "R90010") in pairs
    assert ("K90009", "R90300") in pairs  # kept even though KO entry is missing


def test_organism_gene_ko_table(reactions, kos):
    tables = build_kegg_tables(reactions, kos)
    rows = set(map(tuple, tables["organism_gene_ko"].to_numpy()))
    assert ("aaa", "GENE01", "K90001") in rows
    assert ("bbb", "GENE03", "K90002") in rows
    assert len(rows) == 5


def test_rxn_flags_table(reactions, kos):
    tables = build_kegg_tables(reactions, kos)
    flags = tables["rxn_flags"].set_index("reaction")
    assert bool(flags.loc["R90100", "spontaneous"])
    assert bool(flags.loc["R90200", "undefined_stoich"])
    assert bool(flags.loc["R90300", "general"])
    assert not bool(flags.loc["R90010", "spontaneous"])


# --------------------------------------------------------------------------- #
# Round-trip + orchestrator
# --------------------------------------------------------------------------- #
def test_tables_roundtrip_gzipped_tsv(reactions, kos, tmp_path):
    tables = build_kegg_tables(reactions, kos)
    paths = write_kegg_tables(tables, tmp_path)
    assert all(p.name.endswith(".tsv.gz") for p in paths)
    back = read_kegg_table(tmp_path / "ko_reaction.tsv.gz")
    assert set(map(tuple, back.to_numpy())) == set(map(tuple, tables["ko_reaction"].to_numpy()))


def test_parse_kegg_dump_writes_artefacts(kegg_dump, tmp_path):
    paths = parse_kegg_dump(kegg_dump, tmp_path)
    assert set(paths) >= {
        "ko_reaction", "ko_names", "organism_gene_ko", "rxn_flags", "reference_model"
    }
    assert (tmp_path / "reference_model.yml.gz").is_file()
    # organism_gene_ko is streamed to a sorted, gzipped TSV.
    assert paths["organism_gene_ko"].name == "organism_gene_ko.tsv.gz"
    ogk = read_kegg_table(paths["organism_gene_ko"])
    assert set(ogk.columns) == {"organism", "gene", "ko"}
    assert ("bbb", "GENE03", "K90002") in set(map(tuple, ogk.to_numpy()))
    # Rows are sorted by (organism, gene) — the property that makes them compress.
    keys = list(zip(ogk["organism"], ogk["gene"], strict=True))
    assert keys == sorted(keys)


def test_parse_kegg_dump_version_prefixes_filenames(kegg_dump, tmp_path):
    paths = parse_kegg_dump(kegg_dump, tmp_path, version="kegg116")
    # Dict keys stay logical; the files on disk are version-prefixed.
    assert set(paths) >= {"ko_reaction", "organism_gene_ko", "reference_model"}
    assert paths["organism_gene_ko"].name == "kegg116_organism_gene_ko.tsv.gz"
    assert paths["reference_model"].name == "kegg116_reference_model.yml.gz"
    assert (tmp_path / "kegg116_ko_reaction.tsv.gz").is_file()
    assert read_kegg_table(paths["organism_gene_ko"]).columns.tolist() == ["organism", "gene", "ko"]


def test_stream_organism_gene_ko_external_merge(kegg_dump, tmp_path):
    """A tiny chunk_rows forces multiple sorted runs to be merged; output stays sorted."""
    from raven_toolbox.reconstruction.kegg.parse import stream_organism_gene_ko

    out = tmp_path / "organism_gene_ko.tsv.gz"
    keep = {ko.id for ko in parse_kegg_kos(kegg_dump)}
    names = stream_organism_gene_ko(kegg_dump, keep, out, chunk_rows=1)
    assert out.is_file() and not list(tmp_path.glob("ogk_sort_*"))  # temp dir cleaned up
    ogk = read_kegg_table(out)
    keys = list(zip(ogk["organism"], ogk["gene"], strict=True))
    assert keys == sorted(keys)
    assert ("bbb", "GENE03", "K90002") in set(map(tuple, ogk.to_numpy()))
    assert set(names.columns) == {"ko", "name"}
