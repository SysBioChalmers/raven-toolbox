"""Tests for the KEGG dump parser (reconstruction/kegg/parse.py, step 3b.2)."""
from pathlib import Path

import pytest

from raven_python.reconstruction.kegg import (
    build_kegg_tables,
    build_reference_model,
    parse_kegg_compounds,
    parse_kegg_dump,
    parse_kegg_kos,
    parse_kegg_reactions,
    read_kegg_table,
    write_kegg_tables,
)

DUMP = Path(__file__).parent / "data" / "kegg_dump"


@pytest.fixture(scope="module")
def reactions():
    return parse_kegg_reactions(DUMP)


@pytest.fixture(scope="module")
def compounds():
    return parse_kegg_compounds(DUMP)


@pytest.fixture(scope="module")
def kos():
    linked = {ko for r in parse_kegg_reactions(DUMP) for ko in r.kos}
    return parse_kegg_kos(DUMP, keep=linked)


# --------------------------------------------------------------------------- #
# Reactions
# --------------------------------------------------------------------------- #
def test_reactions_parsed(reactions):
    assert {r.id for r in reactions} == {"R00010", "R00100", "R00200", "R00300", "R00400"}


def test_reaction_fields(reactions):
    r = next(r for r in reactions if r.id == "R00010")
    assert r.name == "alpha,alpha-trehalose glucohydrolase"
    assert r.eccodes == ["3.2.1.28"]
    assert r.kos == ["K01194"]
    # rn01100 is an overview map and must be skipped.
    assert r.pathways == ["rn00500"]


def test_stoichiometry_cached(reactions):
    """parse_kegg_reactions populates the cached stoichiometry so
    build_reference_model doesn't have to re-parse (known_issues.md D2)."""
    r = next(r for r in reactions if r.id == "R00010")
    assert r.stoichiometry  # non-empty
    # Reactants negative, products positive.
    assert all(c != 0 for c in r.stoichiometry.values())
    assert any(c < 0 for c in r.stoichiometry.values())
    assert any(c > 0 for c in r.stoichiometry.values())


def test_spontaneous_flag(reactions):
    assert next(r for r in reactions if r.id == "R00100").spontaneous
    assert not next(r for r in reactions if r.id == "R00010").spontaneous


def test_general_flag(reactions):
    assert next(r for r in reactions if r.id == "R00300").general


def test_undefined_stoich_flag(reactions):
    assert next(r for r in reactions if r.id == "R00200").undefined_stoich
    assert not next(r for r in reactions if r.id == "R00010").undefined_stoich


def test_mapformula_makes_irreversible(reactions):
    # R00100 is drawn one direction in its only map -> irreversible.
    assert not next(r for r in reactions if r.id == "R00100").reversible
    # R00010 is drawn in conflicting directions across maps -> stays reversible.
    assert next(r for r in reactions if r.id == "R00010").reversible


# --------------------------------------------------------------------------- #
# Compounds
# --------------------------------------------------------------------------- #
def test_compound_first_name_only(compounds):
    water = next(c for c in compounds if c.id == "C00001")
    assert water.name == "H2O"
    assert water.chebi == ["CHEBI:15377"]
    assert water.pubchem == ["3303"]


def test_inchi_overrides_formula(compounds):
    glucose = next(c for c in compounds if c.id == "C00031")
    assert glucose.inchi.startswith("InChI=")
    assert glucose.formula == ""  # cleared when an InChI is available
    assert glucose.chebi == ["CHEBI:4167", "CHEBI:17634"]


# --------------------------------------------------------------------------- #
# KOs / genes
# --------------------------------------------------------------------------- #
def test_kos_limited_to_keep(kos):
    # K99999 is unlinked (excluded by keep); K09999 is referenced but absent.
    assert {ko.id for ko in kos} == {"K01194", "K00002"}


def test_ko_genes_lowercased_and_stripped(kos):
    k = next(ko for ko in kos if ko.id == "K01194")
    assert k.name == "alpha,alpha-trehalase [EC:3.2.1.28]"
    assert ("bsu", "BSU31050") in k.genes  # '(gbsB)' suffix stripped, org lowercased
    assert ("hsa", "125") in k.genes


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
    assert "R00400" not in model.reactions  # C00007 <=> C00007 cancels out
    assert "C00007" not in model.metabolites  # and its only metabolite is unused


def test_reaction_bounds_follow_reversibility(reactions, compounds):
    model = build_reference_model(reactions, compounds)
    assert model.reactions.get_by_id("R00010").bounds == (-1000.0, 1000.0)
    assert model.reactions.get_by_id("R00100").bounds == (0.0, 1000.0)


def test_reaction_stoichiometry_and_annotation(reactions, compounds):
    model = build_reference_model(reactions, compounds)
    r = model.reactions.get_by_id("R00010")
    coefs = {m.id: c for m, c in r.metabolites.items()}
    assert coefs == {"C01083": -1.0, "C00001": -1.0, "C00031": 2.0}
    assert r.annotation["kegg.orthology"] == ["K01194"]
    assert r.annotation["ec-code"] == ["3.2.1.28"]


def test_metabolite_annotation(reactions, compounds):
    model = build_reference_model(reactions, compounds)
    glucose = model.metabolites.get_by_id("C00031")
    assert glucose.name == "D-Glucose"
    assert glucose.annotation["inchi"].startswith("InChI=")


# --------------------------------------------------------------------------- #
# Tables
# --------------------------------------------------------------------------- #
def test_ko_reaction_table(reactions, kos):
    tables = build_kegg_tables(reactions, kos)
    pairs = set(map(tuple, tables["ko_reaction"].to_numpy()))
    assert ("K01194", "R00010") in pairs
    assert ("K09999", "R00300") in pairs  # kept even though KO entry is missing


def test_organism_gene_ko_table(reactions, kos):
    tables = build_kegg_tables(reactions, kos)
    rows = set(map(tuple, tables["organism_gene_ko"].to_numpy()))
    assert ("bsu", "BSU31050", "K01194") in rows
    assert ("eco", "b0001", "K00002") in rows
    assert len(rows) == 5


def test_rxn_flags_table(reactions, kos):
    tables = build_kegg_tables(reactions, kos)
    flags = tables["rxn_flags"].set_index("reaction")
    assert bool(flags.loc["R00100", "spontaneous"])
    assert bool(flags.loc["R00200", "undefined_stoich"])
    assert bool(flags.loc["R00300", "general"])
    assert not bool(flags.loc["R00010", "spontaneous"])


# --------------------------------------------------------------------------- #
# Round-trip + orchestrator
# --------------------------------------------------------------------------- #
def test_tables_roundtrip_gzipped_tsv(reactions, kos, tmp_path):
    tables = build_kegg_tables(reactions, kos)
    paths = write_kegg_tables(tables, tmp_path)
    assert all(p.name.endswith(".tsv.gz") for p in paths)
    back = read_kegg_table(tmp_path / "ko_reaction.tsv.gz")
    assert set(map(tuple, back.to_numpy())) == set(map(tuple, tables["ko_reaction"].to_numpy()))


def test_parse_kegg_dump_writes_artefacts(tmp_path):
    paths = parse_kegg_dump(DUMP, tmp_path)
    assert set(paths) >= {
        "ko_reaction", "ko_names", "organism_gene_ko", "rxn_flags", "reference_model"
    }
    assert (tmp_path / "reference_model.yml.gz").is_file()
    # organism_gene_ko is streamed to a sorted, gzipped TSV.
    assert paths["organism_gene_ko"].name == "organism_gene_ko.tsv.gz"
    ogk = read_kegg_table(paths["organism_gene_ko"])
    assert set(ogk.columns) == {"organism", "gene", "ko"}
    assert ("eco", "b0001", "K00002") in set(map(tuple, ogk.to_numpy()))
    # Rows are sorted by (organism, gene) — the property that makes them compress.
    keys = list(zip(ogk["organism"], ogk["gene"], strict=True))
    assert keys == sorted(keys)


def test_parse_kegg_dump_version_prefixes_filenames(tmp_path):
    paths = parse_kegg_dump(DUMP, tmp_path, version="kegg116")
    # Dict keys stay logical; the files on disk are version-prefixed.
    assert set(paths) >= {"ko_reaction", "organism_gene_ko", "reference_model"}
    assert paths["organism_gene_ko"].name == "kegg116_organism_gene_ko.tsv.gz"
    assert paths["reference_model"].name == "kegg116_reference_model.yml.gz"
    assert (tmp_path / "kegg116_ko_reaction.tsv.gz").is_file()
    assert read_kegg_table(paths["organism_gene_ko"]).columns.tolist() == ["organism", "gene", "ko"]


def test_stream_organism_gene_ko_external_merge(tmp_path):
    """A tiny chunk_rows forces multiple sorted runs to be merged; output stays sorted."""
    from raven_python.reconstruction.kegg.parse import stream_organism_gene_ko

    out = tmp_path / "organism_gene_ko.tsv.gz"
    keep = {ko.id for ko in parse_kegg_kos(DUMP)}
    names = stream_organism_gene_ko(DUMP, keep, out, chunk_rows=1)
    assert out.is_file() and not list(tmp_path.glob("ogk_sort_*"))  # temp dir cleaned up
    ogk = read_kegg_table(out)
    keys = list(zip(ogk["organism"], ogk["gene"], strict=True))
    assert keys == sorted(keys)
    assert ("eco", "b0001", "K00002") in set(map(tuple, ogk.to_numpy()))
    assert set(names.columns) == {"ko", "name"}
