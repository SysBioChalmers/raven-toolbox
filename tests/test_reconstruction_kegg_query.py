"""Tests for the KEGG HMM-query path (reconstruction/kegg/query.py, step 3b.5)."""
from pathlib import Path

import pandas as pd
import pytest

from raven_python.reconstruction.kegg import (
    assign_kos,
    build_kegg_tables,
    build_reference_model,
    get_kegg_model_from_sequences,
    parse_hmmscan_tblout,
    parse_kegg_compounds,
    parse_kegg_kos,
    parse_kegg_reactions,
)

DUMP = Path(__file__).parent / "data" / "kegg_dump"

# A minimal hmmscan --tblout excerpt: target(KO) accession query(gene) ... evalue ...
TBLOUT = """\
#                                                               --- full sequence ----
# target name        accession  query name  accession   E-value  score  bias
#------------------- ---------- ----------- ---------- --------- ------ -----
K01194               -          gene1       -          1e-120     400.0   0.0
K01194               -          gene2       -          1e-100     350.0   0.0
K00002               -          gene1       -          1e-10      40.0    0.0
"""


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #
def test_parse_tblout_skips_comments():
    hits = parse_hmmscan_tblout(TBLOUT)
    assert list(hits.columns) == ["ko", "gene", "evalue"]
    assert len(hits) == 3
    assert set(hits["ko"]) == {"K01194", "K00002"}
    assert hits.iloc[0]["evalue"] == 1e-120


def test_parse_tblout_empty():
    assert parse_hmmscan_tblout("# only a header\n").empty


# --------------------------------------------------------------------------- #
# assign_kos scoring/filters
# --------------------------------------------------------------------------- #
def test_cutoff_excludes_weak_hits():
    hits = parse_hmmscan_tblout(TBLOUT)
    # gene1->K00002 has evalue 1e-10, above the default cutoff 1e-30: dropped.
    assigned = assign_kos(hits)
    assert "K00002" not in assigned
    assert set(assigned["K01194"]) == {"gene1", "gene2"}


def test_loose_cutoff_keeps_hit():
    hits = parse_hmmscan_tblout(TBLOUT)
    assigned = assign_kos(hits, cutoff=1e-5, min_score_ratio_g=0.0, min_score_ratio_ko=0.0)
    assert assigned.get("K00002") == ["gene1"]


def test_min_score_ratio_ko_prunes_weak_member():
    # In one KO: best 1e-200, weak 1e-20. log(1e-20)/log(1e-200)=0.1 < 0.3 -> pruned.
    hits = pd.DataFrame(
        [("K1", "strong", 1e-200), ("K1", "weak", 1e-20)],
        columns=["ko", "gene", "evalue"],
    )
    assigned = assign_kos(hits, cutoff=1e-5, min_score_ratio_ko=0.3, min_score_ratio_g=0.0)
    assert assigned["K1"] == ["strong"]


def test_min_score_ratio_g_keeps_gene_in_best_ko_only():
    # gene g hits K1 strongly (1e-200) and K2 weakly (1e-20).
    # For the gene: log(1e-20)/log(1e-200)=0.1 < 0.8 -> K2 assignment dropped.
    hits = pd.DataFrame(
        [("K1", "g", 1e-200), ("K2", "g", 1e-20)],
        columns=["ko", "gene", "evalue"],
    )
    assigned = assign_kos(hits, cutoff=1e-5, min_score_ratio_ko=0.0, min_score_ratio_g=0.8)
    assert assigned == {"K1": ["g"]}


def test_zero_evalue_does_not_crash():
    hits = pd.DataFrame([("K1", "g", 0.0)], columns=["ko", "gene", "evalue"])
    assert assign_kos(hits) == {"K1": ["g"]}


def test_cutoff_ge_one_rejected():
    """cutoff >= 1 would let log(best_evalue)=0 through and ZeroDivisionError later
    (known_issues.md A6). Reject up front with a clear message."""
    hits = pd.DataFrame([("K1", "g", 0.5)], columns=["ko", "gene", "evalue"])
    with pytest.raises(ValueError, match="cutoff must be < 1"):
        assign_kos(hits, cutoff=1.0)


# --------------------------------------------------------------------------- #
# Model assembly via the HMM path (hmmscan mocked)
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def reference_and_tables():
    reactions = parse_kegg_reactions(DUMP)
    compounds = parse_kegg_compounds(DUMP)
    linked = {ko for r in reactions for ko in r.kos}
    kos = parse_kegg_kos(DUMP, keep=linked)
    return build_reference_model(reactions, compounds), build_kegg_tables(reactions, kos)


def test_get_model_from_sequences(reference_and_tables, monkeypatch):
    model_ref, tables = reference_and_tables
    # Mock the HMM search: K01194 -> myGeneA/myGeneB (-> R00010).
    monkeypatch.setattr(
        "raven_python.reconstruction.kegg.query.run_hmmscan",
        lambda *a, **k: (
            "K01194 - myGeneA - 1e-120 400 0\n"
            "K01194 - myGeneB - 1e-110 380 0\n"
        ),
    )
    model = get_kegg_model_from_sequences(
        "ignored.fasta",
        model_ref,
        tables["ko_reaction"],
        "ignored.hmm",
        rxn_flags=tables["rxn_flags"],
        model_id="myorg",
    )
    assert model.id == "myorg"
    r = model.reactions.get_by_id("R00010")
    assert set(r.gene_reaction_rule.split(" or ")) == {"myGeneA", "myGeneB"}
    assert r.notes["note"].endswith("(using HMMs)")
    # R00200/R00300 had no matched KOs and are not spontaneous -> absent.
    assert "R00200" not in model.reactions
