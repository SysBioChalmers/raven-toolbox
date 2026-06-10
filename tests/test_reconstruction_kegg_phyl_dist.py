"""Tests for the KEGG phylogenetic distance generator (RAVEN getPhylDist port)."""
import gzip

import numpy as np

from raven_python.reconstruction.kegg.taxonomy import parse_taxonomy_records, phyl_dist

# Tiny taxonomy: two prokaryotes sharing a lineage, a 3-deep mammal, a 2-deep fungus.
TAXONOMY = (
    "# Prokaryotes\n"
    "## Bacteria\n"
    "T1\tbsu\tT1\tBacillus subtilis\n"
    "T2\teco\tT2\tEscherichia coli\n"
    "# Eukaryotes\n"
    "## Animals\n"
    "### Mammals\n"
    "T3\thsa\tT3\tHomo sapiens (human)\n"
    "## Fungi\n"
    "T4\tsce\tT4\tSaccharomyces cerevisiae\n"
)


def _write(tmp_path):
    p = tmp_path / "taxonomy"
    p.write_text(TAXONOMY, encoding="utf-8")
    return p


def test_parse_records_ids_names_lineages(tmp_path):
    recs = parse_taxonomy_records(_write(tmp_path))
    assert [r[0] for r in recs] == ["bsu", "eco", "hsa", "sce"]
    # Names keep RAVEN's trailing parenthetical.
    assert [r[1] for r in recs] == [
        "Bacillus subtilis",
        "Escherichia coli",
        "Homo sapiens (human)",
        "Saccharomyces cerevisiae",
    ]
    lineage = {r[0]: r[2] for r in recs}
    assert lineage["hsa"] == ["Eukaryotes", "Animals", "Mammals"]
    assert lineage["sce"] == ["Eukaryotes", "Fungi"]


def test_phyl_dist_matches_raven_formula(tmp_path):
    pd = phyl_dist(_write(tmp_path))
    assert pd.ids == ["bsu", "eco", "hsa", "sce"]
    # Hand-computed from RAVEN's distMat = (Li-Lj) + min(Li,Lj) - k.
    expected = np.array(
        [
            [0, 0, 0, 1],
            [0, 0, 0, 1],
            [2, 2, 0, 2],
            [1, 1, 0, 0],
        ]
    )
    np.testing.assert_array_equal(pd.dist_matrix, expected)
    assert np.all(np.diag(pd.dist_matrix) == 0)  # self-distance is 0
    assert pd.dist_matrix[2, 0] != pd.dist_matrix[0, 2]  # RAVEN metric is asymmetric


def test_phyl_dist_only_in_kingdom_blocks_cross_domain(tmp_path):
    pd = phyl_dist(_write(tmp_path), only_in_kingdom=True)
    assert np.isinf(pd.dist_matrix[0, 2]) and np.isinf(pd.dist_matrix[0, 3])  # prok vs euk
    assert np.isfinite(pd.dist_matrix[2, 3])  # both eukaryotes


def test_phyl_dist_reads_gzip(tmp_path):
    p = tmp_path / "taxonomy.gz"
    p.write_bytes(gzip.compress(TAXONOMY.encode("utf-8")))
    assert phyl_dist(p).ids == ["bsu", "eco", "hsa", "sce"]
