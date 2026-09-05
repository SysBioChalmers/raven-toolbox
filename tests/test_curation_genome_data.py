"""Tests for genome_data (curation/genome_data.py, getGeneData/downloadGenomeData port)."""
import io
import zipfile
from unittest.mock import patch

import pytest

from raven_toolbox.curation import (
    download_genome_data,
    get_gene_data,
    parse_gff_gene_table,
)

EUKARYOTE_GFF = """\
##gff-version 3
chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene1;Name=rbcL;locus_tag=Cre01.g000001;old_locus_tag=CreOld1;Dbxref=GeneID:5723799,UniProtKB:P12345
chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=rna1;Parent=gene1
chr1\tsrc\tCDS\t1\t50\t.\t+\t0\tID=cds1;Parent=rna1;protein_id=XP_001698190.2
chr1\tsrc\tCDS\t51\t100\t.\t+\t0\tID=cds1;Parent=rna1;protein_id=XP_001698190.2
"""

PROKARYOTE_GFF = """\
##gff-version 3
chr1\tsrc\tgene\t200\t300\t.\t+\t.\tID=gene2;locus_tag=b0002
chr1\tsrc\tCDS\t200\t300\t.\t+\t0\tID=cds2;Parent=gene2;protein_id=WP_000000001.1
"""


@pytest.fixture
def eukaryote_gff(tmp_path):
    f = tmp_path / "euk.gff"
    f.write_text(EUKARYOTE_GFF)
    return f


@pytest.fixture
def prokaryote_gff(tmp_path):
    f = tmp_path / "prok.gff"
    f.write_text(PROKARYOTE_GFF)
    return f


def test_eukaryote_cds_mrna_gene_chain_resolved(eukaryote_gff):
    table = parse_gff_gene_table(eukaryote_gff)
    assert len(table) == 1  # two CDS (exons) dedup to one row
    row = table.iloc[0]
    assert row["locus_tag"] == "Cre01.g000001"
    assert row["old_locus_tag"] == "CreOld1"
    assert row["GeneID"] == "5723799"
    assert row["gene_name"] == "rbcL"
    assert row["GenBank_protein"] == "XP_001698190.2"
    assert row["UniProt"] == "P12345"


def test_prokaryote_cds_gene_direct_parent(prokaryote_gff):
    table = parse_gff_gene_table(prokaryote_gff)
    assert len(table) == 1
    row = table.iloc[0]
    assert row["locus_tag"] == "b0002"
    assert row["GenBank_protein"] == "WP_000000001.1"


def test_protein_id_falls_back_to_dbxref_genbank(tmp_path):
    gff = tmp_path / "fallback.gff"
    gff.write_text(
        "chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=g1;locus_tag=tag1\n"
        "chr1\tsrc\tCDS\t1\t10\t.\t+\t0\tID=c1;Parent=g1;Dbxref=GenBank:WP_999.1\n"
    )
    table = parse_gff_gene_table(gff)
    assert table.iloc[0]["GenBank_protein"] == "WP_999.1"


def test_cds_with_no_resolvable_gene_gets_empty_gene_fields(tmp_path):
    gff = tmp_path / "orphan.gff"
    gff.write_text("chr1\tsrc\tCDS\t1\t10\t.\t+\t0\tID=c1;Parent=missing_gene;protein_id=XP_1.1\n")
    table = parse_gff_gene_table(gff)
    assert len(table) == 1
    row = table.iloc[0]
    assert row["locus_tag"] == ""
    assert row["GenBank_protein"] == "XP_1.1"


def test_cds_missing_parent_or_protein_id_is_skipped(tmp_path):
    gff = tmp_path / "skip.gff"
    gff.write_text(
        "chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=g1;locus_tag=tag1\n"
        "chr1\tsrc\tCDS\t1\t10\t.\t+\t0\tID=c1;Parent=g1\n"  # no protein_id anywhere
        "chr1\tsrc\tCDS\t1\t10\t.\t+\t0\tID=c2;protein_id=XP_2.1\n"  # no Parent
    )
    table = parse_gff_gene_table(gff)
    assert table.empty


def test_percent_encoded_attribute_is_decoded(tmp_path):
    gff = tmp_path / "encoded.gff"
    gff.write_text(
        "chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=g1;Name=alpha%2Cbeta;locus_tag=tag1\n"
        "chr1\tsrc\tCDS\t1\t10\t.\t+\t0\tID=c1;Parent=g1;protein_id=XP_1.1\n"
    )
    table = parse_gff_gene_table(gff)
    assert table.iloc[0]["gene_name"] == "alpha,beta"


def test_dedup_key_excludes_uniprot(tmp_path):
    """getGeneData.m's own dedup key omits UniProt; two rows differing only
    in UniProt collapse to one, keeping the first's UniProt value."""
    gff = tmp_path / "dupuniprot.gff"
    gff.write_text(
        "chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=g1;locus_tag=tag1;Dbxref=UniProtKB:P1\n"
        "chr1\tsrc\tCDS\t1\t5\t.\t+\t0\tID=c1;Parent=g1;protein_id=XP_1.1\n"
    )
    table = parse_gff_gene_table(gff)
    assert len(table) == 1
    assert table.iloc[0]["UniProt"] == "P1"


def test_comment_and_blank_lines_ignored(tmp_path):
    gff = tmp_path / "comments.gff"
    gff.write_text(
        "##gff-version 3\n\n# a comment\n"
        "chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=g1;locus_tag=tag1\n"
        "chr1\tsrc\tCDS\t1\t10\t.\t+\t0\tID=c1;Parent=g1;protein_id=XP_1.1\n"
    )
    table = parse_gff_gene_table(gff)
    assert len(table) == 1


# --- get_gene_data -----------------------------------------------------

def test_get_gene_data_uses_local_file_directly(eukaryote_gff):
    table = get_gene_data(eukaryote_gff)
    assert len(table) == 1


def test_get_gene_data_writes_output_file(eukaryote_gff, tmp_path):
    out = tmp_path / "out.tsv"
    get_gene_data(eukaryote_gff, output_file=out)
    assert out.is_file()
    assert "locus_tag" in out.read_text()


def test_get_gene_data_rejects_bad_input():
    with pytest.raises(ValueError, match="neither an existing file"):
        get_gene_data("not_a_real_path_or_accession")


def test_get_gene_data_downloads_when_given_accession(eukaryote_gff, tmp_path):
    with patch("raven_toolbox.curation.genome_data.download_genome_data") as mock_dl:
        mock_dl.return_value = (eukaryote_gff, tmp_path / "unused.faa")
        table = get_gene_data("GCF_000002595.2", download_dir=tmp_path)
    mock_dl.assert_called_once_with("GCF_000002595.2", output_dir=tmp_path)
    assert len(table) == 1


# --- download_genome_data ------------------------------------------------

def test_download_rejects_bad_accession_prefix():
    with pytest.raises(ValueError, match="GCF_.*GCA_"):
        download_genome_data("XYZ_12345.1")


def test_download_skips_if_files_already_present(tmp_path):
    gff = tmp_path / "GCF_1_genomic.gff"
    faa = tmp_path / "GCF_1_protein.faa"
    gff.write_text("x")
    faa.write_text("y")
    with patch("raven_toolbox.curation.genome_data.requests.get") as mock_get:
        result = download_genome_data("GCF_1", output_dir=tmp_path, verbose=False)
    mock_get.assert_not_called()
    assert result == (gff, faa)


def _zip_bytes(files: dict[str, str]) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for name, content in files.items():
            zf.writestr(name, content)
    return buf.getvalue()


def test_download_extracts_and_moves_files(tmp_path):
    archive = _zip_bytes({
        "ncbi_dataset/data/GCF_1.2/genomic.gff": "gff-content",
        "ncbi_dataset/data/GCF_1.2/protein.faa": "faa-content",
    })

    class FakeResponse:
        content = archive
        def raise_for_status(self): pass

    with patch("raven_toolbox.curation.genome_data.requests.get", return_value=FakeResponse()):
        gff_file, faa_file = download_genome_data("GCF_1", output_dir=tmp_path, verbose=False)

    assert gff_file.read_text() == "gff-content"
    assert faa_file.read_text() == "faa-content"


def test_download_raises_when_archive_missing_annotation(tmp_path):
    archive = _zip_bytes({"ncbi_dataset/data/GCF_1.2/genomic.gff": "gff-only"})

    class FakeResponse:
        content = archive
        def raise_for_status(self): pass

    with patch("raven_toolbox.curation.genome_data.requests.get", return_value=FakeResponse()):
        with pytest.raises(RuntimeError, match="did not return both"):
            download_genome_data("GCF_1", output_dir=tmp_path, verbose=False)
