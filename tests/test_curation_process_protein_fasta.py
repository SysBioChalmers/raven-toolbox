"""Tests for process_protein_fasta_file (processProteinFastaFile port)."""
import pandas as pd
import pytest

from raven_toolbox.curation import process_protein_fasta_file

FASTA = ">NP_111.1 protein one\nMAAAAAAAAA\n>NP_222.1 protein two\nMBBBBBBBBB\n>NP_999.9 unmatched\nMZZZ\n"


@pytest.fixture
def faa_file(tmp_path):
    f = tmp_path / "proteins.faa"
    f.write_text(FASTA)
    return f


@pytest.fixture
def gene_table():
    return pd.DataFrame({
        "GenBank_protein": ["NP_111.1", "NP_222.1"],
        "locus_tag": ["LT1", "LT2"],
    })


def test_renames_matched_headers_keeps_unmatched(faa_file, gene_table, tmp_path):
    out = process_protein_fasta_file(faa_file, gene_table, "locus_tag", output_dir=tmp_path)
    assert out == tmp_path / "proteins_processed.faa"
    content = out.read_text()
    assert ">LT1\nMAAAAAAAAA\n" in content
    assert ">LT2\nMBBBBBBBBB\n" in content
    assert ">NP_999.9 unmatched\nMZZZ\n" in content  # unmatched kept as-is


def test_accepts_tsv_path_for_gene_table(faa_file, gene_table, tmp_path):
    tsv = tmp_path / "genes.tsv"
    gene_table.to_csv(tsv, sep="\t", index=False)
    out = process_protein_fasta_file(faa_file, tsv, "locus_tag", output_dir=tmp_path)
    assert ">LT1" in out.read_text()


def test_missing_faa_file_raises(gene_table, tmp_path):
    with pytest.raises(FileNotFoundError):
        process_protein_fasta_file(tmp_path / "nope.faa", gene_table, "locus_tag")


def test_missing_column_raises(faa_file, gene_table):
    with pytest.raises(ValueError, match="nonexistent_column"):
        process_protein_fasta_file(faa_file, gene_table, "nonexistent_column")


def test_empty_protein_id_row_skipped(faa_file, tmp_path):
    table = pd.DataFrame({
        "GenBank_protein": ["NP_111.1", ""],
        "locus_tag": ["LT1", "LT2"],
    })
    out = process_protein_fasta_file(faa_file, table, "locus_tag", output_dir=tmp_path)
    content = out.read_text()
    assert ">LT1" in content
    assert ">NP_222.1 protein two" in content  # row with empty protein id had no effect


def test_empty_header_col_value_still_renames_to_empty_string(faa_file, tmp_path):
    """processProteinFastaFile.m only checks GenBank_protein for emptiness,
    not header_col -- a row with a valid accession but an empty target
    value still renames that sequence's header, to an empty string."""
    table = pd.DataFrame({
        "GenBank_protein": ["NP_111.1"],
        "locus_tag": [""],
    })
    out = process_protein_fasta_file(faa_file, table, "locus_tag", output_dir=tmp_path)
    assert out.read_text().startswith(">\nMAAAAAAAAA\n")


def test_sequence_wrapped_at_80_columns(faa_file, gene_table, tmp_path):
    long_seq = "A" * 200
    faa_file.write_text(f">NP_111.1\n{long_seq}\n")
    out = process_protein_fasta_file(faa_file, gene_table, "locus_tag", output_dir=tmp_path)
    lines = out.read_text().splitlines()
    assert lines[0] == ">LT1"
    assert lines[1] == "A" * 80
    assert lines[2] == "A" * 80
    assert lines[3] == "A" * 40
