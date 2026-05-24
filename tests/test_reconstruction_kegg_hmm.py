"""Tests for KEGG HMM-library construction (taxonomy + hmm, step 3b.3)."""
from pathlib import Path

import pandas as pd
import pytest

from ravengem.reconstruction.kegg import (
    build_ko_fastas,
    organism_domains,
    organisms_in_domain,
    parse_taxonomy,
)
from ravengem.reconstruction.kegg.hmm import (
    _cdhit_cmd,
    _cdhit_word_size,
    _hmmbuild_cmd,
    _mafft_cmd,
    _subsample_fasta,
    build_ko_hmm,
)

DUMP = Path(__file__).parent / "data" / "kegg_dump"


@pytest.fixture
def organism_gene_ko():
    return pd.DataFrame(
        [
            ("bsu", "BSU31050", "K01194"),
            ("bsu", "BSU31060", "K01194"),
            ("hsa", "124", "K01194"),
            ("hsa", "125", "K01194"),
            ("eco", "b0001", "K00002"),
        ],
        columns=["organism", "gene", "ko"],
    )


# --------------------------------------------------------------------------- #
# Taxonomy
# --------------------------------------------------------------------------- #
def test_parse_taxonomy_lineages():
    cats = parse_taxonomy(DUMP / "taxonomy")
    assert cats["bsu"] == ["Prokaryotes", "Bacteria", "Firmicutes"]
    assert cats["hsa"][0] == "Eukaryotes"
    assert cats["eco"][1] == "Bacteria"


def test_organism_domains():
    assert organism_domains(DUMP / "taxonomy") == {
        "bsu": "Prokaryotes",
        "eco": "Prokaryotes",
        "hsa": "Eukaryotes",
    }


def test_organisms_in_domain_prefix_match():
    assert organisms_in_domain(DUMP / "taxonomy", "prok") == {"bsu", "eco"}
    assert organisms_in_domain(DUMP / "taxonomy", "Eukaryotes") == {"hsa"}


# --------------------------------------------------------------------------- #
# build_ko_fastas (constructMultiFasta)
# --------------------------------------------------------------------------- #
def test_build_ko_fastas_groups_by_ko(organism_gene_ko, tmp_path):
    written = build_ko_fastas(organism_gene_ko, DUMP / "genes.pep", tmp_path)
    assert set(written) == {"K01194", "K00002"}
    k01194 = (tmp_path / "K01194.fa").read_text()
    assert k01194.count(">") == 4  # bsu x2 + hsa x2
    assert ">bsu:BSU31050" in k01194
    assert ">xxx:unused" not in k01194  # gene not in any KO is excluded


def test_build_ko_fastas_domain_filter(organism_gene_ko, tmp_path):
    prok = organisms_in_domain(DUMP / "taxonomy", "prokaryotes")
    written = build_ko_fastas(organism_gene_ko, DUMP / "genes.pep", tmp_path, organisms=prok)
    # Only prokaryote genes: K01194 keeps bsu (2), K00002 keeps eco (1).
    assert (tmp_path / "K01194.fa").read_text().count(">") == 2
    assert ">hsa:" not in (tmp_path / "K01194.fa").read_text()
    assert set(written) == {"K01194", "K00002"}


def test_build_ko_fastas_sequences_intact(organism_gene_ko, tmp_path):
    build_ko_fastas(organism_gene_ko, DUMP / "genes.pep", tmp_path)
    text = (tmp_path / "K00002.fa").read_text()
    assert text.startswith(">eco:b0001")
    assert "MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDA" in text


# --------------------------------------------------------------------------- #
# Command builders / CD-HIT word size (pure)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "identity, expected",
    [(0.9, "5"), (0.7, "4"), (0.65, "4"), (0.55, "3"), (0.45, "2")],
)
def test_cdhit_word_size(identity, expected):
    assert _cdhit_word_size(identity) == expected


def test_cdhit_word_size_out_of_range():
    with pytest.raises(ValueError, match="seq_identity"):
        _cdhit_word_size(0.3)


def test_command_builders():
    cd = _cdhit_cmd("cd-hit", Path("in.fa"), Path("out.fa"), 0.9, 4)
    assert cd[:3] == ["cd-hit", "-i", "in.fa"]
    assert "-c" in cd and "0.9" in cd and "-n" in cd and "5" in cd
    # Default is fast progressive (FFT-NS-2), not --auto.
    assert _mafft_cmd("mafft", Path("in.fa"), 2) == [
        "mafft", "--retree", "2", "--maxiterate", "0", "--anysymbol", "--thread", "2", "in.fa"
    ]
    assert _mafft_cmd("mafft", Path("in.fa"), 2, fast=False)[:2] == ["mafft", "--auto"]
    assert "--parttree" in _mafft_cmd("mafft", Path("in.fa"), 2, parttree=True)
    assert _hmmbuild_cmd("hmmbuild", Path("o.hmm"), Path("a.fa"), 3) == [
        "hmmbuild", "--cpu", "3", "o.hmm", "a.fa"
    ]


# --------------------------------------------------------------------------- #
# build_ko_hmm orchestration (binaries mocked)
# --------------------------------------------------------------------------- #
def test_build_ko_hmm_multi_sequence_runs_full_pipeline(tmp_path, monkeypatch):
    fasta = tmp_path / "K01194.fa"
    fasta.write_text(">a\nMKV\n>b\nMRV\n")
    calls = []

    monkeypatch.setattr(
        "ravengem.reconstruction.kegg.hmm.resolve_binary",
        lambda exe, binary=None: binary or exe,
    )

    def fake_run(cmd, *, stdout_path=None):
        calls.append(Path(cmd[0]).name)
        # Emulate each tool producing its expected output file.
        if stdout_path is not None:
            Path(stdout_path).write_text(">a\nMKV\n>b\nMRV\n")
        if Path(cmd[0]).name == "cd-hit":
            Path(cmd[cmd.index("-o") + 1]).write_text(">a\nMKV\n>b\nMRV\n")
        if Path(cmd[0]).name == "hmmbuild":
            Path(cmd[-2]).write_text("HMM\n")
        return ""

    monkeypatch.setattr("ravengem.reconstruction.kegg.hmm._run", fake_run)
    out = build_ko_hmm(fasta, tmp_path / "K01194.hmm")
    assert calls == ["cd-hit", "mafft", "hmmbuild"]
    assert out.read_text() == "HMM\n"


def test_build_ko_hmm_single_sequence_skips_align(tmp_path, monkeypatch):
    fasta = tmp_path / "K9.fa"
    fasta.write_text(">only\nMKV\n")
    calls = []
    monkeypatch.setattr(
        "ravengem.reconstruction.kegg.hmm.resolve_binary",
        lambda exe, binary=None: binary or exe,
    )

    def fake_run(cmd, *, stdout_path=None):
        calls.append(Path(cmd[0]).name)
        if Path(cmd[0]).name == "hmmbuild":
            Path(cmd[-2]).write_text("HMM\n")
        return ""

    monkeypatch.setattr("ravengem.reconstruction.kegg.hmm._run", fake_run)
    build_ko_hmm(fasta, tmp_path / "K9.hmm")
    assert calls == ["hmmbuild"]  # no cd-hit / mafft for a lone sequence


def test_subsample_fasta_caps_and_spreads(tmp_path):
    src = tmp_path / "in.fa"
    src.write_text("".join(f">s{i}\nSEQ{i}\n" for i in range(10)))
    dst = tmp_path / "out.fa"
    assert _subsample_fasta(src, dst, 4) == 4
    kept = [ln[1:].strip() for ln in dst.read_text().splitlines() if ln.startswith(">")]
    assert kept == ["s0", "s2", "s5", "s7"]  # evenly spaced across the 10


def test_subsample_fasta_noop_when_under_cap(tmp_path):
    src = tmp_path / "in.fa"
    src.write_text(">a\nMK\n>b\nMR\n")
    dst = tmp_path / "out.fa"
    assert _subsample_fasta(src, dst, 5) == 2
    assert dst.read_text() == src.read_text()


def test_build_ko_hmm_caps_before_mafft(tmp_path, monkeypatch):
    fasta = tmp_path / "K.fa"
    fasta.write_text("".join(f">g{i}\nMKV\n" for i in range(20)))
    monkeypatch.setattr(
        "ravengem.reconstruction.kegg.hmm.resolve_binary", lambda exe, binary=None: binary or exe
    )
    seen = {}

    def fake_run(cmd, *, stdout_path=None):
        name = Path(cmd[0]).name
        if name == "cd-hit":
            # CD-HIT "keeps" all 20 here.
            Path(cmd[cmd.index("-o") + 1]).write_text(fasta.read_text())
        if name == "mafft":
            seen["mafft_input_seqs"] = Path(cmd[-1]).read_text().count(">")
            Path(stdout_path).write_text(fasta.read_text())
        if name == "hmmbuild":
            Path(cmd[-2]).write_text("HMM\n")
        return ""

    monkeypatch.setattr("ravengem.reconstruction.kegg.hmm._run", fake_run)
    build_ko_hmm(fasta, tmp_path / "K.hmm", max_sequences=5)
    assert seen["mafft_input_seqs"] == 5  # capped from 20 -> 5 before MAFFT


def test_build_ko_hmm_verbose_logs_each_stage(tmp_path, monkeypatch, caplog):
    fasta = tmp_path / "K01194.fa"
    fasta.write_text(">a\nMKV\n>b\nMRV\n")
    monkeypatch.setattr(
        "ravengem.reconstruction.kegg.hmm.resolve_binary", lambda exe, binary=None: binary or exe
    )

    def fake_run(cmd, *, stdout_path=None):
        if stdout_path is not None:
            Path(stdout_path).write_text(">a\nMKV\n>b\nMRV\n")
        if Path(cmd[0]).name == "cd-hit":
            Path(cmd[cmd.index("-o") + 1]).write_text(">a\nMKV\n>b\nMRV\n")
        if Path(cmd[0]).name == "hmmbuild":
            Path(cmd[-2]).write_text("HMM\n")
        return ""

    monkeypatch.setattr("ravengem.reconstruction.kegg.hmm._run", fake_run)
    with caplog.at_level("INFO", logger="ravengem.reconstruction.kegg.hmm"):
        build_ko_hmm(fasta, tmp_path / "K01194.hmm", verbose=True)
    text = caplog.text
    # Each stage is logged, labelled with the KO id.
    assert "[K01194] start: 2 sequences" in text
    assert "[K01194] CD-HIT" in text
    assert "[K01194] MAFFT" in text
    assert "running hmmbuild" in text
    assert "[K01194] complete" in text


def test_build_ko_hmm_quiet_by_default(tmp_path, monkeypatch, caplog):
    fasta = tmp_path / "K9.fa"
    fasta.write_text(">only\nMKV\n")
    monkeypatch.setattr(
        "ravengem.reconstruction.kegg.hmm.resolve_binary", lambda exe, binary=None: binary or exe
    )
    monkeypatch.setattr(
        "ravengem.reconstruction.kegg.hmm._run",
        lambda cmd, *, stdout_path=None: Path(cmd[-2]).write_text("HMM\n") and "",
    )
    with caplog.at_level("INFO", logger="ravengem.reconstruction.kegg.hmm"):
        build_ko_hmm(fasta, tmp_path / "K9.hmm")  # verbose defaults False
    assert caplog.text == ""


def test_build_ko_hmm_empty_fasta_raises(tmp_path):
    fasta = tmp_path / "empty.fa"
    fasta.write_text("")
    with pytest.raises(ValueError, match="no sequences"):
        build_ko_hmm(fasta, tmp_path / "empty.hmm")
