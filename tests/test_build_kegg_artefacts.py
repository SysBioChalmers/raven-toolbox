"""Tests for scripts/build_kegg_artefacts.py — the resumable/idempotent build flow.

Covers the parse + core-bundle path (no ``--hmms``, so no HMMER/MAFFT/CD-HIT
needed); the per-KO HMM resume is unit-tested in test_reconstruction_kegg_hmm.py.
The ``kegg_dump`` fixture (tests/conftest.py) is a small, fully fictional dump.
"""
import importlib.util
from pathlib import Path

# scripts/ is not a package; load the module directly by path (as test_scripts_registry does).
_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_kegg_artefacts.py"
_spec = importlib.util.spec_from_file_location("build_kegg_artefacts", _SCRIPT)
bka = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bka)


def _argv(kegg_dump, out):
    return ["--keggdb", str(kegg_dump), "--out", str(out), "--version", "kegg999"]


def test_first_run_parses_and_bundles(kegg_dump, tmp_path, capsys):
    out = tmp_path / "art"
    bka.main(_argv(kegg_dump, out))
    cap = capsys.readouterr().out
    assert "Parsing KEGG dump" in cap
    assert (out / "kegg999_core.tar.gz").is_file()
    assert (out / "kegg999_taxonomy.gz").is_file()
    # Loose core members are removed once bundled.
    assert not (out / "kegg999_organism_gene_ko.tsv.gz").exists()


def test_rerun_skips_completed_stages(kegg_dump, tmp_path, capsys):
    out = tmp_path / "art"
    argv = _argv(kegg_dump, out)
    bka.main(argv)
    capsys.readouterr()

    bka.main(argv)  # nothing left to do
    cap = capsys.readouterr().out
    assert "Parsing KEGG dump" not in cap
    assert "skipping parse" in cap
    assert "core bundle" in cap and "exists; skipped" in cap


def test_force_rebuilds_despite_existing_outputs(kegg_dump, tmp_path, capsys):
    out = tmp_path / "art"
    argv = _argv(kegg_dump, out)
    bka.main(argv)
    capsys.readouterr()

    bka.main([*argv, "--force"])
    cap = capsys.readouterr().out
    assert "Parsing KEGG dump" in cap
    assert (out / "kegg999_core.tar.gz").is_file()


def test_core_paths_keys_match_parse_output(tmp_path):
    paths = bka._core_paths(tmp_path, "kegg999_")
    assert set(paths) == set(bka._CORE_NAMES)
    assert paths["organism_gene_ko"] == tmp_path / "kegg999_organism_gene_ko.tsv.gz"
    assert paths["reference_model"] == tmp_path / "kegg999_reference_model.yml.gz"
