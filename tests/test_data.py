"""Tests for ensure_data (data.py). Uses file:// URLs to avoid the network."""
import gzip
import hashlib
import subprocess
from pathlib import Path

import pytest

from raven_python import data
from raven_python.data import (
    CORE_KEGG_FILES,
    ensure_data_file,
    ensure_kegg_data,
    ensure_kegg_hmm_library,
    ensure_kegg_taxonomy,
)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


@pytest.fixture
def served(tmp_path, monkeypatch):
    """A fake registry served from local files, with the cache pointed at tmp.

    Published assets are version-prefixed (``v1_<base>``), matching the real
    registry shape that ``ensure_kegg_data`` / ``ensure_kegg_hmm_library`` build.
    ``payloads`` is keyed by the published (prefixed) name.
    """
    src = tmp_path / "src"
    src.mkdir()
    ver = "v1"
    bases = {
        "reference_model.yml.gz": b"!!omap model bytes",
        "ko_reaction.tsv.gz": b"ko\treaction\n",
        "ko_names.tsv.gz": b"ko\tname\n",
        "organism_gene_ko.tsv.gz": b"organism\tgene\tko\n",
        "rxn_flags.tsv.gz": b"reaction\tspontaneous\n",
        "taxonomy.gz": b"# Prokaryotes\n",
    }
    files = {}
    payloads = {}
    for base, data_bytes in bases.items():
        name = f"{ver}_{base}"
        (src / name).write_bytes(data_bytes)
        files[name] = {"url": (src / name).as_uri(), "sha256": _sha256(data_bytes)}
        payloads[name] = data_bytes
    registry = {"kegg": {"version": ver, "files": files}}

    cache = tmp_path / "cache"
    monkeypatch.setenv("XDG_CACHE_HOME", str(cache))
    return registry, cache, payloads


def test_ensure_data_file_downloads_and_caches(served):
    registry, cache, payloads = served
    path = ensure_data_file("kegg", "v1_ko_reaction.tsv.gz", registry=registry)
    assert path == cache / "raven_python" / "data" / "kegg-v1" / "v1_ko_reaction.tsv.gz"
    assert path.read_bytes() == payloads["v1_ko_reaction.tsv.gz"]


def test_ensure_data_file_reuses_cache(served):
    registry, _, _ = served
    first = ensure_data_file("kegg", "v1_ko_names.tsv.gz", registry=registry)
    # Break the URL: a second call must hit the cache, not re-download.
    registry["kegg"]["files"]["v1_ko_names.tsv.gz"]["url"] = "file:///nonexistent"
    second = ensure_data_file("kegg", "v1_ko_names.tsv.gz", registry=registry)
    assert first == second and second.exists()


def test_sha256_mismatch_rejected(served):
    registry, cache, _ = served
    registry["kegg"]["files"]["v1_rxn_flags.tsv.gz"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        ensure_data_file("kegg", "v1_rxn_flags.tsv.gz", registry=registry)
    # The corrupt partial download must not be left behind.
    assert not (cache / "raven_python" / "data" / "kegg-v1" / "v1_rxn_flags.tsv.gz").exists()


def test_unknown_dataset_actionable_error(served):
    registry, _, _ = served
    with pytest.raises(FileNotFoundError, match="No data artefacts registered"):
        ensure_data_file("metacyc", "x", registry=registry)


def test_unknown_file_lists_available(served):
    registry, _, _ = served
    with pytest.raises(FileNotFoundError, match="not registered"):
        ensure_data_file("kegg", "missing.tsv.gz", registry=registry)


def test_ensure_kegg_data_fetches_core_set(served):
    registry, cache, _ = served
    out = ensure_kegg_data(registry=registry)
    assert out == cache / "raven_python" / "data" / "kegg-v1"
    # CORE_KEGG_FILES are base names; the cached files are version-prefixed.
    for base in CORE_KEGG_FILES:
        assert (out / f"v1_{base}").is_file()


def test_ensure_kegg_hmm_library_decompresses_and_presses(served, tmp_path, monkeypatch):
    registry, cache, _ = served
    raw = b"HMMER3/f [3.4]\nNAME  K00001\n//\n"
    blob = gzip.compress(raw, mtime=0)
    gz = tmp_path / "src" / "v1_prokaryotes.hmm.gz"
    gz.write_bytes(blob)
    registry["kegg"]["files"]["v1_prokaryotes.hmm.gz"] = {
        "url": gz.as_uri(),
        "sha256": _sha256(blob),
    }

    presses: list[Path] = []

    def fake_run(cmd, capture_output, text):
        library = Path(cmd[-1])
        presses.append(library)
        for suffix in (".h3f", ".h3i", ".h3m", ".h3p"):
            library.with_name(library.name + suffix).write_bytes(b"")
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(data, "resolve_binary", lambda name, binary=None: "hmmpress")
    monkeypatch.setattr(data.subprocess, "run", fake_run)

    library = ensure_kegg_hmm_library("prokaryotes", registry=registry)
    assert library.name == "v1_prokaryotes.hmm"
    assert library.read_bytes() == raw  # decompressed in place
    assert library.with_name(library.name + ".h3m").exists()
    assert len(presses) == 1

    # Second call: library + sidecars already cached, so no re-decompress/-press.
    again = ensure_kegg_hmm_library("prokaryotes", registry=registry)
    assert again == library
    assert len(presses) == 1


def test_ensure_kegg_taxonomy(served):
    registry, cache, _ = served
    path = ensure_kegg_taxonomy(registry=registry)
    assert path == cache / "raven_python" / "data" / "kegg-v1" / "v1_taxonomy.gz"
    assert path.is_file()


def test_unregistered_dataset_raises():
    # An unpublished dataset still raises an actionable error against the shipped registry.
    with pytest.raises(FileNotFoundError, match="No data artefacts registered"):
        ensure_data_file("metacyc", "x")


def test_shipped_registry_matches_resolver_names():
    # The published registry keys must equal what ensure_kegg_data /
    # ensure_kegg_hmm_library construct as f"{version}_{base}", or fetches 404.
    kegg = data._DATA_REGISTRY["kegg"]
    ver = kegg["version"]
    names = set(kegg["files"])
    assert {f"{ver}_{base}" for base in CORE_KEGG_FILES} <= names
    assert {f"{ver}_{d}.hmm.gz" for d in ("prokaryotes", "eukaryotes")} <= names
    assert f"{ver}_taxonomy.gz" in names
