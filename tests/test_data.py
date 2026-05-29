"""Tests for ensure_data (data.py). Uses file:// URLs to avoid the network."""
import hashlib

import pytest

from raven_python.data import (
    CORE_KEGG_FILES,
    ensure_data_file,
    ensure_kegg_data,
)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


@pytest.fixture
def served(tmp_path, monkeypatch):
    """A fake registry served from local files, with the cache pointed at tmp."""
    src = tmp_path / "src"
    src.mkdir()
    payloads = {
        "reference_model.yml.gz": b"!!omap model bytes",
        "ko_reaction.tsv.gz": b"ko\treaction\n",
        "ko_names.tsv.gz": b"ko\tname\n",
        "organism_gene_ko.tsv.xz": b"organism\tgene\tko\n",
        "rxn_flags.tsv.gz": b"reaction\tspontaneous\n",
    }
    files = {}
    for name, data in payloads.items():
        path = src / name
        path.write_bytes(data)
        files[name] = {"url": path.as_uri(), "sha256": _sha256(data)}
    registry = {"kegg": {"version": "v1", "files": files}}

    cache = tmp_path / "cache"
    monkeypatch.setenv("XDG_CACHE_HOME", str(cache))
    return registry, cache, payloads


def test_ensure_data_file_downloads_and_caches(served):
    registry, cache, payloads = served
    path = ensure_data_file("kegg", "ko_reaction.tsv.gz", registry=registry)
    assert path == cache / "raven_python" / "data" / "kegg-v1" / "ko_reaction.tsv.gz"
    assert path.read_bytes() == payloads["ko_reaction.tsv.gz"]


def test_ensure_data_file_reuses_cache(served, monkeypatch):
    registry, _, _ = served
    first = ensure_data_file("kegg", "ko_names.tsv.gz", registry=registry)
    # Break the URL: a second call must hit the cache, not re-download.
    registry["kegg"]["files"]["ko_names.tsv.gz"]["url"] = "file:///nonexistent"
    second = ensure_data_file("kegg", "ko_names.tsv.gz", registry=registry)
    assert first == second and second.exists()


def test_sha256_mismatch_rejected(served):
    registry, cache, _ = served
    registry["kegg"]["files"]["rxn_flags.tsv.gz"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        ensure_data_file("kegg", "rxn_flags.tsv.gz", registry=registry)
    # The corrupt partial download must not be left behind.
    assert not (cache / "raven_python" / "data" / "kegg-v1" / "rxn_flags.tsv.gz").exists()


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
    for name in CORE_KEGG_FILES:
        assert (out / name).is_file()


def test_empty_registry_raises():
    # The shipped registry is empty until artefacts are published.
    with pytest.raises(FileNotFoundError, match="No data artefacts registered"):
        ensure_data_file("kegg", "ko_reaction.tsv.gz")
