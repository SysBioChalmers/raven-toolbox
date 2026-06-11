"""Tests for ensure_data (data.py). Uses file:// URLs to avoid the network."""
import gzip
import hashlib
import io
import tarfile

import pytest

from raven_python import data
from raven_python.data import (
    ensure_data_file,
    ensure_kegg_data,
    ensure_kegg_hmm_library,
    ensure_kegg_taxonomy,
)


def _sha256(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


@pytest.fixture
def served(tmp_path, monkeypatch):
    """A fake registry served from local files, with the cache pointed at tmp.

    The core artefacts are delivered as a single version-prefixed bundle
    (``v1_core.tar.gz``) that ``ensure_kegg_data`` extracts; taxonomy is a separate
    file. ``core`` maps each bundled member name -> its raw payload.
    """
    src = tmp_path / "src"
    src.mkdir()
    ver = "v1"
    core = {
        f"{ver}_reference_model.yml.gz": b"!!omap model bytes",
        f"{ver}_ko_reaction.tsv.gz": b"ko\treaction\n",
        f"{ver}_ko_names.tsv.gz": b"ko\tname\n",
        f"{ver}_organism_gene_ko.tsv.gz": b"organism\tgene\tko\n",
        f"{ver}_rxn_flags.tsv.gz": b"reaction\tspontaneous\n",
    }
    bundle = src / f"{ver}_core.tar.gz"
    with tarfile.open(bundle, "w:gz") as tar:
        for name, payload in core.items():
            info = tarfile.TarInfo(name=name)
            info.size = len(payload)
            tar.addfile(info, io.BytesIO(payload))
    taxonomy = src / f"{ver}_taxonomy.gz"
    taxonomy.write_bytes(b"# Prokaryotes\n")

    files = {
        f"{ver}_core.tar.gz": {"url": bundle.as_uri(), "sha256": _sha256(bundle.read_bytes())},
        f"{ver}_taxonomy.gz": {"url": taxonomy.as_uri(), "sha256": _sha256(taxonomy.read_bytes())},
    }
    registry = {"kegg": {"version": ver, "files": files}}

    cache = tmp_path / "cache"
    monkeypatch.setenv("XDG_CACHE_HOME", str(cache))
    return registry, cache, core


def test_ensure_data_file_downloads_and_caches(served):
    registry, cache, _ = served
    path = ensure_data_file("kegg", "v1_taxonomy.gz", registry=registry)
    assert path == cache / "raven_python" / "data" / "kegg-v1" / "v1_taxonomy.gz"
    assert path.read_bytes() == b"# Prokaryotes\n"


def test_ensure_data_file_reuses_cache(served):
    registry, _, _ = served
    first = ensure_data_file("kegg", "v1_taxonomy.gz", registry=registry)
    # Break the URL: a second call must hit the cache, not re-download.
    registry["kegg"]["files"]["v1_taxonomy.gz"]["url"] = "file:///nonexistent"
    second = ensure_data_file("kegg", "v1_taxonomy.gz", registry=registry)
    assert first == second and second.exists()


def test_sha256_mismatch_rejected(served):
    registry, cache, _ = served
    registry["kegg"]["files"]["v1_taxonomy.gz"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        ensure_data_file("kegg", "v1_taxonomy.gz", registry=registry)
    # The corrupt partial download must not be left behind.
    assert not (cache / "raven_python" / "data" / "kegg-v1" / "v1_taxonomy.gz").exists()


def test_unknown_dataset_actionable_error(served):
    registry, _, _ = served
    with pytest.raises(FileNotFoundError, match="No data artefacts registered"):
        ensure_data_file("metacyc", "x", registry=registry)


def test_unknown_file_lists_available(served):
    registry, _, _ = served
    with pytest.raises(FileNotFoundError, match="not registered"):
        ensure_data_file("kegg", "missing.tsv.gz", registry=registry)


def test_ensure_kegg_data_extracts_core_bundle(served):
    registry, cache, core = served
    out = ensure_kegg_data(registry=registry)
    assert out == cache / "raven_python" / "data" / "kegg-v1"
    # The single bundle is fetched and its members extracted into the cache dir.
    for name, payload in core.items():
        member = out / name
        assert member.is_file() and member.read_bytes() == payload


def test_ensure_kegg_taxonomy(served):
    registry, cache, _ = served
    path = ensure_kegg_taxonomy(registry=registry)
    assert path == cache / "raven_python" / "data" / "kegg-v1" / "v1_taxonomy.gz"
    assert path.is_file()


def test_ensure_kegg_hmm_library_decompresses(served, tmp_path):
    registry, cache, _ = served
    raw = b"HMMER3/f [3.4]\nNAME  K00001\n//\n"
    blob = gzip.compress(raw, mtime=0)
    gz = tmp_path / "src" / "v1_prokaryotes.hmm.gz"
    gz.write_bytes(blob)
    registry["kegg"]["files"]["v1_prokaryotes.hmm.gz"] = {
        "url": gz.as_uri(),
        "sha256": _sha256(blob),
    }
    library = ensure_kegg_hmm_library("prokaryotes", registry=registry)
    assert library.name == "v1_prokaryotes.hmm"
    assert library.read_bytes() == raw  # decompressed flatfile, no hmmpress

    # Second call: already decompressed, returns the same cached library.
    assert ensure_kegg_hmm_library("prokaryotes", registry=registry) == library


def test_unregistered_dataset_raises():
    # An unpublished dataset still raises an actionable error against the shipped registry.
    with pytest.raises(FileNotFoundError, match="No data artefacts registered"):
        ensure_data_file("metacyc", "x")


def test_shipped_registry_has_expected_assets():
    # The published registry holds the core bundle, both HMM libraries, and taxonomy.
    kegg = data._DATA_REGISTRY["kegg"]
    ver = kegg["version"]
    names = set(kegg["files"])
    assert f"{ver}_core.tar.gz" in names
    assert {f"{ver}_{d}.hmm.gz" for d in ("prokaryotes", "eukaryotes")} <= names
    assert f"{ver}_taxonomy.gz" in names
