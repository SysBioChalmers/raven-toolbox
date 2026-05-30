"""Tests for the data/binary manifest loader (manifest.py) and its wiring into the
resolvers. Uses file:// URLs + a tmp manifest to avoid the network."""
import hashlib
import json
from pathlib import Path

import pytest

from raven_python import binaries, data, manifest

REPO_ROOT = Path(__file__).resolve().parents[1]


def _sha256(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


@pytest.fixture
def clean_registries():
    """Snapshot and restore the live registries so a test's loads don't leak."""
    data_snap = dict(data._DATA_REGISTRY)
    bin_snap = dict(binaries._REGISTRY)
    data._DATA_REGISTRY.clear()
    binaries._REGISTRY.clear()
    yield
    data._DATA_REGISTRY.clear()
    data._DATA_REGISTRY.update(data_snap)
    binaries._REGISTRY.clear()
    binaries._REGISTRY.update(bin_snap)


def _write_manifest(tmp_path: Path, payload: dict) -> Path:
    p = tmp_path / "manifest.json"
    p.write_text(json.dumps(payload), encoding="utf-8")
    return p


def test_converters_strip_to_registry_shapes():
    m = {
        "manifest_version": 1,
        "data": {
            "kegg": {
                "version": "kegg116",
                "license": "metadata that the registry should drop",
                "files": {"x.gz": {"url": "https://e/x.gz", "sha256": "a" * 64, "bytes": 5}},
            }
        },
        "binaries": {
            "diamond": {
                "version": "2.1.9",
                "provides": ["diamond"],
                "license": "GPL-3.0-only",
                "platforms": {"linux-x86_64": {"url": "https://e/d.zip", "sha256": "b" * 64, "bytes": 9}},
            }
        },
    }
    assert manifest.to_data_registry(m) == {
        "kegg": {"version": "kegg116", "files": {"x.gz": {"url": "https://e/x.gz", "sha256": "a" * 64}}}
    }
    assert manifest.to_binary_registry(m) == {
        "diamond": {
            "version": "2.1.9",
            "provides": ["diamond"],
            "platforms": {"linux-x86_64": {"url": "https://e/d.zip", "sha256": "b" * 64}},
        }
    }


def test_load_manifest_rejects_unknown_version(tmp_path):
    p = _write_manifest(tmp_path, {"manifest_version": 999})
    with pytest.raises(ValueError, match="manifest_version"):
        manifest.load_manifest(p)


def test_load_manifest_requires_a_source(monkeypatch):
    monkeypatch.delenv(manifest.ENV_MANIFEST, raising=False)
    with pytest.raises(ValueError, match="No manifest source"):
        manifest.load_manifest()


def test_load_into_registries_populates_both(tmp_path, clean_registries):
    p = _write_manifest(
        tmp_path,
        {
            "manifest_version": 1,
            "data": {"kegg": {"version": "v1", "files": {"a": {"url": "https://e/a", "sha256": "0" * 64}}}},
            "binaries": {"diamond": {"version": "2", "provides": ["diamond"], "platforms": {}}},
        },
    )
    manifest.load_into_registries(p)
    assert data._DATA_REGISTRY["kegg"]["version"] == "v1"
    assert binaries._REGISTRY["diamond"]["provides"] == ["diamond"]


def test_resolver_lazy_autoload_via_env(tmp_path, monkeypatch, clean_registries):
    # A real artefact file served over file://, registered through the manifest.
    artefact = tmp_path / "reference_model.yml.gz"
    artefact.write_bytes(b"hello kegg")
    payload = {
        "manifest_version": 1,
        "data": {
            "kegg": {
                "version": "kegg-test",
                "files": {
                    artefact.name: {"url": artefact.as_uri(), "sha256": _sha256(artefact.read_bytes())}
                },
            }
        },
    }
    manifest_path = _write_manifest(tmp_path, payload)
    monkeypatch.setenv(manifest.ENV_MANIFEST, str(manifest_path))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))

    # _DATA_REGISTRY is empty; ensure_data_file must lazily load the manifest and fetch.
    got = data.ensure_data_file("kegg", artefact.name)
    assert got.read_bytes() == b"hello kegg"
    assert data._DATA_REGISTRY["kegg"]["version"] == "kegg-test"


@pytest.mark.parametrize("name", ["manifest.json", "manifest.example.json"])
def test_repo_manifests_are_valid(name):
    m = manifest.load_manifest(REPO_ROOT / "data" / name)
    # Both must convert cleanly to the runtime registry shapes.
    manifest.to_data_registry(m)
    manifest.to_binary_registry(m)
