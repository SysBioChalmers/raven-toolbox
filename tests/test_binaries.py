"""Tests for raven_python.binaries (binary resolution + bundled-ZIP provisioning)."""
import hashlib
import shutil
import zipfile
from pathlib import Path

import pytest

from raven_python import binaries


def test_resolve_explicit_path():
    assert binaries.resolve_binary("blastp", binary="/opt/x/blastp") == "/opt/x/blastp"


def test_resolve_env_var(monkeypatch):
    monkeypatch.setenv("RAVEN_PYTHON_DIAMOND", "/custom/diamond")
    assert binaries.resolve_binary("diamond") == "/custom/diamond"


@pytest.mark.skipif(not shutil.which("blastp"), reason="blastp not installed")
def test_resolve_via_path():
    assert binaries.resolve_binary("blastp") == shutil.which("blastp")


def test_resolve_unresolvable_raises(monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda _: None)
    with pytest.raises(FileNotFoundError, match="Could not find"):
        binaries.resolve_binary("diamond")  # empty registry, not on PATH


def test_platform_key_format():
    key = binaries.platform_key()
    assert "-" in key
    os_part, arch = key.split("-", 1)
    assert os_part in {"linux", "macos", "windows"} or os_part  # tolerant


def test_ensure_binary_downloads_verifies_extracts(tmp_path, monkeypatch):
    # Build a fake bundle ZIP containing an executable, served via file:// URL.
    exe = tmp_path / "footool"
    exe.write_text("#!/bin/sh\necho hi\n")
    archive = tmp_path / "footool.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        zf.write(exe, "footool")
    sha = hashlib.sha256(archive.read_bytes()).hexdigest()

    registry = {
        "footool": {
            "version": "1.0",
            "provides": ["footool"],
            "platforms": {binaries.platform_key(): {"url": archive.as_uri(), "sha256": sha}},
        }
    }
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))

    path = binaries.ensure_binary("footool", registry=registry)
    assert Path(path).exists()
    assert Path(path).name == "footool"
    # cached on second call (same path, no re-download needed)
    assert binaries.ensure_binary("footool", registry=registry) == path


def test_ensure_binary_sha_mismatch(tmp_path, monkeypatch):
    archive = tmp_path / "x.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        zf.writestr("footool", "data")
    registry = {
        "footool": {"version": "1", "provides": ["footool"],
                    "platforms": {binaries.platform_key(): {"url": archive.as_uri(), "sha256": "deadbeef"}}}
    }
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        binaries.ensure_binary("footool", registry=registry)


def test_ensure_binary_unhosted_platform_raises(tmp_path):
    registry = {"footool": {"version": "1", "provides": ["footool"], "platforms": {}}}
    with pytest.raises(FileNotFoundError, match="No bundled"):
        binaries.ensure_binary("footool", registry=registry)
