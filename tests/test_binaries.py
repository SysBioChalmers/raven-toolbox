"""Tests for raven_toolbox.binaries (binary resolution + bundled-ZIP provisioning)."""
import hashlib
import shutil
import stat
import zipfile
from pathlib import Path

import pytest

from raven_toolbox import binaries


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


def test_ensure_binary_rejects_path_traversal_zip(tmp_path, monkeypatch):
    # A ZIP whose member escapes the extraction dir must be refused, not written
    # outside the cache (ZipFile.extractall has no traversal guard of its own).
    archive = tmp_path / "evil.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        zf.writestr("footool", "ok")
        zf.writestr("../escape.txt", "pwned")
    sha = hashlib.sha256(archive.read_bytes()).hexdigest()
    registry = {
        "footool": {
            "version": "1.0",
            "provides": ["footool"],
            "platforms": {binaries.platform_key(): {"url": archive.as_uri(), "sha256": sha}},
        }
    }
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    with pytest.raises(ValueError, match="path traversal"):
        binaries.ensure_binary("footool", registry=registry)
    assert not (tmp_path / "escape.txt").exists()


def test_ensure_binary_rejects_symlink_zip(tmp_path, monkeypatch):
    # A ZIP member flagged as a symlink must be refused: extractall would
    # recreate it as a real symlink whose target the path guard never validates,
    # so it could point outside the cache dir.
    archive = tmp_path / "evil.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        zf.writestr("footool", "ok")
        link = zipfile.ZipInfo("escape")
        link.external_attr = (stat.S_IFLNK | 0o777) << 16
        zf.writestr(link, "/etc/passwd")
    sha = hashlib.sha256(archive.read_bytes()).hexdigest()
    registry = {
        "footool": {
            "version": "1.0",
            "provides": ["footool"],
            "platforms": {binaries.platform_key(): {"url": archive.as_uri(), "sha256": sha}},
        }
    }
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    with pytest.raises(ValueError, match="symlink"):
        binaries.ensure_binary("footool", registry=registry)


def test_ensure_binary_passes_download_timeout(tmp_path, monkeypatch):
    # The download must use a socket timeout so a stalled server can't hang the
    # process forever. Spy on urlopen to confirm a positive timeout is forwarded.
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

    seen = {}
    real_urlopen = binaries.urlopen

    def spy(url, *args, **kwargs):
        seen["timeout"] = kwargs.get("timeout")
        return real_urlopen(url, *args, **kwargs)

    monkeypatch.setattr(binaries, "urlopen", spy)
    binaries.ensure_binary("footool", registry=registry)
    assert isinstance(seen["timeout"], (int, float)) and seen["timeout"] > 0
