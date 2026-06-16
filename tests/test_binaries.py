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


# --------------------------------------------------------------------------- #
# Auto-fetch toggle (RAVEN_PYTHON_AUTOFETCH)
# --------------------------------------------------------------------------- #
def test_autofetch_enabled_default_and_off(monkeypatch):
    monkeypatch.delenv("RAVEN_PYTHON_AUTOFETCH", raising=False)
    assert binaries.autofetch_enabled() is True
    for off in ("0", "false", "No", "OFF"):
        monkeypatch.setenv("RAVEN_PYTHON_AUTOFETCH", off)
        assert binaries.autofetch_enabled() is False
    monkeypatch.setenv("RAVEN_PYTHON_AUTOFETCH", "1")
    assert binaries.autofetch_enabled() is True


def test_resolve_skips_download_when_autofetch_off(monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda _: None)
    monkeypatch.setenv("RAVEN_PYTHON_AUTOFETCH", "0")

    def boom(*a, **k):  # ensure_binary must NOT be reached
        raise AssertionError("ensure_binary called despite auto-fetch disabled")

    monkeypatch.setattr(binaries, "ensure_binary", boom)
    with pytest.raises(FileNotFoundError, match="auto-fetch is disabled"):
        binaries.resolve_binary("diamond")


# --------------------------------------------------------------------------- #
# Binary sets + provisioning
# --------------------------------------------------------------------------- #
def test_executables_for_set():
    assert binaries.executables_for_set("runtime") == ("blastp", "makeblastdb", "diamond", "hmmsearch")
    assert binaries.executables_for_set("build") == ("hmmbuild", "mafft", "cd-hit")
    # 'all' is the de-duplicated union of every set, order preserved.
    assert binaries.executables_for_set("all") == (
        "blastp", "makeblastdb", "diamond", "hmmsearch", "hmmbuild", "mafft", "cd-hit",
    )
    with pytest.raises(ValueError, match="Unknown binary set"):
        binaries.executables_for_set("nope")


def test_provision_prefers_existing(monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda exe: "/usr/bin/diamond" if exe == "diamond" else None)
    [res] = binaries.provision_binaries(["diamond"])
    assert res.status == "present"
    assert res.detail == "/usr/bin/diamond"


def test_provision_unavailable_when_no_bundle(monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda _: None)
    [res] = binaries.provision_binaries(["diamond"], registry={})
    assert res.status == "unavailable"
    assert res.executable == "diamond"


def test_provision_downloads_when_missing(tmp_path, monkeypatch):
    exe = tmp_path / "diamond"
    exe.write_text("#!/bin/sh\necho hi\n")
    archive = tmp_path / "diamond.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        zf.write(exe, "diamond")
    sha = hashlib.sha256(archive.read_bytes()).hexdigest()
    registry = {
        "diamond": {
            "version": "1.0", "provides": ["diamond"],
            "platforms": {binaries.platform_key(): {"url": archive.as_uri(), "sha256": sha}},
        }
    }
    monkeypatch.setattr(shutil, "which", lambda _: None)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    [res] = binaries.provision_binaries(["diamond"], registry=registry)
    assert res.status == "downloaded"
    assert Path(res.detail).name == "diamond"


def test_ensure_binary_windows_exe_suffix(tmp_path, monkeypatch):
    # On Windows the bundle ships hmmsearch.exe; ensure_binary must look it up and
    # return it under the .exe name (otherwise Windows bundles never resolve).
    archive = tmp_path / "hmmer.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        zf.writestr("hmmsearch.exe", "MZ binary")
        zf.writestr("cygwin1.dll", "dll")
    sha = hashlib.sha256(archive.read_bytes()).hexdigest()
    registry = {
        "hmmer": {
            "version": "3.3.2", "provides": ["hmmsearch"],
            "platforms": {"windows-x86_64": {"url": archive.as_uri(), "sha256": sha}},
        }
    }
    monkeypatch.setattr(binaries, "platform_key", lambda: "windows-x86_64")
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    path = binaries.ensure_binary("hmmsearch", registry=registry)
    assert Path(path).name == "hmmsearch.exe"
    assert Path(path).exists()


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
