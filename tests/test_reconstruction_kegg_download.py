"""Tests for the KEGG download/arrange tooling (reconstruction/kegg/download.py).

The network fetch needs a paid KEGG subscription, so it is not exercised here.
We test credential resolution and the network-free extract/arrange core against
hand-built fake archives.
"""
import gzip
import io
import tarfile
from pathlib import Path

import pytest

from raven_python.reconstruction.kegg.download import (
    _resolve_auth,
    extract_kegg_dump,
)


def _make_targz(path: Path, members: dict[str, bytes]) -> None:
    with tarfile.open(path, "w:gz") as tar:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))


def _make_gz(path: Path, data: bytes) -> None:
    with gzip.open(path, "wb") as fh:
        fh.write(data)


# --------------------------------------------------------------------------- #
# Credentials
# --------------------------------------------------------------------------- #
def test_resolve_auth_explicit_wins():
    assert _resolve_auth("ftp.kegg.net", auth=("u", "p")) == ("u", "p")


def test_resolve_auth_from_netrc(tmp_path):
    netrc_file = tmp_path / ".netrc"
    netrc_file.write_text("machine ftp.kegg.net login alice password s3cret\n")
    netrc_file.chmod(0o600)
    assert _resolve_auth("ftp.kegg.net", netrc_path=netrc_file) == ("alice", "s3cret")


def test_resolve_auth_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError, match="does not exist"):
        _resolve_auth("ftp.kegg.net", netrc_path=tmp_path / "nope")


def test_resolve_auth_host_absent(tmp_path):
    netrc_file = tmp_path / ".netrc"
    netrc_file.write_text("machine other.host login a password b\n")
    netrc_file.chmod(0o600)
    with pytest.raises(ValueError, match="No credentials for"):
        _resolve_auth("ftp.kegg.net", netrc_path=netrc_file)


def test_resolve_auth_malformed_netrc(tmp_path):
    netrc_file = tmp_path / ".netrc"
    netrc_file.write_text("this is not a valid netrc line\n")
    netrc_file.chmod(0o600)
    with pytest.raises(ValueError, match="Could not read credentials"):
        _resolve_auth("ftp.kegg.net", netrc_path=netrc_file)


# --------------------------------------------------------------------------- #
# Extract / arrange
# --------------------------------------------------------------------------- #
@pytest.fixture
def fake_dump(tmp_path):
    """A tmp dir populated with fake KEGG archives, as fetch would leave them."""
    _make_targz(
        tmp_path / "reaction.tar.gz",
        {
            "reaction/reaction": b"RXN_ENTRIES\n",
            "reaction/reaction.lst": b"R00010: A <=> B\n",
            "reaction/reaction_mapformula.lst": b"R00010: 00010: A => B\n",
            "reaction/reaction.name": b"discard me\n",  # extra file, not lifted
        },
    )
    _make_targz(
        tmp_path / "compound.tar.gz",
        {"compound/compound": b"CPD\n", "compound/compound.inchi": b"C00031\tInChI=x\n"},
    )
    _make_targz(tmp_path / "glycan.tar.gz", {"glycan/glycan": b"GLY\n"})
    _make_targz(tmp_path / "ko.tar.gz", {"ko/ko": b"KO\n"})
    _make_gz(tmp_path / "eukaryotes.pep.gz", b">euk\nMKV\n")
    _make_gz(tmp_path / "prokaryotes.pep.gz", b">prok\nMAA\n")
    (tmp_path / "taxonomy").write_text("tax\n")
    return tmp_path


def test_extract_produces_flat_layout(fake_dump):
    result = extract_kegg_dump(fake_dump)
    expected = {
        "reaction",
        "reaction.lst",
        "reaction_mapformula.lst",
        "compound",
        "compound.inchi",
        "ko",
        "genes.pep",
        "taxonomy",
    }
    assert set(result) == expected
    assert all(p.is_file() for p in result.values())


def test_extract_concatenates_compound_and_glycan(fake_dump):
    extract_kegg_dump(fake_dump)
    assert (fake_dump / "compound").read_bytes() == b"CPD\nGLY\n"


def test_extract_concatenates_proteomes(fake_dump):
    extract_kegg_dump(fake_dump)
    assert (fake_dump / "genes.pep").read_bytes() == b">euk\nMKV\n>prok\nMAA\n"


def test_extract_removes_subdirs_and_archives(fake_dump):
    extract_kegg_dump(fake_dump)
    assert not list(fake_dump.glob("*.tar.gz"))
    assert not list(fake_dump.glob("*.gz"))
    for subdir in ("reaction", "compound", "glycan", "ko"):
        assert not (fake_dump / subdir).is_dir()
    assert not (fake_dump / "reaction.name").exists()  # extra file discarded


def test_extract_requires_core_archives(tmp_path):
    _make_targz(tmp_path / "compound.tar.gz", {"compound/compound": b"CPD\n"})
    with pytest.raises(FileNotFoundError, match="required file"):
        extract_kegg_dump(tmp_path)
