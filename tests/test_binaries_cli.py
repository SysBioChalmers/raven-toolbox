"""Tests for the raven-toolbox-binaries console script (raven_toolbox.binaries_cli)."""
import pytest

from raven_toolbox import binaries, binaries_cli


def test_list_prints_sets(capsys):
    rc = binaries_cli.main(["--list"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "runtime" in out and "build" in out
    assert "hmmsearch" in out and "mafft" in out


def test_unknown_set_errors():
    # argparse parser.error exits with code 2.
    with pytest.raises(SystemExit) as exc:
        binaries_cli.main(["--set", "bogus"])
    assert exc.value.code == 2


def test_set_runtime_reports_and_returns_zero(monkeypatch, capsys):
    def fake_provision(executables, *, prefer_existing=True):
        return [
            binaries.BinaryStatus("blastp", "present", "/usr/bin/blastp"),
            binaries.BinaryStatus("diamond", "downloaded", "/cache/diamond"),
            binaries.BinaryStatus("hmmsearch", "unavailable", "no bundle"),
        ]

    monkeypatch.setattr(binaries, "provision_binaries", fake_provision)
    rc = binaries_cli.main(["--set", "runtime"])
    out = capsys.readouterr().out
    assert rc == 0  # 'unavailable' alone is not a failure (OS limitation)
    assert "blastp" in out and "diamond" in out
    assert "no bundle for" in out  # the conda/WSL2 hint block


def test_download_error_returns_one(monkeypatch, capsys):
    def fake_provision(executables, *, prefer_existing=True):
        return [binaries.BinaryStatus("diamond", "error", "SHA256 mismatch")]

    monkeypatch.setattr(binaries, "provision_binaries", fake_provision)
    rc = binaries_cli.main(["--set", "runtime"])
    out = capsys.readouterr().out
    assert rc == 1
    assert "failed" in out and "SHA256 mismatch" in out
