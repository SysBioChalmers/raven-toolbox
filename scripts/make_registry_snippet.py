#!/usr/bin/env python
"""Emit ready-to-paste registry entries for published artefacts / binary ZIPs.

Computes the SHA256 of each file and prints the Python/JSON entry to merge into
``ravengem.data._DATA_REGISTRY`` (data artefacts) or ``ravengem.binaries._REGISTRY``
(binary bundles). Run once per release, after uploading the files to the release.

Examples
--------
Data artefacts (KEGG reference model + tables + HMM libraries) for one release::

    python scripts/make_registry_snippet.py data \\
        --dataset kegg --version kegg116 --dir artefacts \\
        --base-url https://github.com/ORG/ravengem/releases/download/kegg-data-kegg116

Binary bundle (one ZIP per platform, named ``<bundle>-<version>-<os>-<arch>.zip``)::

    python scripts/make_registry_snippet.py binary \\
        --bundle blast --version 2.16.0 --provides blastp makeblastdb --dir zips \\
        --base-url https://github.com/ORG/ravengem/releases/download/blast-2.16.0

The SHA256 helper is shared with the runtime resolvers (``ravengem.binaries``), so
published checksums always match what ``ensure_data`` / ``ensure_binary`` verify.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from ravengem.binaries import _sha256


def _files_in(directory: Path) -> list[Path]:
    """Regular, non-hidden files in ``directory``, sorted by name."""
    return sorted(p for p in directory.iterdir() if p.is_file() and not p.name.startswith("."))


def data_entry(dataset: str, version: str, base_url: str, directory: Path) -> dict:
    """Build the ``_DATA_REGISTRY[dataset]`` entry for every file in ``directory``."""
    base = base_url.rstrip("/")
    files = {
        p.name: {"url": f"{base}/{p.name}", "sha256": _sha256(p)} for p in _files_in(directory)
    }
    if not files:
        raise SystemExit(f"No files found in {directory}")
    return {"version": version, "files": files}


def binary_entry(
    bundle: str, version: str, provides: list[str], base_url: str, directory: Path
) -> dict:
    """Build the ``_REGISTRY[bundle]`` entry from ``<bundle>-<version>-<os>-<arch>.zip``."""
    base = base_url.rstrip("/")
    prefix = f"{bundle}-{version}-"
    platforms = {}
    for zip_path in directory.glob(f"{prefix}*.zip"):
        platform = zip_path.name[len(prefix) : -len(".zip")]
        platforms[platform] = {"url": f"{base}/{zip_path.name}", "sha256": _sha256(zip_path)}
    if not platforms:
        raise SystemExit(f"No {prefix}*.zip files found in {directory}")
    return {"version": version, "provides": provides, "platforms": dict(sorted(platforms.items()))}


def render(key: str, entry: dict) -> str:
    """Render ``{key: entry}`` as an indented JSON block (valid Python to paste)."""
    return json.dumps({key: entry}, indent=4)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="kind", required=True)

    d = sub.add_parser("data", help="data-artefact registry entry (ravengem.data)")
    d.add_argument("--dataset", required=True, help="dataset key, e.g. 'kegg'")
    d.add_argument("--version", required=True)
    d.add_argument("--dir", required=True, type=Path, help="directory of uploaded artefacts")
    d.add_argument("--base-url", required=True, help="release download URL prefix")

    b = sub.add_parser("binary", help="binary-bundle registry entry (ravengem.binaries)")
    b.add_argument("--bundle", required=True, help="bundle key, e.g. 'blast'")
    b.add_argument("--version", required=True)
    b.add_argument("--provides", nargs="+", required=True, help="executables the bundle provides")
    b.add_argument("--dir", required=True, type=Path, help="directory of uploaded ZIPs")
    b.add_argument("--base-url", required=True, help="release download URL prefix")

    args = parser.parse_args(argv)
    if args.kind == "data":
        key, entry = args.dataset, data_entry(args.dataset, args.version, args.base_url, args.dir)
        target = "ravengem/data.py  _DATA_REGISTRY"
    else:
        key = args.bundle
        entry = binary_entry(args.bundle, args.version, args.provides, args.base_url, args.dir)
        target = "ravengem/binaries.py  _REGISTRY"

    print(f"# Merge into {target}:", file=sys.stderr)
    print(render(key, entry))


if __name__ == "__main__":
    main()
