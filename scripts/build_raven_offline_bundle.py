#!/usr/bin/env python
"""Build per-platform offline binary bundles for RAVEN (air-gapped install).

RAVEN normally downloads its external binaries (BLAST+/DIAMOND/HMMER/WoLFPSORT)
on demand from the raven-data repository. For an offline or air-gapped machine,
this produces **one ZIP per platform** whose contents drop straight into RAVEN's
``software/`` directory, so the binaries are present without any network access.

It fetches the published per-platform ZIPs from raven-data and repacks them under
``<tool>/`` prefixes — preserving the execute bits recorded in the source ZIPs'
metadata, so a Linux/macOS ``unzip`` yields runnable binaries without RAVEN's
on-demand ``chmod`` step. Writes ``RAVEN-binaries-<os>-<arch>.zip`` (+ SHA256s and
an INSTALL note) to ``--out``.

Publish: attach each ZIP as an asset on the matching RAVEN release. An offline
user then unzips the one for their platform into ``<RAVEN>/software/``::

    unzip RAVEN-binaries-linux-x86_64.zip -d /path/to/RAVEN/software

Usage::

    python scripts/build_raven_offline_bundle.py                       # all platforms
    python scripts/build_raven_offline_bundle.py --platform linux-x86_64
"""
from __future__ import annotations

import argparse
import hashlib
import zipfile
from pathlib import Path
from urllib.request import urlopen

RAVEN_DATA = "https://github.com/SysBioChalmers/raven-data/releases/download"
PLATFORMS = ("linux-x86_64", "macos-arm64", "windows-x86_64")


def _tools(plat: str) -> list[tuple[str, str, str]]:
    """(asset filename, release tag, destination prefix inside software/) per tool."""
    hmmer = "3.3.2" if plat == "windows-x86_64" else "3.4.0"  # native Windows HMMER is 3.3.2
    return [
        (f"blast-2.17.0-{plat}.zip",   "blast-2.17.0",   "blast+/"),
        (f"diamond-2.1.17-{plat}.zip", "diamond-2.1.17", "diamond/"),
        (f"hmmer-{hmmer}-{plat}.zip",  f"hmmer-{hmmer}",  "hmmer/"),
        ("wolfpsort-0.2.zip",          "wolfpsort-0.2",   ""),  # zip already holds WoLFPSORT/
    ]


def _download(tag: str, asset: str, cache: Path) -> Path:
    """Fetch a raven-data asset into the cache (skip if already complete)."""
    dest = cache / asset
    url = f"{RAVEN_DATA}/{tag}/{asset}"
    with urlopen(url, timeout=120) as resp:  # noqa: S310 (pinned raven-data URLs)
        total = int(resp.headers.get("Content-Length", 0)) or None
        if dest.exists() and total is not None and dest.stat().st_size == total:
            return dest
        data = resp.read()
    if total is not None and len(data) != total:
        raise RuntimeError(f"short read for {url}: {len(data)} of {total} bytes")
    dest.write_bytes(data)
    return dest


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build(out: Path, platforms: tuple[str, ...]) -> None:
    cache = out / "_src"
    cache.mkdir(parents=True, exist_ok=True)
    rows: list[tuple[str, str, int]] = []

    for plat in platforms:
        bundle = out / f"RAVEN-binaries-{plat}.zip"
        print(f"  building {bundle.name} ...")
        with zipfile.ZipFile(bundle, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as dst:
            for asset, tag, prefix in _tools(plat):
                src_path = _download(tag, asset, cache)
                with zipfile.ZipFile(src_path) as src:
                    for info in src.infolist():
                        if info.is_dir():
                            continue
                        zi = zipfile.ZipInfo(prefix + info.filename, date_time=(1980, 1, 1, 0, 0, 0))
                        zi.external_attr = info.external_attr  # preserve mode (exec bits)
                        zi.compress_type = zipfile.ZIP_DEFLATED
                        dst.writestr(zi, src.read(info))
        sha, size = _sha256(bundle), bundle.stat().st_size
        rows.append((bundle.name, sha, size))

    (out / "checksums.txt").write_text(
        "".join(f"{sha}  {name}  ({size} bytes)\n" for name, sha, size in rows), encoding="utf-8"
    )
    (out / "INSTALL.txt").write_text(
        "Offline RAVEN binaries.\n\n"
        "On the target machine, unzip the file for your platform into RAVEN's\n"
        "software/ directory, e.g.:\n\n"
        "    unzip RAVEN-binaries-linux-x86_64.zip -d /path/to/RAVEN/software\n\n"
        "RAVEN then finds the binaries already present and never tries to download\n"
        "them. (WoLFPSORT is a Linux build; on Windows it runs via WSL.)\n",
        encoding="utf-8",
    )

    print(f"\n{len(rows)} bundle(s) written to {out}")
    for name, sha, size in rows:
        print(f"  {size/1e6:7.1f} MB  {sha[:12]}…  {name}")
    print(f"\nSHA256s: {out / 'checksums.txt'}\nInstall note: {out / 'INSTALL.txt'}")
    print("\nUpload each ZIP as an asset on the matching RAVEN release.")


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", type=Path, default=Path("dist/raven-offline"), help="output folder")
    p.add_argument("--platform", choices=PLATFORMS, action="append",
                   help="limit to this platform (repeatable); default is all")
    args = p.parse_args(argv)
    build(args.out, tuple(args.platform) if args.platform else PLATFORMS)


if __name__ == "__main__":
    main()
