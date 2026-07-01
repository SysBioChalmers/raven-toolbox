#!/usr/bin/env python3
"""Build the transporter reference databases for evidence-aware transport scoring.

Produces two raven-data artefacts (host them under a ``transporters-<date>`` release, then run
``make_registry_snippet.py`` + ``publish_to_raven_data.py`` per docs/maintenance/artefact_hosting.md):

* ``transporter_pfam.hmm`` — the HMMER database of the transporter Pfam families in
  :data:`raven_toolbox.localization.transporter_tables.PFAM_TRANSPORTERS`, fetched one HMM at a time
  from the InterPro API and concatenated (used directly by ``hmmsearch``; no ``hmmpress`` needed).
* ``tcdb.dmnd`` — a DIAMOND database of the TCDB protein sequences.

Every family accession is verified against the live InterPro API; a stale/typo'd accession is reported
and skipped rather than silently shipped. Needs network access (InterPro + TCDB) and the bundled
``hmmpress`` / ``diamond`` binaries (auto-provisioned via ``raven_toolbox.binaries``).

Usage: ``python scripts/build_transporter_data.py --out dist/transporters [--only-pfam|--only-tcdb]``
ASCII-only output.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import subprocess
import sys
import urllib.request
from pathlib import Path

from raven_toolbox.binaries import resolve_binary
from raven_toolbox.localization.transporter_tables import PFAM_TRANSPORTERS

INTERPRO_HMM = "https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc}?annotation=hmm"
TCDB_FASTA = "https://tcdb.org/public/tcdb"
UA = {"User-Agent": "raven-toolbox transporter-data build"}


def _get(url: str, timeout: int = 60) -> bytes:
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
        data = r.read()
    return gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data  # InterPro serves gzip


def build_pfam_hmm(out_dir: Path) -> Path:
    hmm_path = out_dir / "transporter_pfam.hmm"
    fetched, missing = 0, []
    with hmm_path.open("wb") as fh:
        for acc, (fam, _classes) in sorted(PFAM_TRANSPORTERS.items()):
            try:
                hmm = _get(INTERPRO_HMM.format(acc=acc))
            except Exception as exc:  # noqa: BLE001
                missing.append(f"{acc} ({fam}): fetch failed - {exc}")
                continue
            if not hmm.lstrip().startswith(b"HMMER3"):
                missing.append(f"{acc} ({fam}): no HMMER3 model returned")
                continue
            fh.write(hmm if hmm.endswith(b"\n") else hmm + b"\n")
            fetched += 1
            print(f"  + {acc} {fam}", flush=True)
    print(f"fetched {fetched}/{len(PFAM_TRANSPORTERS)} transporter HMMs -> {hmm_path}")
    if missing:
        print("MISSING/INVALID (fix the accession in transporter_tables.py):")
        for m in missing:
            print(f"  ! {m}")
    # No hmmpress: the backend uses `hmmsearch` (query = these HMMs, target = proteome), which reads a
    # plain multi-model .hmm directly; pressing is only needed for hmmscan's random access.
    return hmm_path


def build_tcdb_diamond(out_dir: Path) -> Path:
    fasta = out_dir / "tcdb.fasta"
    print(f"downloading TCDB sequences -> {fasta}", flush=True)
    fasta.write_bytes(_get(TCDB_FASTA, timeout=180))
    n = fasta.read_text(errors="ignore").count(">")
    print(f"  {n} sequences")
    diamond = resolve_binary("diamond")
    db = out_dir / "tcdb"
    subprocess.run([diamond, "makedb", "--in", str(fasta), "--db", str(db)], check=True)
    return db.with_suffix(".dmnd")


def _checksums(out_dir: Path) -> None:
    lines = []
    for p in sorted(out_dir.iterdir()):
        if p.is_file() and p.suffix not in {".fasta"}:
            lines.append(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}")
    (out_dir / "SHA256SUMS").write_text("\n".join(lines) + "\n")
    print(f"wrote {out_dir/'SHA256SUMS'} ({len(lines)} files)")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=Path("dist/transporters"))
    ap.add_argument("--only-pfam", action="store_true")
    ap.add_argument("--only-tcdb", action="store_true")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    if not args.only_tcdb:
        build_pfam_hmm(args.out)
    if not args.only_pfam:
        build_tcdb_diamond(args.out)
    _checksums(args.out)
    print("done.", file=sys.stderr)


if __name__ == "__main__":
    main()
