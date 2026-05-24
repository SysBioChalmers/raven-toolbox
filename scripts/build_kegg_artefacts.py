#!/usr/bin/env python
"""Build the publishable KEGG artefact set for one release (maintainer-side).

Runs the maintainer pipeline against an arranged KEGG dump (see
``download_kegg_dump`` / ``fetch_keggdb``):

* 3b.2 — ``parse_kegg_dump`` → ``reference_model.yml.gz`` + the gzipped-TSV tables;
* 3b.3 — ``build_hmm_library`` per domain → a pressed ``<domain>.hmm`` (+ hmmpress
  sidecars), named so :func:`ravengem.data.ensure_kegg_hmm_library` can fetch them.

Everything lands in ``--out`` ready to upload as release assets; feed that
directory to ``scripts/make_registry_snippet.py data`` to emit the registry entry.

Examples
--------
Tables + reference model only (fast, no binaries)::

    python scripts/build_kegg_artefacts.py --keggdb keggdb --out artefacts

Full build incl. both HMM libraries (slow; needs HMMER/MAFFT/CD-HIT)::

    python scripts/build_kegg_artefacts.py --keggdb keggdb --out artefacts \\
        --hmms --threads 8
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from ravengem.reconstruction.kegg import (
    build_hmm_library,
    parse_kegg_dump,
    read_kegg_table,
)

# hmmpress sidecar extensions, alongside the .hmm.
_HMM_SIDECARS = (".h3f", ".h3i", ".h3m", ".h3p")


def _publish_library(work: dict, out_dir: Path, domain: str) -> Path:
    """Copy a built ``library.hmm`` (+ sidecars) to ``out_dir/<domain>.hmm``."""
    library = work["library"]
    if library is None:
        raise SystemExit(f"No HMMs built for {domain!r}; nothing to publish.")
    target = out_dir / f"{domain}.hmm"
    shutil.copyfile(library, target)
    for suffix in _HMM_SIDECARS:
        sidecar = library.with_name(library.name + suffix)
        if sidecar.exists():
            shutil.copyfile(sidecar, target.with_name(target.name + suffix))
    return target


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--keggdb", required=True, type=Path, help="arranged KEGG dump directory")
    parser.add_argument("--out", required=True, type=Path, help="artefact output directory")
    parser.add_argument("--hmms", action="store_true", help="also build the HMM libraries")
    parser.add_argument(
        "--domains", nargs="+", default=["prokaryotes", "eukaryotes"], help="HMM domains to build"
    )
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--seq-identity", type=float, default=0.9, help="CD-HIT identity (-1 skips)")
    parser.add_argument(
        "--max-sequences", type=int, default=None,
        help="cap sequences per KO after CD-HIT (bounds MAFFT time/memory on huge KOs)",
    )
    args = parser.parse_args(argv)

    args.out.mkdir(parents=True, exist_ok=True)
    print(">>> Parsing KEGG dump (3b.2)...")
    paths = parse_kegg_dump(args.keggdb, args.out)
    for name, path in paths.items():
        print(f"    {name}: {path}")

    if args.hmms:
        ogk = read_kegg_table(args.out / "organism_gene_ko.tsv.gz")
        genes_pep = args.keggdb / "genes.pep"
        taxonomy = args.keggdb / "taxonomy"
        for domain in args.domains:
            print(f">>> Building HMM library for {domain} (3b.3)...")
            work = build_hmm_library(
                ogk, genes_pep, taxonomy, args.out / f"_hmms-{domain}",
                domain=domain, seq_identity=args.seq_identity,
                max_sequences=args.max_sequences, threads=args.threads,
            )
            published = _publish_library(work, args.out, domain)
            print(f"    {domain}: {published} ({len(work['hmms'])} profiles)")

    print(f"\n>>> Done. Upload the contents of {args.out} as release assets, then run:")
    print("    python scripts/make_registry_snippet.py data --dataset kegg "
          f"--version <VER> --dir {args.out} --base-url <RELEASE_URL>")


if __name__ == "__main__":
    main()
