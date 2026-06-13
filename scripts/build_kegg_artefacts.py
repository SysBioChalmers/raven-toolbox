#!/usr/bin/env python
"""Build the publishable KEGG artefact set for one release (maintainer-side).

Runs the maintainer pipeline against an arranged KEGG dump (see
``download_kegg_dump`` / ``fetch_keggdb``):

* 3b.2 — ``parse_kegg_dump`` → ``reference_model.yml.gz`` + the gzipped-TSV tables;
* 3b.3 — ``build_hmm_library`` per domain → a gzipped concatenated flatfile
  ``<version>_<domain>.hmm.gz`` (the client decompresses it and searches with
  ``hmmsearch``), named so :func:`raven_toolbox.data.ensure_kegg_hmm_library` can fetch it.

Pass ``--version`` (e.g. ``kegg116``) to version-prefix every output filename, matching
the published release assets. Everything lands in ``--out`` ready to upload; feed that
directory to ``scripts/make_registry_snippet.py data`` to emit the registry entry.

Examples
--------
Tables + reference model only (fast, no binaries)::

    python scripts/build_kegg_artefacts.py --keggdb keggdb --out artefacts --version kegg116

Full build incl. both HMM libraries (slow; needs HMMER/MAFFT/CD-HIT)::

    python scripts/build_kegg_artefacts.py --keggdb keggdb --out artefacts \\
        --version kegg116 --hmms --threads 8
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import tarfile
from pathlib import Path

from raven_toolbox.reconstruction.kegg import (
    build_hmm_library,
    parse_kegg_dump,
    read_kegg_table,
)


def _publish_library(work: dict, out_dir: Path, domain: str, prefix: str = "") -> Path:
    """Gzip a built ``library.hmm`` to ``out_dir/<prefix><domain>.hmm.gz``.

    Only the concatenated flatfile is published (gzip, portable across HMMER
    versions); the client decompresses it and searches with ``hmmsearch``, so the
    same artefact also serves MATLAB RAVEN.
    """
    library = work["library"]
    if library is None:
        raise SystemExit(f"No HMMs built for {domain!r}; nothing to publish.")
    target = out_dir / f"{prefix}{domain}.hmm.gz"
    with open(library, "rb") as src, gzip.open(target, "wb") as out:
        shutil.copyfileobj(src, out)
    return target


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--keggdb", required=True, type=Path, help="arranged KEGG dump directory")
    parser.add_argument("--out", required=True, type=Path, help="artefact output directory")
    parser.add_argument(
        "--version", default=None,
        help="KEGG version (e.g. kegg116); version-prefixes every output filename",
    )
    parser.add_argument("--hmms", action="store_true", help="also build the HMM libraries")
    parser.add_argument(
        "--domains", nargs="+", default=["prokaryotes", "eukaryotes"], help="HMM domains to build"
    )
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--seq-identity", type=float, default=0.9, help="CD-HIT identity (-1 skips)")
    parser.add_argument(
        "--parttree-residues", type=int, default=None,
        help="total-residue budget above which MAFFT uses PartTree (default 1M, tuned "
             "for ~7 GB RAM; raise on machines with more memory)",
    )
    args = parser.parse_args(argv)

    args.out.mkdir(parents=True, exist_ok=True)
    prefix = f"{args.version}_" if args.version else ""
    print(">>> Parsing KEGG dump (3b.2)...")
    paths = parse_kegg_dump(args.keggdb, args.out, version=args.version)
    for name, path in paths.items():
        print(f"    {name}: {path}")

    # Publish the taxonomy file too: domain split, plus the source for phyl_dist
    # (RAVEN's keggPhylDist, used by GECKO). It is a raw dump file, so gzip it as-is.
    tax_src = args.keggdb / "taxonomy"
    if tax_src.is_file():
        tax_out = args.out / f"{prefix}taxonomy.gz"
        with open(tax_src, "rb") as src, gzip.open(tax_out, "wb") as out:
            shutil.copyfileobj(src, out)
        print(f"    taxonomy: {tax_out}")

    if args.hmms:
        ogk = read_kegg_table(paths["organism_gene_ko"])
        genes_pep = args.keggdb / "genes.pep"
        taxonomy = args.keggdb / "taxonomy"
        for domain in args.domains:
            print(f">>> Building HMM library for {domain} (3b.3)...")
            work = build_hmm_library(
                ogk, genes_pep, taxonomy, args.out / f"_hmms-{domain}",
                domain=domain, seq_identity=args.seq_identity,
                parttree_residues=args.parttree_residues, threads=args.threads,
            )
            published = _publish_library(work, args.out, domain, prefix)
            print(f"    {domain}: {published} ({len(work['hmms'])} profiles)")

    # Bundle the core model artefacts (reference model + tables) into one archive that
    # ensure_kegg_data fetches and extracts. HMMs and taxonomy stay separate assets.
    # (After the HMM step, which reads organism_gene_ko.)
    core_members = [
        paths[n]
        for n in ("reference_model", "ko_reaction", "ko_names", "organism_gene_ko", "rxn_flags")
    ]
    bundle = args.out / f"{prefix}core.tar.gz"
    with tarfile.open(bundle, "w:gz") as tar:
        for member in core_members:
            tar.add(member, arcname=member.name)
    for member in core_members:
        member.unlink()
    print(f"    core bundle: {bundle} ({len(core_members)} files)")

    print(f"\n>>> Done. Upload the contents of {args.out} as release assets, then run:")
    print("    python scripts/make_registry_snippet.py data --dataset kegg "
          f"--version {args.version or '<VER>'} --dir {args.out} --base-url <RELEASE_URL>")


if __name__ == "__main__":
    main()
