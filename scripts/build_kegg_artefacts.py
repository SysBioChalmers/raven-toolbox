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

Resuming: the build is idempotent. If it fails partway, re-run the **same
command** — each stage is skipped when its output already exists (parsed tables,
taxonomy, per-domain HMM library, core bundle), and the per-KO HMM build picks up
where it left off. Pass ``--force`` to rebuild everything from scratch.
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


# The parse outputs later bundled into <prefix>core.tar.gz; their presence (all of
# them) marks the parse step as already done, for the resumable re-run path.
_CORE_NAMES = ("reference_model", "ko_reaction", "ko_names", "organism_gene_ko", "rxn_flags")


def _core_paths(out_dir: Path, prefix: str) -> dict[str, Path]:
    """Expected loose parse-output paths (keys match ``parse_kegg_dump``'s return)."""
    names = {
        "reference_model": f"{prefix}reference_model.yml.gz",
        "ko_reaction": f"{prefix}ko_reaction.tsv.gz",
        "ko_names": f"{prefix}ko_names.tsv.gz",
        "organism_gene_ko": f"{prefix}organism_gene_ko.tsv.gz",
        "rxn_flags": f"{prefix}rxn_flags.tsv.gz",
    }
    return {key: out_dir / name for key, name in names.items()}


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
        "--force", action="store_true",
        help="rebuild everything, ignoring existing outputs. By default each stage "
             "is skipped when its output already exists, so a re-run after a failure "
             "continues instead of starting over.",
    )
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
    force = args.force

    # Which requested HMM libraries are not yet published? (Their absence tells us
    # whether the parse outputs are still needed as build input.)
    domains_to_build = (
        [d for d in args.domains if force or not (args.out / f"{prefix}{d}.hmm.gz").exists()]
        if args.hmms
        else []
    )

    core = _core_paths(args.out, prefix)
    bundle = args.out / f"{prefix}core.tar.gz"
    parse_done = all(p.exists() for p in core.values())

    # 3b.2 — parse. Idempotent: skip when the loose outputs already exist. Re-parse
    # only if forcing, or if they're missing and still needed (as HMM build input,
    # or because the core bundle hasn't been written yet).
    if force or (not parse_done and (domains_to_build or not bundle.exists())):
        print(">>> Parsing KEGG dump (3b.2)...")
        paths = parse_kegg_dump(args.keggdb, args.out, version=args.version, progress=True)
        for name, path in paths.items():
            print(f"    {name}: {path}")
    else:
        paths = core
        print(
            ">>> Parse outputs present; skipping parse (use --force to redo)."
            if parse_done
            else ">>> Core bundle present and nothing to rebuild; skipping parse."
        )

    # Publish the taxonomy file too: domain split, plus the source for phyl_dist
    # (RAVEN's keggPhylDist, used by GECKO). It is a raw dump file, so gzip it as-is.
    tax_src = args.keggdb / "taxonomy"
    tax_out = args.out / f"{prefix}taxonomy.gz"
    if tax_src.is_file() and (force or not tax_out.exists()):
        with open(tax_src, "rb") as src, gzip.open(tax_out, "wb") as out:
            shutil.copyfileobj(src, out)
        print(f"    taxonomy: {tax_out}")
    elif tax_out.exists():
        print(f"    taxonomy: {tax_out} (exists; skipped)")

    if args.hmms:
        genes_pep = args.keggdb / "genes.pep"
        taxonomy = args.keggdb / "taxonomy"
        ogk = None  # loaded lazily, only when a domain actually needs building
        for domain in args.domains:
            pub_path = args.out / f"{prefix}{domain}.hmm.gz"
            if pub_path.exists() and not force:
                print(f">>> HMM library for {domain} present; skipping (use --force to redo).")
                continue
            if ogk is None:
                ogk = read_kegg_table(paths["organism_gene_ko"])
            print(f">>> Building HMM library for {domain} (3b.3)...")
            work = build_hmm_library(
                ogk, genes_pep, taxonomy, args.out / f"_hmms-{domain}",
                domain=domain, seq_identity=args.seq_identity,
                parttree_residues=args.parttree_residues, threads=args.threads,
                progress=True, force=force,
            )
            published = _publish_library(work, args.out, domain, prefix)
            print(f"    {domain}: {published} ({len(work['hmms'])} profiles)")

    # Bundle the core model artefacts (reference model + tables) into one archive that
    # ensure_kegg_data fetches and extracts. HMMs and taxonomy stay separate assets.
    # (After the HMM step, which reads organism_gene_ko.) Idempotent: skip if the
    # bundle already exists; the loose members are removed once bundled.
    if force or not bundle.exists():
        members = [paths[n] for n in _CORE_NAMES]
        if all(m.exists() for m in members):
            with tarfile.open(bundle, "w:gz") as tar:
                for member in members:
                    tar.add(member, arcname=member.name)
            for member in members:
                member.unlink()
            print(f"    core bundle: {bundle} ({len(members)} files)")
        else:
            missing = [m.name for m in members if not m.exists()]
            print(f"    core bundle: skipped (loose members missing: {missing})")
    else:
        print(f"    core bundle: {bundle} (exists; skipped)")

    print(f"\n>>> Done. Upload the contents of {args.out} as release assets, then run:")
    print("    python scripts/make_registry_snippet.py data --dataset kegg "
          f"--version {args.version or '<VER>'} --dir {args.out} --tag <RELEASE_TAG>")


if __name__ == "__main__":
    main()
