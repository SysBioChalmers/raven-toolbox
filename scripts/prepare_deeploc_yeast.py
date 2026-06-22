#!/usr/bin/env python3
"""Prepare a DeepLoc 2.1 input FASTA of protein sequences for yeast-GEM's genes.

DeepLoc 2.1 has no public batch API, so the workflow is: write this FASTA, run DeepLoc 2.1 on it
yourself (web server https://services.healthtech.dtu.dk/services/DeepLoc-2.1/, max 500 sequences per
submission — the FASTA is chunked accordingly; or the downloadable standalone), then load the
result with :func:`raven_toolbox.localization.load_deeploc` and feed it to ``predict_localization``.

Sequences are fetched from UniProtKB keyed by ordered-locus (ORF) name, which is exactly yeast-GEM's
gene id, so the FASTA headers — and therefore DeepLoc's ``Protein_ID`` output column — line up with
the model with no remapping.

Usage
-----
    python scripts/prepare_deeploc_yeast.py \\
        --yeast-gem ~/github/yeast-GEM/model/yeast-GEM.xml \\
        --out deeploc_yeast.fasta
"""
from __future__ import annotations

import argparse
from pathlib import Path

import cobra

from raven_toolbox.localization import prepare_deeploc_input


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--yeast-gem", required=True, type=Path,
                    help="path to the yeast-GEM SBML file")
    ap.add_argument("--out", default=Path("deeploc_yeast.fasta"), type=Path,
                    help="output FASTA path (chunked as <stem>_001.fasta, … past --max-per-file)")
    ap.add_argument("--organism", default=559292, type=int,
                    help="UniProt taxon id (default 559292, S. cerevisiae S288C)")
    ap.add_argument("--max-per-file", default=500, type=int,
                    help="sequences per FASTA file (500 = DeepLoc 2.1 web limit; 0 = no chunking)")
    ap.add_argument("--include-unreviewed", action="store_true",
                    help="also search UniProt TrEMBL (not just curated Swiss-Prot)")
    args = ap.parse_args()

    model = cobra.io.read_sbml_model(str(args.yeast_gem))
    print(f"{model.id or 'model'}: {len(model.genes)} genes", flush=True)

    res = prepare_deeploc_input(
        model, args.organism, args.out,
        max_records_per_file=(args.max_per_file or None),
        reviewed=not args.include_unreviewed,
    )

    print(res)
    print("wrote:")
    for p in res.paths:
        print(f"  {p}")
    if res.missing:
        print(f"{len(res.missing)} gene(s) had no reviewed UniProt sequence, e.g. "
              f"{', '.join(res.missing[:10])}"
              f"{' …' if len(res.missing) > 10 else ''}")
        print("  (try --include-unreviewed, or fetch those sequences from SGD/another source)")
    print("\nNext: run DeepLoc 2.1 on the FASTA, then "
          "load_deeploc(<output.csv>, compartment_map=DEFAULT_COMPARTMENT_MAP).")


if __name__ == "__main__":
    main()
