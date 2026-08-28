#!/usr/bin/env python3
"""Cross-check DeepLoc 2.1 against UniProt and yeast-GEM (gene-level, organelle-collapsed).

yeast-GEM's curated localisation is a fair indication, not indisputable ground truth, so this
triangulates with UniProt's curated `Subcellular location` annotation: where DeepLoc disagrees with
yeast-GEM but agrees with UniProt, yeast-GEM is the likely outlier. Also reports how well DeepLoc's
own confidence predicts being corroborated - the basis for `load_deeploc(min_confidence=...)`.

Needs network (UniProt REST).

Usage:
    python scripts/compare_localization_sources.py --yeast-gem <path>/yeast-GEM.xml \
        --csv data/deeploc/yeast-GEM_deeploc_*.csv --organism 559292
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import cobra
import numpy as np
import pandas as pd

from raven_toolbox.localization import fetch_uniprot_localization

ORG = {"Cytoplasm": "c", "Nucleus": "n", "Extracellular": "e", "Cell membrane": "ce",
       "Peroxisome": "p", "Mitochondrion": "m", "Endoplasmic reticulum": "er",
       "Golgi apparatus": "g", "Lysosome/Vacuole": "v"}
COLLAPSE = {"erm": "er", "mm": "m", "gm": "g", "vm": "v"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--yeast-gem", required=True, type=Path)
    ap.add_argument("--csv", required=True, nargs="+", type=Path)
    ap.add_argument("--organism", default=559292, type=int)
    args = ap.parse_args()

    # DeepLoc: gene -> top organelle (raw probabilities, before normalisation, for confidence)
    raw = pd.concat([pd.read_csv(c) for c in args.csv]).set_index("Protein_ID")
    org_cols = [c for c in ORG if c in raw.columns]
    dl_top = {g: ORG[raw.loc[g, org_cols].astype(float).idxmax()] for g in raw.index}
    dl_conf = raw[org_cols].astype(float).max(axis=1)

    model = cobra.io.read_sbml_model(str(args.yeast_gem))
    ygem: dict[str, set[str]] = defaultdict(set)
    for r in model.reactions:
        if r.boundary:
            continue
        comps = {COLLAPSE.get(m.compartment, m.compartment) for m in r.metabolites if m.compartment}
        if len(comps) == 1:
            for gene in r.genes:
                ygem[gene.id].add(next(iter(comps)))

    up_df = fetch_uniprot_localization(args.organism).df
    up = {g: {c for c in up_df.columns if up_df.loc[g, c] > 0} for g in up_df.index}
    up = {g: s for g, s in up.items() if s}

    print(f"coverage: yeast-GEM {len(ygem)} genes, DeepLoc {len(dl_top)}, "
          f"UniProt-with-location {len(up)}\n")

    def rate(pairs):
        ok = sum(pairs)
        return f"{ok}/{len(pairs)} = {ok / len(pairs):.1%}" if pairs else "n/a"

    print("Pairwise agreement (gene-level, organelle-collapsed):")
    print(f"  DeepLoc  -> yeast-GEM   {rate([dl_top[g] in ygem[g] for g in dl_top if g in ygem])}")
    print(f"  DeepLoc  -> UniProt     {rate([dl_top[g] in up[g] for g in dl_top if g in up])}")
    print(f"  UniProt <-> yeast-GEM   "
          f"{rate([bool(up[g] & ygem[g]) for g in up if g in ygem])}  (two curated sources)")

    common = [g for g in dl_top if g in ygem and g in up]
    cat = defaultdict(int)
    for g in common:
        cat[(dl_top[g] in ygem[g], dl_top[g] in up[g])] += 1
    print(f"\nThree-way on {len(common)} genes annotated by all sources (DeepLoc top vs each):")
    print(f"  agrees with BOTH          {cat[(True, True)]:4} ({cat[(True, True)] / len(common):.1%})")
    print(f"  UniProt only (yGEM outlier){cat[(False, True)]:4} ({cat[(False, True)] / len(common):.1%})"
          "  <- DeepLoc likely right")
    print(f"  yeast-GEM only            {cat[(True, False)]:4} ({cat[(True, False)] / len(common):.1%})")
    print(f"  agrees NEITHER            {cat[(False, False)]:4} ({cat[(False, False)] / len(common):.1%})"
          "  <- DeepLoc likely wrong")

    anyref = [g for g in dl_top if g in ygem or g in up]
    conf = np.array([dl_conf[g] for g in anyref])
    corr = np.array([dl_top[g] in (ygem.get(g, set()) | up.get(g, set())) for g in anyref])
    print(f"\nDeepLoc confidence vs corroboration (top organelle in yeast-GEM OR UniProt), "
          f"n={len(anyref)}:")
    for lo, hi in [(0.0, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.01)]:
        mask = (conf >= lo) & (conf < hi)
        if mask.sum():
            print(f"  confidence [{lo:.2f},{hi:.2f}): n={int(mask.sum()):4}  "
                  f"corroborated={corr[mask].mean():.1%}")


if __name__ == "__main__":
    main()
