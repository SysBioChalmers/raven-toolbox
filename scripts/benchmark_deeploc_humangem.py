#!/usr/bin/env python3
"""Benchmark DeepLoc 2.1 against Human-GEM's gene localisations (gene-level, circularity-aware).

Human-GEM is a **lenient positive control**: human is the core of DeepLoc's training, so high
agreement is expected and does not by itself prove generalisation (pair it with the stringent plant
model -- see the AraCore benchmark). The benchmark is *gene-level*: Human-GEM ships a per-gene
localisation in ``model/genes.tsv`` (the ``compartments`` column, DeepLoc-style labels, possibly
several per gene), with a provenance column ``compDataSource``.

**Circularity guard.** 439/2848 of those gene compartments were assigned **by DeepLoc 2 itself**
(``compDataSource == "DeepLoc2"``). Scoring DeepLoc against them grades DeepLoc on its own output, so
they are dropped; only the SwissProt / CellAtlas-sourced genes are scored. The script reports the
*with-DeepLoc2* number too, to show how much the circular rows inflate it.

A gene's annotation is a set (e.g. ``Nucleus;Cytosol``); DeepLoc predicts the single dominant
location, so a call counts as correct when the predicted compartment is **in** the annotated set.
ASCII-only output (Windows console is cp1252).

Usage:
    python scripts/benchmark_deeploc_humangem.py --genes <path>/Human-GEM/model/genes.tsv \
        --csv data/deeploc/humangem/Human-GEM_deeploc_*.csv --doc /tmp/deeploc_humangem.md
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd

#: DeepLoc label and Human-GEM genes.tsv label -> a common compartment id. Human has every DeepLoc
#: compartment; "Inner mitochondria" folds into mitochondrion.
CANON = {
    # DeepLoc organelle labels
    "Cytoplasm": "c", "Nucleus": "n", "Extracellular": "e", "Cell membrane": "ce",
    "Mitochondrion": "m", "Endoplasmic reticulum": "er", "Lysosome/Vacuole": "ly",
    "Golgi apparatus": "g", "Peroxisome": "p",
    # Human-GEM genes.tsv labels not already covered
    "Cytosol": "c", "Mitochondria": "m", "Inner mitochondria": "m", "Lysosome": "ly",
}
DEEPLOC_ORG = ["Cytoplasm", "Nucleus", "Extracellular", "Cell membrane", "Mitochondrion",
               "Endoplasmic reticulum", "Lysosome/Vacuole", "Golgi apparatus", "Peroxisome"]
COMP_ORDER = ["c", "n", "er", "g", "m", "ce", "ly", "e", "p"]


def deeploc_top(csvs: list[Path]) -> tuple[pd.Series, pd.Series]:
    df = pd.concat([pd.read_csv(c) for c in csvs]).set_index("Protein_ID")
    cols = [c for c in DEEPLOC_ORG if c in df.columns]
    return df[cols].astype(float).idxmax(axis=1).map(CANON), df[cols].astype(float).max(axis=1)


def truth_sets(genes: pd.DataFrame) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for gid, comp in zip(genes["genes"], genes["compartments"].fillna(""), strict=True):
        s = {CANON.get(p.strip()) for p in str(comp).split(";") if p.strip()}
        s.discard(None)
        if s:
            out[gid] = s
    return out


def agreement(top: pd.Series, truth: dict[str, set[str]],
              addressable: set[str] | None = None) -> tuple[int, int]:
    """Agreement = DeepLoc top compartment is in the gene's annotated set. ``addressable`` restricts
    to genes whose predicted compartment exists in the truth vocabulary at all (else the prediction
    is structurally uncorroboratable, not wrong)."""
    ok = tot = 0
    for gid, t in truth.items():
        if gid in top.index and pd.notna(top[gid]):
            if addressable is not None and top[gid] not in addressable:
                continue
            tot += 1
            ok += top[gid] in t
    return ok, tot


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", required=True, type=Path, help="Human-GEM model/genes.tsv")
    ap.add_argument("--csv", required=True, nargs="+", type=Path)
    ap.add_argument("--doc", type=Path)
    args = ap.parse_args()

    genes = pd.read_csv(args.genes, sep="\t")
    top, conf = deeploc_top(args.csv)
    src = genes["compDataSource"].fillna("")
    indep = genes[~src.str.contains("DeepLoc2") & genes["compartments"].notna()]
    circ = genes[src.str.contains("DeepLoc2") & genes["compartments"].notna()]

    truth_all = truth_sets(genes[genes["compartments"].notna()])
    truth_indep = truth_sets(indep)
    truth_circ = truth_sets(circ)

    # Compartments the independent truth can actually express (DeepLoc may predict more).
    addressable = {c for s in truth_indep.values() for c in s}
    ok_i, n_i = agreement(top, truth_indep)
    ok_ia, n_ia = agreement(top, truth_indep, addressable=addressable)
    ok_a, n_a = agreement(top, truth_all)
    ok_c, n_c = agreement(top, truth_circ)
    # labels DeepLoc predicts that the independent truth never contains
    off_vocab = sorted({top[g] for g in truth_indep if g in top.index and pd.notna(top[g])}
                       - addressable)

    out = ["# DeepLoc 2.1 vs Human-GEM gene localisations", "",
           f"Gene-level benchmark against `genes.tsv` `compartments` (a gene's annotated location "
           f"set). {len(genes)} genes; {len(top)} DeepLoc predictions. **Lenient positive control** "
           "-- human is core DeepLoc training, so high agreement is expected (pair with the "
           "stringent [AraCore benchmark](deeploc_aracore_benchmark.md)).", "",
           "## Circularity guard", "",
           f"Human-GEM assigned {len(circ)} gene compartments with DeepLoc 2 itself "
           f"(`compDataSource == DeepLoc2`); scoring against them is circular, so they are dropped. "
           f"Only the {len(indep)} SwissProt/CellAtlas-sourced genes are scored.", "",
           "| gene set | n scored | agreement (DeepLoc top in annotated set) |",
           "|---|--:|--:|",
           f"| **independent** (SwissProt/CellAtlas) | {n_i} | **{ok_i / n_i:.1%}** |",
           f"| independent, addressable only | {n_ia} | **{ok_ia / n_ia:.1%}** |",
           f"| circular (DeepLoc2-sourced) | {n_c} | {ok_c / n_c:.1%} |",
           f"| all (incl. circular) | {n_a} | {ok_a / n_a:.1%} |", "",
           f"The DeepLoc2-sourced rows score {ok_c / n_c:.1%} (DeepLoc grading itself) and inflate "
           f"the naive all-genes number by {(ok_a / n_a - ok_i / n_i) * 100:+.1f} pp -- which is why "
           "they must be excluded.", "",
           f"**Addressability.** Dropping the circular rows also removes the entire `"
           f"{'`/`'.join(off_vocab)}` truth vocabulary: in Human-GEM, *every* gene annotated to those "
           f"compartments was DeepLoc2-sourced (the independent sources annotate no metabolic gene "
           f"there). DeepLoc still routes {n_i - n_ia} independent genes to them, which cannot be "
           f"corroborated and count as misses. On the {len(addressable)} compartments the independent "
           f"truth can express, agreement is **{ok_ia / n_ia:.1%}** -- the fair DeepLoc-skill number; "
           f"the 75% headline is the conservative floor.", "",
           "## Per-compartment (independent genes)", "",
           "Precision by predicted label: among genes DeepLoc calls a compartment, how often that "
           "compartment is in the gene's annotated set.", "",
           "| predicted | n | in annotated set | precision |", "|---|--:|--:|--:|"]
    per: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for gid, t in truth_indep.items():
        if gid in top.index and pd.notna(top[gid]):
            per[top[gid]][0] += 1
            per[top[gid]][1] += top[gid] in t
    for c in COMP_ORDER + [c for c in per if c not in COMP_ORDER]:
        if c in per:
            n, ok = per[c]
            out.append(f"| {c} | {n} | {ok} | {ok / n:.1%} |")
    out.append("")

    out += ["## Accuracy by DeepLoc confidence (independent genes)", "",
            "| confidence bin | n | correct | accuracy |", "|---|--:|--:|--:|"]
    bins = [(0.0, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.01)]
    rows = {b: [0, 0] for b in bins}
    for gid, t in truth_indep.items():
        if gid in top.index and pd.notna(top[gid]):
            cf = float(conf[gid])
            for b in bins:
                if b[0] <= cf < b[1]:
                    rows[b][0] += 1
                    rows[b][1] += top[gid] in t
                    break
    for b in bins:
        n, ok = rows[b]
        if n:
            out.append(f"| [{b[0]:.1f}, {min(b[1], 1.0):.1f}] | {n} | {ok} | {ok / n:.1%} |")
    out.append("")

    text = "\n".join(out) + "\n"
    print(text)
    if args.doc:
        args.doc.write_text(text, encoding="utf-8", newline="")
        print(f"wrote {args.doc}")


if __name__ == "__main__":
    main()
