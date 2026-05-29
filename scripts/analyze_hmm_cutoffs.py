#!/usr/bin/env python3
"""Cut-off sensitivity for the KEGG HMM query path (step 3b.5).

Cross-validates ``assign_kos`` against an organism's *real* KEGG gene→KO
annotation (from the ``organism_gene_ko`` table) and sweeps the E-value cut-off
and the two score-ratio filters. Produces the tables in
``docs/kegg_hmm_cutoff_calibration.md``.

Usage
-----
    python scripts/analyze_hmm_cutoffs.py \
        --artefacts ~/keggdb_artefacts \
        --proteome /path/to/org.pep \
        --org sce --library ~/keggdb_artefacts/eukaryotes.hmm

``--proteome`` is the organism's protein FASTA (headers ``>org:gene ...``, e.g.
extracted from KEGG ``genes.pep``). ``--tblout`` may be given instead of
``--library`` to reuse a cached ``hmmscan --tblout`` file. Requires ``hmmscan``
on PATH or via ``RAVEN_PYTHON_HMMER`` when ``--library`` is used.

Caveat: organisms present in the library's training set give an upper bound on
recall; the comparison is relative (see the doc).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from raven_python.reconstruction.kegg.parse import read_kegg_table
from raven_python.reconstruction.kegg.query import (
    assign_kos,
    parse_hmmscan_tblout,
    run_hmmscan,
)

CUTOFFS = (1e-10, 1e-20, 1e-30, 1e-50, 1e-70, 1e-100)
KO_RATIOS = (0.0, 0.3, 0.5)
G_RATIOS = (0.5, 0.8, 0.95)


def load_ko2rxn(artefacts: Path) -> dict[str, set[str]]:
    tbl = read_kegg_table(artefacts / "ko_reaction.tsv.gz")
    ko2rxn: dict[str, set[str]] = {}
    for ko, rxn in zip(tbl["ko"], tbl["reaction"], strict=True):
        ko2rxn.setdefault(ko, set()).add(rxn)
    return ko2rxn


def ground_truth(artefacts: Path, org: str, ko2rxn) -> tuple[set, set]:
    ogk = read_kegg_table(artefacts / "organism_gene_ko.tsv.xz")
    rows = ogk[ogk["organism"].str.lower() == org]
    pairs = set(zip(rows["gene"], rows["ko"], strict=True))
    rxns = {r for _, ko in pairs for r in ko2rxn.get(ko, ())}
    return pairs, rxns


def predicted_pairs(hits: pd.DataFrame, **kw) -> set:
    out = set()
    for ko, genes in assign_kos(hits, **kw).items():
        for g in genes:
            out.add((g.split(":", 1)[1] if ":" in g else g, ko))
    return out


def prf(pred: set, truth: set) -> tuple[float, float, float]:
    tp = len(pred & truth)
    rec = tp / len(truth) if truth else 0.0
    prec = tp / len(pred) if pred else 0.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    return prec, rec, f1


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--artefacts", type=Path, required=True)
    ap.add_argument("--org", required=True, help="KEGG organism code, e.g. sce")
    ap.add_argument("--proteome", type=Path, help="protein FASTA (headers >org:gene)")
    ap.add_argument("--library", type=Path, help="pressed HMM library for hmmscan")
    ap.add_argument("--tblout", type=Path, help="cached hmmscan --tblout (skip hmmscan)")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args(argv)

    if args.tblout:
        text = args.tblout.read_text()
    elif args.library and args.proteome:
        text = run_hmmscan(args.proteome, args.library, threads=args.threads)
    else:
        ap.error("give --tblout, or --library and --proteome")

    org = args.org.lower()
    hits = parse_hmmscan_tblout(text)
    hits = hits[hits["gene"].str.startswith(f"{org}:")].reset_index(drop=True)
    ko2rxn = load_ko2rxn(args.artefacts)
    gt_pairs, gt_rxns = ground_truth(args.artefacts, org, ko2rxn)

    print(f"\n{'='*70}\n{org}: {hits['gene'].nunique()} query genes with hits, "
          f"{len(gt_pairs)} true gene->KO pairs, {len(gt_rxns)} true reactions\n{'='*70}")

    best: dict[tuple, float] = {}
    for ko, gene, e in zip(hits["ko"], hits["gene"], hits["evalue"], strict=True):
        key = (gene.split(":", 1)[1], ko)
        if key not in best or e < best[key]:
            best[key] = e
    matched = np.array([e for k, e in best.items() if k in gt_pairs])
    novel = np.array([e for k, e in best.items() if k not in gt_pairs])

    def logq(arr, q):
        if not len(arr):
            return float("nan")
        v = np.quantile(arr, q)
        return np.log10(v) if v > 0 else -300.0

    print("\nlog10(E-value) percentiles  [matched=in annotation, novel=not]:")
    print(f"  {'group':<8}{'n':>7}{'p50':>8}{'p90':>8}{'p95':>8}{'p99':>8}")
    for name, arr in (("matched", matched), ("novel", novel)):
        print(f"  {name:<8}{len(arr):>7}{logq(arr,.5):>8.0f}{logq(arr,.9):>8.0f}"
              f"{logq(arr,.95):>8.0f}{logq(arr,.99):>8.0f}")

    print("\ncutoff sweep (min_score_ratio_ko=0.3, min_score_ratio_g=0.8):")
    print(f"  {'cutoff':>8}{'gKO_prec':>9}{'gKO_rec':>8}{'gKO_F1':>8}{'rxn_rec':>9}{'rxn_novel':>10}")
    for cutoff in CUTOFFS:
        pred = predicted_pairs(hits, cutoff=cutoff)
        prec, rec, f1 = prf(pred, gt_pairs)
        pred_rxns = {r for _, ko in pred for r in ko2rxn.get(ko, ())}
        rrec = len(pred_rxns & gt_rxns) / len(gt_rxns) if gt_rxns else 0.0
        print(f"  {cutoff:>8.0e}{prec:>9.2f}{rec:>8.2f}{f1:>8.2f}{rrec:>9.2f}"
              f"{len(pred_rxns - gt_rxns):>10}")

    print("\nratio sweep (cutoff=1e-50):")
    print(f"  {'ko_ratio':>9}{'g_ratio':>8}{'gKO_prec':>9}{'gKO_rec':>8}{'gKO_F1':>8}")
    for rko in KO_RATIOS:
        for rg in G_RATIOS:
            pred = predicted_pairs(hits, cutoff=1e-50,
                                   min_score_ratio_ko=rko, min_score_ratio_g=rg)
            prec, rec, f1 = prf(pred, gt_pairs)
            print(f"  {rko:>9.1f}{rg:>8.2f}{prec:>9.2f}{rec:>8.2f}{f1:>8.2f}")


if __name__ == "__main__":
    main()
