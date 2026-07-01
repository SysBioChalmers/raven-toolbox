#!/usr/bin/env python3
"""Finetune the DeepLoc-loading hyperparameters against yeast-GEM (slow / ProtT5 data).

The hyperparameters the DeepLoc-vs-yeast-GEM data determines directly:

* ``load_deeploc(membrane_threshold=...)`` — the ``1 - P(Soluble)`` cut that routes
  ``Mitochondrion -> mm``. Tuned to the ROC-optimal point (Youden's J) on curated matrix-vs-membrane.
* ``load_deeploc(min_confidence=...)`` — the top-organelle-probability gate. A coverage/accuracy
  trade-off; reported against yeast-GEM and (if ``--organism`` and network) against the consensus
  (yeast-GEM **or** UniProt), the study's overfitting guard — never tune to yeast-GEM alone.
* triage ``DEEPLOC_COMPARTMENT_TRUST`` — per-compartment reliability, refreshed from the slow run's
  organelle-collapsed accuracy.

Reuses the benchmark helpers in ``benchmark_deeploc.py``. ASCII-only output (Windows console cp1252).

Usage:
    python scripts/finetune_localization_yeast.py --model <path>/yeast-GEM.xml \
        --csv data/deeploc/yeast-GEM_deeploc_*.csv [--organism 559292]
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import cobra
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_deeploc import (  # noqa: E402
    SPECIES,
    DeepLoc,
    build_truth,
    collapse_membranes,
    organelle_accuracy,
)


def membrane_threshold_sweep(truth, rxn_genes, dl: DeepLoc, lum: str, mem: str) -> tuple[list[str], float]:
    """Sweep the membrane cut for an organelle's lumen/membrane split; ROC-optimal by Youden's J.

    Restricted to reactions DeepLoc placed in the correct organelle (the rule only fires there).
    """
    sig, lab = [], []
    for rid, c in truth.items():
        if c not in (lum, mem):
            continue
        genes = [g for g in rxn_genes[rid] if dl.covers(g)]
        if not genes or dl.rxn_organelle(genes) != lum:
            continue
        sig.append(dl.rxn_membrane(genes))
        lab.append(1 if c == mem else 0)
    sig_a, lab_a = np.array(sig), np.array(lab)
    n_mem, n_lum = int(lab_a.sum()), int((lab_a == 0).sum())
    lines = [f"### Mitochondrion matrix vs membrane ({lum}/{mem}) -- routing threshold",
             "", f"Among {len(sig)} reactions DeepLoc placed in `{lum}` "
             f"({n_lum} curated `{lum}`, {n_mem} curated `{mem}`): predict `{mem}` when "
             "`1 - P(Soluble) >= threshold`. Balanced accuracy = mean of the two recalls.", "",
             "| threshold | mm recall | m recall | balanced acc | Youden J |",
             "|---|--:|--:|--:|--:|"]
    best_t, best_j = 0.5, -1.0
    for t in [round(x, 2) for x in np.arange(0.30, 0.91, 0.05)]:
        pred_mem = sig_a >= t
        mm_recall = pred_mem[lab_a == 1].mean() if n_mem else float("nan")
        m_recall = (~pred_mem[lab_a == 0]).mean() if n_lum else float("nan")
        bal = np.nanmean([mm_recall, m_recall])
        j = mm_recall + m_recall - 1
        mark = ""
        if j > best_j:
            best_j, best_t = j, t
        lines.append(f"| {t:.2f}{mark} | {mm_recall:.1%} | {m_recall:.1%} | {bal:.1%} | {j:.3f} |")
    # mark the winner row
    lines = [ln.replace(f"| {best_t:.2f} |", f"| **{best_t:.2f}** |") for ln in lines]
    lines.append("")
    lines.append(f"**Optimal `membrane_threshold` = {best_t:.2f}** (Youden's J = {best_j:.3f}); "
                 f"the current default is 0.50.")
    lines.append("")
    return lines, best_t


def min_confidence_sweep(truth, rxn_genes, dl: DeepLoc, consensus: dict | None) -> list[str]:
    """Coverage/accuracy of the confidence gate vs yeast-GEM and (optionally) the consensus."""
    # per-reaction: top-gene confidence, predicted organelle, curated compartment
    recs = []
    for rid, c in truth.items():
        genes = [g for g in rxn_genes[rid] if dl.covers(g)]
        if not genes:
            continue
        conf = max(dl.org_conf[g] for g in genes)
        recs.append((conf, dl.rxn_organelle(genes) == c, rid))
    recs.sort(reverse=True)
    total = len(recs)
    lines = ["### `min_confidence` gate -- coverage vs accuracy (organelle-collapsed)", "",
             "Drop reactions whose top gene is below the threshold; accuracy is on what remains.", "",
             "| min_confidence | reactions kept | coverage | accuracy vs yeast-GEM |"]
    head = "|---|--:|--:|--:|"
    if consensus is not None:
        lines[-1] += " corroborated (yGEM or UniProt) |"
        head = "|---|--:|--:|--:|--:|"
    lines.append(head)
    for t in [0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        kept = [r for r in recs if r[0] >= t]
        if not kept:
            continue
        acc = sum(ok for _, ok, _ in kept) / len(kept)
        row = (f"| {t:.1f} | {len(kept)} | {len(kept) / total:.0%} | {acc:.1%} |")
        if consensus is not None:
            corr = np.mean([consensus.get(rid, False) for _, _, rid in kept])
            row += f" {corr:.1%} |"
        lines.append(row)
    lines.append("")
    return lines


def trust_table(model: cobra.Model, dl: DeepLoc, membrane_to_lumen: dict) -> list[str]:
    """Per-compartment organelle-collapsed accuracy -> refreshed DEEPLOC_COMPARTMENT_TRUST."""
    cmodel = collapse_membranes(model, membrane_to_lumen)
    ctruth, crxn_genes = build_truth(cmodel)
    per = organelle_accuracy(ctruth, crxn_genes, dl)
    trust: dict[str, float] = {}
    for c, (n, ok) in per.items():
        trust[c] = round(ok / n, 2)
    lines = ["### Refreshed `DEEPLOC_COMPARTMENT_TRUST` (slow, organelle-collapsed)", "",
             "Per-compartment reliability of a DeepLoc organelle call, from the slow run. Only `mm` "
             "inherits its lumen's trust (the mitochondrial split is the one validated routing, AUC "
             "~0.93); the other membranes (`erm`/`gm`/`vm`) and `ce`/`lp` stay 0 -- DeepLoc cannot "
             "reach them reliably.", "",
             "| compartment | n | accuracy (trust) |", "|---|--:|--:|"]
    for c in sorted(per, key=lambda c: -per[c][1] / per[c][0]):
        n, ok = per[c]
        lines.append(f"| {c} | {n} | {ok / n:.2f} |")
    lines.append("")
    # emit a copy-pasteable dict: organelle trusts (incl. ce/lp from the collapsed run) + mm
    # (validated mito split inherits m) + the non-validated membranes at 0.
    full = dict(trust)
    full["mm"] = trust.get("m", 0.0)
    for z in ("erm", "gm", "vm"):
        full[z] = 0.0
    full.setdefault("ce", 0.0)
    full.setdefault("lp", 0.0)
    ordered = sorted(full.items(), key=lambda kv: -kv[1])
    lines.append("```python")
    lines.append("DEEPLOC_COMPARTMENT_TRUST = {")
    lines.append("    " + ", ".join(f'"{c}": {v:.2f}' for c, v in ordered))
    lines.append("}")
    lines.append("```")
    lines.append("")
    return lines


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True, type=Path)
    ap.add_argument("--csv", required=True, nargs="+", type=Path)
    ap.add_argument("--organism", type=int, default=None,
                    help="UniProt organism id (e.g. 559292) to corroborate min_confidence; needs network")
    ap.add_argument("--doc", type=Path)
    args = ap.parse_args()

    cfg = SPECIES["yeast"]
    org = cfg["org"]
    membrane_to_lumen = {mem: lum for lum, mem in org.values() if mem}
    dl = DeepLoc(pd.concat([pd.read_csv(c) for c in args.csv]).set_index("Protein_ID"), org)
    model = cobra.io.read_sbml_model(str(args.model))
    truth, rxn_genes = build_truth(model)

    # Optional consensus (yeast-GEM OR UniProt) per reaction, for the min_confidence guard.
    consensus = None
    if args.organism:
        from raven_toolbox.localization import fetch_uniprot_localization
        up_df = fetch_uniprot_localization(args.organism).df
        collapse = {"erm": "er", "mm": "m", "gm": "g", "vm": "v"}
        up = {g: {c for c in up_df.columns if up_df.loc[g, c] > 0} for g in up_df.index}
        ygem: dict[str, set[str]] = defaultdict(set)
        for r in model.reactions:
            if r.boundary:
                continue
            comps = {collapse.get(m.compartment, m.compartment) for m in r.metabolites if m.compartment}
            if len(comps) == 1:
                for gene in r.genes:
                    ygem[gene.id].add(next(iter(comps)))
        consensus = {}
        for rid in truth:
            genes = [g for g in rxn_genes[rid] if dl.covers(g)]
            if not genes:
                continue
            pred = dl.rxn_organelle(genes)
            ref = set()
            for g in genes:
                ref |= ygem.get(g, set()) | up.get(g, set())
            consensus[rid] = pred in ref

    out = ["# Finetuning localisation hyperparameters on yeast-GEM (slow DeepLoc)", "",
           f"Model `{model.id or args.model.name}`, {len(dl.df)} slow-model DeepLoc predictions, "
           f"{len(truth)} single-compartment GPR reactions. Tunes the DeepLoc-loading "
           "hyperparameters against the curated compartmentalisation.", ""]
    out += ["## 1. Mitochondrial membrane routing (`membrane_threshold`)", ""]
    mem_lines, best_t = membrane_threshold_sweep(truth, rxn_genes, dl, "m", "mm")
    out += mem_lines
    out += ["## 2. Confidence gate (`min_confidence`)", ""]
    out += min_confidence_sweep(truth, rxn_genes, dl, consensus)
    out += ["## 3. Per-compartment trust (triage)", ""]
    out += trust_table(model, dl, membrane_to_lumen)

    text = "\n".join(out) + "\n"
    print(text)
    if args.doc:
        args.doc.write_text(text, encoding="utf-8", newline="")
        print(f"wrote {args.doc}")


if __name__ == "__main__":
    main()
