#!/usr/bin/env python3
"""Benchmark :func:`ravengem.localization.predict_localization` on yeast-GEM.

Treats yeast-GEM's curated compartmentalisation as ground truth, flattens the model with
:func:`merge_compartments` to a single compartment (so the algorithm cannot lean on
metabolite-topology evidence), then asks ``predict_localization`` to place every
GPR-annotated reaction back into a compartment given a per-gene score table.

The reference score table is derived directly from yeast-GEM (each gene scores 1.0 in
the compartments where its reactions actually live). Noise can be added — a configurable
fraction of genes have a random other compartment swapped in as the best score — to see
how the algorithm degrades with imperfect predictor evidence. With ``--scores-csv`` the
reference table is replaced by a real predictor output (WoLF PSORT / DeepLoc / hand-built
``gene_id × compartment`` CSV).

Outputs a per-noise-level accuracy summary and, optionally, a markdown table to a doc.

Usage
-----
    python scripts/benchmark_localization_yeast.py \\
        --yeast-gem ~/github/pcSecYeastSpecies/Model/yeastGEM.xml \\
        --noise 0,0.1,0.25,0.5 \\
        --doc /tmp/yeast_localization_benchmark.md
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import cobra
import numpy as np
import pandas as pd

from ravengem.localization import LocalizationScores, predict_localization
from ravengem.manipulation.compartments import merge_compartments

# --------------------------------------------------------------------------- inputs

def build_truth(model: cobra.Model) -> dict[str, str]:
    """For each single-compartment GPR-annotated reaction, ``{rxn_id: compartment}``.

    Boundary reactions and multi-compartment transports are excluded — those aren't
    placeable by the algorithm and shouldn't enter the benchmark.
    """
    truth: dict[str, str] = {}
    for r in model.reactions:
        if r.boundary or not r.genes:
            continue
        comps = {m.compartment for m in r.metabolites if m.compartment}
        if len(comps) != 1:
            continue
        truth[r.id] = next(iter(comps))
    return truth


def derive_scores_from_model(model: cobra.Model) -> LocalizationScores:
    """Each gene scores 1.0 in every compartment where its reactions actually live.

    For genes shared across compartments (dual-localised in the curation), all of those
    compartments get the top score — which is exactly the situation
    ``multi_compartment_penalty`` is designed to handle.
    """
    rows: dict[str, dict[str, float]] = {}
    for g in model.genes:
        seen: set[str] = set()
        for r in g.reactions:
            for m in r.metabolites:
                if m.compartment:
                    seen.add(m.compartment)
        if seen:
            rows[g.id] = {c: 1.0 for c in seen}
    df = pd.DataFrame.from_dict(rows, orient="index").fillna(0.0)
    df.index.name = "gene_id"
    return LocalizationScores(df)


def add_noise(scores: LocalizationScores, fraction: float, seed: int) -> LocalizationScores:
    """For ``fraction`` of genes, replace their score row with a single 1.0 in a random
    *wrong* compartment (everything else 0). Simulates "predictor is confidently wrong".
    """
    if fraction <= 0:
        return scores
    rng = np.random.default_rng(seed)
    df = scores.df.copy()
    compartments = list(df.columns)
    n_to_noise = int(round(fraction * len(df)))
    targets = rng.choice(df.index, size=n_to_noise, replace=False)
    for g in targets:
        # find a wrong compartment (any non-top one) to confidently mis-predict
        true_top = df.loc[g].idxmax() if df.loc[g].max() > 0 else compartments[0]
        candidates = [c for c in compartments if c != true_top]
        wrong = rng.choice(candidates)
        df.loc[g, :] = 0.0
        df.at[g, wrong] = 1.0
    return LocalizationScores(df)


def load_csv_scores(path: Path) -> LocalizationScores:
    """Load a ``gene_id × compartment`` CSV (first column = gene_id)."""
    df = pd.read_csv(path, index_col=0)
    df.index.name = "gene_id"
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return LocalizationScores(df)


# --------------------------------------------------------------------------- benchmark

def run_one_test(
    model_orig: cobra.Model,
    truth: dict[str, str],
    scores: LocalizationScores,
    *,
    default_compartment: str,
    transport_cost: float,
    multi_compartment_penalty: float,
    mip_gap: float | None,
    time_limit: float | None,
) -> dict:
    """One MILP solve + accuracy summary.

    Flattens the model to a single compartment (using the curated default as the merged
    id, so reactions truly *in* the default appear unmoved when correctly predicted),
    runs ``predict_localization`` on every truth-set reaction, and returns metrics +
    per-reaction predictions.
    """
    flat, _, _ = merge_compartments(
        model_orig, merged_id=default_compartment, merged_name=default_compartment,
        drop_single_metabolite_reactions=False, deduplicate_reactions=False,
    )
    # The flattened model may have lost some reactions if their net stoichiometry
    # cancelled after the merge — restrict the truth set to surviving reactions.
    surviving = {r.id for r in flat.reactions}
    relevant = {rid: c for rid, c in truth.items() if rid in surviving}

    t = time.time()
    proposal = predict_localization(
        flat, scores, list(relevant),
        default_compartment=default_compartment,
        transport_cost=transport_cost,
        multi_compartment_penalty=multi_compartment_penalty,
        apply=False, mip_gap=mip_gap, time_limit=time_limit,
    )
    elapsed = time.time() - t

    # `moved` only lists reactions whose chosen compartment differs from the flattened
    # `from_compartment` (i.e. `default_compartment`). Anything not in `moved` was
    # placed in the default — record it as such.
    moved_to = dict(zip(proposal.moved["rxn_id"], proposal.moved["to_compartment"], strict=True))
    predictions = {rid: moved_to.get(rid, default_compartment) for rid in relevant}

    correct = sum(predictions[rid] == c for rid, c in relevant.items())
    unplaced = set(proposal.unplaced_reactions) & set(relevant)
    return {
        "seconds": elapsed,
        "n_total": len(relevant),
        "n_correct": correct,
        "n_unplaced": len(unplaced),
        "accuracy": correct / len(relevant) if relevant else 0.0,
        "predictions": predictions,
        "truth": relevant,
    }


def confusion_matrix(predictions: dict[str, str], truth: dict[str, str]) -> pd.DataFrame:
    """Tidy `true × predicted` count matrix."""
    rows = pd.DataFrame({
        "true": [truth[r] for r in predictions],
        "predicted": list(predictions.values()),
    })
    cm = rows.groupby(["true", "predicted"]).size().unstack(fill_value=0)
    return cm.sort_index().sort_index(axis=1)


# --------------------------------------------------------------------------- main

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--yeast-gem", type=Path,
                    default=Path.home() / "github" / "pcSecYeastSpecies" / "Model" / "yeastGEM.xml")
    ap.add_argument("--scores-csv", type=Path,
                    help="optional gene_id × compartment CSV; defaults to from-model scores")
    ap.add_argument("--noise", default="0,0.1,0.25,0.5",
                    help="comma-separated noise fractions to sweep (ignored with --scores-csv)")
    ap.add_argument("--default-compartment", default="c")
    ap.add_argument("--transport-cost", type=float, default=0.5)
    ap.add_argument("--multi-compartment-penalty", type=float, default=0.5)
    ap.add_argument("--mip-gap", type=float, default=0.01)
    ap.add_argument("--time-limit", type=float, default=900)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--max-reactions", type=int, default=None,
                    help="optionally subsample the truth set to N reactions (keeps the "
                         "compartment distribution, drawn deterministically with --seed)")
    ap.add_argument("--doc", type=Path, help="write a markdown summary here")
    args = ap.parse_args()

    print(f"loading {args.yeast_gem} ...", flush=True)
    model = cobra.io.read_sbml_model(str(args.yeast_gem))
    truth = build_truth(model)
    print(f"yeast-GEM: {len(model.reactions)} reactions, {len(model.genes)} genes, "
          f"{len(model.compartments)} compartments; truth set: {len(truth)} reactions",
          flush=True)
    if args.max_reactions and args.max_reactions < len(truth):
        # Stratified subsample: keep the original compartment distribution.
        rng = np.random.default_rng(args.seed)
        by_comp: dict[str, list[str]] = {}
        for rid, c in truth.items():
            by_comp.setdefault(c, []).append(rid)
        keep: list[str] = []
        for rids in by_comp.values():
            n = max(1, round(args.max_reactions * len(rids) / len(truth)))
            keep += list(rng.choice(rids, size=min(n, len(rids)), replace=False))
        truth = {rid: truth[rid] for rid in keep}
        print(f"subsampled truth set to {len(truth)} reactions "
              f"(--max-reactions={args.max_reactions})", flush=True)

    base_scores: LocalizationScores
    if args.scores_csv:
        print(f"loading scores from {args.scores_csv} ...", flush=True)
        base_scores = load_csv_scores(args.scores_csv)
        noise_levels = [0.0]  # external scores: no synthetic noise sweep
    else:
        print("deriving reference scores from yeast-GEM ...", flush=True)
        base_scores = derive_scores_from_model(model)
        noise_levels = [float(x) for x in args.noise.split(",")]

    results: list[dict] = []
    for noise in noise_levels:
        scores = add_noise(base_scores, noise, args.seed) if noise > 0 else base_scores
        print(f"\n=== noise={noise:.2f} ({int(noise * len(base_scores.df))} genes "
              f"confidently mis-scored) ===", flush=True)
        r = run_one_test(
            model, truth, scores,
            default_compartment=args.default_compartment,
            transport_cost=args.transport_cost,
            multi_compartment_penalty=args.multi_compartment_penalty,
            mip_gap=args.mip_gap, time_limit=args.time_limit,
        )
        r["noise"] = noise
        results.append(r)
        print(f"  solved in {r['seconds']:.0f}s — accuracy {r['n_correct']}/{r['n_total']} = "
              f"{r['accuracy']:.3f} ({r['n_unplaced']} unplaced)", flush=True)

    # --- Reporting -------------------------------------------------------------
    lines: list[str] = []
    lines += ["# yeast-GEM localisation benchmark", "",
              f"Model: `{args.yeast_gem.name}` — {len(model.reactions)} reactions, "
              f"{len(model.genes)} genes, {len(model.compartments)} compartments. "
              f"Truth set: {len(truth)} single-compartment GPR-annotated reactions. "
              f"Default compartment for the merged model: `{args.default_compartment}`. "
              f"`transport_cost={args.transport_cost}`, "
              f"`multi_compartment_penalty={args.multi_compartment_penalty}`, "
              f"`mip_gap={args.mip_gap}`, `time_limit={args.time_limit}s`.", "",
              "## Accuracy vs. predictor noise", "",
              "| noise | seconds | n_total | n_correct | n_unplaced | accuracy |",
              "|------:|--------:|--------:|----------:|-----------:|---------:|"]
    for r in results:
        lines.append(
            f"| {r['noise']:.2f} | {r['seconds']:.0f} | {r['n_total']} | "
            f"{r['n_correct']} | {r['n_unplaced']} | {r['accuracy']:.3f} |"
        )
    lines.append("")

    # Confusion matrix for the lowest-noise run (typically the most informative).
    best = min(results, key=lambda x: x["noise"])
    cm = confusion_matrix(best["predictions"], best["truth"])
    lines += [f"## Confusion matrix at noise={best['noise']:.2f}", "",
              "Rows = curated (true) compartment; columns = predicted.", ""]
    lines.append("| true \\ pred | " + " | ".join(str(c) for c in cm.columns) + " |")
    lines.append("|---" + "|---" * len(cm.columns) + "|")
    for true_c, row in cm.iterrows():
        lines.append(f"| **{true_c}** | " + " | ".join(str(int(v)) for v in row) + " |")
    lines.append("")

    # Per-compartment accuracy at the lowest-noise run.
    per_comp: dict[str, tuple[int, int]] = {}
    for rid, true_c in best["truth"].items():
        n_true, n_correct = per_comp.get(true_c, (0, 0))
        per_comp[true_c] = (n_true + 1, n_correct + (best["predictions"][rid] == true_c))
    lines += [f"## Per-compartment accuracy at noise={best['noise']:.2f}", "",
              "| compartment | n | n_correct | accuracy |",
              "|---|--:|--:|--:|"]
    for c in sorted(per_comp):
        n, ok = per_comp[c]
        lines.append(f"| {c} | {n} | {ok} | {ok / n:.3f} |")
    lines.append("")

    text = "\n".join(lines) + "\n"
    print("\n" + text)
    if args.doc:
        args.doc.write_text(text)
        print(f"wrote {args.doc}", flush=True)


if __name__ == "__main__":
    main()
