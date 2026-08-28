#!/usr/bin/env python3
"""Comparison 2 — is our evidence-aware assignment better than CarveFungi?

Two angles the user asked for, both using only the shipped public API on our side and CarveFungi's own
model on the other side (produced by `build_carvefungi_yeast_model.py`, which drives CarveFungi's real
`bin/` code):

**Leg A — both methods vs curated yeast-GEM, at the gene level.** Both CarveFungi and our approach turn
the *same* yeast DeepLoc evidence into a per-gene compartment call. Score each against curated
yeast-GEM's own gene->compartment truth. Because CarveFungi only resolves four categories (ER /
mitochondrion / peroxisome / other), all three sides are collapsed to that granularity for a fair
head-to-head. This is the direct "which method places genes more accurately" answer.

**Leg B — our placement vs CarveFungi's own placement.** Flatten CarveFungi's carved model with the
real `merge_compartments`, run our pipeline on it, and measure how closely we reproduce CarveFungi's
own per-reaction compartment choice (same UF##### reaction namespace). This is a "how close are the two
methods to each other" angle, plus a functional-connectivity read on our placement of CarveFungi's
network.

Our Leg-A gene calls are reused from Comparison 1's saved placements
(`.research_tmp/replicate_yeast_gem_placements.json`) so the genome-scale yeast-GEM MILP is not re-run.

Usage::

    python scripts/benchmark_vs_carvefungi.py \\
        --carvefungi-model .research_tmp/carvefungi_yeast_model.sbml \\
        --yeast-gem C:/Work/GitHub/yeast-GEM/model/yeast-GEM.xml

ASCII-only output.
"""
from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import cobra  # noqa: E402
import pandas as pd  # noqa: E402
from benchmark_replicate_yeast_gem import annotate_yeast_proteome  # noqa: E402
from cobra.flux_analysis import find_blocked_reactions  # noqa: E402

from raven_toolbox.localization import (  # noqa: E402
    LocalizationScores,
    SubstrateOntology,
    apply_assignment,
    assign_compartments,
    default_substrate_of,
    evidence_aware_transport_cost,
    load_deeploc,
)
from raven_toolbox.manipulation.compartments import merge_compartments  # noqa: E402

# CarveFungi's universal-DB compartment codes (c, m, r=ER, x=peroxisome, n, e) -- distinct from
# yeast-GEM's. DeepLoc labels with no CarveFungi home (Golgi, vacuole, ...) map to nothing and are
# dropped by load_deeploc.
CF_COMPARTMENT_MAP = {
    "cytoplasm": "c", "cytosol": "c", "nucleus": "n", "nucleoplasm": "n",
    "mitochondrion": "m", "mitochondria": "m", "mitochondrial": "m",
    "peroxisome": "x", "endoplasmic reticulum": "r",
    "extracellular": "e", "extracellular space": "e", "extracellular region": "e", "secreted": "e",
}

# Four coarse categories CarveFungi natively resolves; the fair granularity for the gene head-to-head.
_YEASTGEM_COARSE = {"er": "ER", "erm": "ER", "m": "mito", "mm": "mito", "p": "peroxisome"}
_CF_COARSE = {"r": "ER", "m": "mito", "x": "peroxisome"}


def _coarse(compartment: str, table: dict[str, str]) -> str:
    return table.get(compartment, "other")


def _sole_compartment(rxn: cobra.Reaction) -> str | None:
    comps = {m.compartment for m in rxn.metabolites}
    return next(iter(comps)) if len(comps) == 1 else None


def gene_compartments_coarse(model: cobra.Model, table: dict[str, str]) -> dict[str, set[str]]:
    """gene -> the coarse compartments of the single-compartment reactions it catalyses."""
    out: dict[str, set[str]] = {}
    for r in model.reactions:
        if r.boundary:
            continue
        c = _sole_compartment(r)
        if c is None:
            continue
        for g in r.genes:
            out.setdefault(g.id, set()).add(_coarse(c, table))
    return out


def _accuracy(pred: dict[str, set[str]], truth: dict[str, set[str]], genes: set[str]) -> float:
    ok = sum(1 for g in genes if pred[g] & truth[g])
    return ok / len(genes) if genes else float("nan")


def leg_a(cf_model, yeast, our_placements: dict) -> dict:
    """Gene-level accuracy vs curated yeast-GEM, ours vs CarveFungi, at ER/mito/peroxisome/other."""
    truth = gene_compartments_coarse(yeast, _YEASTGEM_COARSE)
    cf = gene_compartments_coarse(cf_model, _CF_COARSE)
    # our per-gene call comes from Comparison 1's saved yeast-GEM-code placements
    ours = {g: {_coarse(c, _YEASTGEM_COARSE) for c in cs}
            for g, cs in our_placements["gene_compartments"].items() if cs}

    common = set(truth) & set(cf) & set(ours)
    our_acc = _accuracy(ours, truth, common)
    cf_acc = _accuracy(cf, truth, common)
    # head-to-head on the genes where the two methods disagree
    disagree = {g for g in common if ours[g] != cf[g]}
    ours_win = sum(1 for g in disagree if (ours[g] & truth[g]) and not (cf[g] & truth[g]))
    cf_win = sum(1 for g in disagree if (cf[g] & truth[g]) and not (ours[g] & truth[g]))

    print("\n=== Leg A: gene-level accuracy vs curated yeast-GEM (ER/mito/peroxisome/other) ===")
    print(f"  common gene set: {len(common)}")
    print(f"  ours:       {our_acc:.1%}")
    print(f"  CarveFungi: {cf_acc:.1%}")
    print(f"  where they disagree ({len(disagree)} genes): ours right/CF wrong {ours_win}, "
          f"CF right/ours wrong {cf_win}")
    return {"common_genes": len(common), "our_accuracy": our_acc, "carvefungi_accuracy": cf_acc,
            "disagree_genes": len(disagree), "ours_win": ours_win, "carvefungi_win": cf_win}


def leg_b(cf_model, data_dir: Path, time_limit: float) -> dict:
    """Our pipeline on CarveFungi's own network: how closely do we reproduce its placement, and is our
    placement functional?"""
    print("\n=== Leg B: our placement vs CarveFungi's own (same reaction set) ===")
    cf_growth = cf_model.slim_optimize()
    biomass_id = next((r.id for r in cf_model.reactions if r.objective_coefficient != 0), "BIOMASS")
    # CarveFungi's own single-compartment placement (ground truth for this leg)
    cf_place = {r.id: _sole_compartment(r) for r in cf_model.reactions
                if not r.boundary and _sole_compartment(r) is not None}
    cf_transports = sum(1 for r in cf_model.reactions
                        if not r.boundary and r.id != biomass_id and _sole_compartment(r) is None)

    draft, _d, _u = merge_compartments(cf_model, merged_id="c", merged_name="cytoplasm",
                                       drop_single_metabolite_reactions=False,
                                       deduplicate_reactions=False)
    if biomass_id in draft.reactions:
        draft.objective = biomass_id
    print(f"  CarveFungi model: {len(cf_model.reactions)} reactions, growth={cf_growth:.4f}; "
          f"flattened draft: {len(draft.reactions)} reactions, growth={draft.slim_optimize():.4f}")

    frames = [load_deeploc(c, compartment_map=CF_COMPARTMENT_MAP).df
              for c in sorted(data_dir.glob("yeast-GEM_deeploc_*.csv"))]
    scores = LocalizationScores(pd.concat(frames))
    ann = annotate_yeast_proteome(data_dir)
    gene_comps = {g: {scores.df.loc[g].astype(float).idxmax()} for g in scores.df.index if g in ann}
    cost = evidence_aware_transport_cost(
        draft, ann, gene_comps, substrate_of=default_substrate_of, ontology=SubstrateOntology.load(),
        base_cost=0.5)

    relocate = [r.id for r in draft.reactions if not r.boundary and r.id != biomass_id]
    min_growth = 0.5 * (draft.slim_optimize() or 0.0)
    t0 = time.monotonic()
    proposal = assign_compartments(draft, scores, relocate, transport_cost=cost,
                                   default_compartment="c", biomass_reaction=biomass_id,
                                   min_growth=min_growth, time_limit=time_limit)
    wall = time.monotonic() - t0
    result = apply_assignment(draft, proposal, default_compartment="c")
    if biomass_id in result.reactions:
        result.objective = biomass_id

    res_place = {rid: cs[0] for rid, cs in proposal.placements.items() if cs}
    common = set(res_place) & set(cf_place)
    agree = sum(1 for rid in common if res_place[rid] == cf_place[rid])
    rate = agree / len(common) if common else float("nan")
    growth = result.slim_optimize()
    blocked = set(find_blocked_reactions(result)) & set(relocate)
    blocked_rate = len(blocked) / len(relocate) if relocate else float("nan")

    print(f"  status={proposal.status} objective={proposal.objective:.3f} solve={wall:.1f}s")
    print(f"  agreement with CarveFungi's placement: {rate:.1%} ({agree}/{len(common)})")
    print(f"  transports added: {len(proposal.added_transports)} (CarveFungi: {cf_transports}); "
          f"blocked re-placed: {blocked_rate:.1%}; growth {growth}")
    return {"status": proposal.status, "solve_seconds": wall,
            "placement_agreement_rate": rate, "placement_agreement_n": len(common),
            "transports_added": len(proposal.added_transports), "carvefungi_transports": cf_transports,
            "blocked_replaced_rate": blocked_rate, "growth": growth,
            "carvefungi_growth": cf_growth}


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--carvefungi-model", required=True, type=Path,
                    help="CarveFungi's carved yeast model (build_carvefungi_yeast_model.py)")
    ap.add_argument("--yeast-gem", type=Path,
                    default=Path("C:/Work/GitHub/yeast-GEM/model/yeast-GEM.xml"))
    ap.add_argument("--our-placements", type=Path,
                    default=Path(".research_tmp/replicate_yeast_gem_results_placements.json"),
                    help="our yeast-GEM gene placements from Comparison 1")
    ap.add_argument("--data-dir", type=Path, default=Path("data/deeploc"))
    ap.add_argument("--time-limit", type=float, default=900.0, help="Leg B per-solve seconds")
    ap.add_argument("--out", type=Path, default=Path(".research_tmp/vs_carvefungi_results.json"))
    args = ap.parse_args(argv)

    cf_model = cobra.io.read_sbml_model(str(args.carvefungi_model))
    yeast = cobra.io.read_sbml_model(str(args.yeast_gem))
    print(f"CarveFungi model: {len(cf_model.reactions)} reactions, {len(cf_model.genes)} genes")
    print(f"curated yeast-GEM: {len(yeast.reactions)} reactions, {len(yeast.genes)} genes")

    results = {}
    if args.our_placements.is_file():
        our_placements = json.loads(args.our_placements.read_text())
        results["leg_a"] = leg_a(cf_model, yeast, our_placements)
    else:
        print(f"\n[skip Leg A] our placements not found at {args.our_placements}; run "
              "benchmark_replicate_yeast_gem.py first")

    results["leg_b"] = leg_b(cf_model, args.data_dir, args.time_limit)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(results, indent=2, default=str))
    print(f"\nwritten -> {args.out}")


if __name__ == "__main__":
    main()
