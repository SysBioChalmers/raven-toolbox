#!/usr/bin/env python3
"""Authoritative, reproducible yeast benchmark for ``assign_compartments``.

Regenerates the numbers the earlier studies produced with the now-superseded method. Two parts, one
JSON output (``.research_tmp/certified_yeast.json``):

1. **Comparison 1 (vs curated yeast-GEM):** reaction- and gene-level compartment agreement,
   added-transport count, materialised growth, and functional connectivity (blocked-reaction fraction).
2. **Comparison 2 Leg A (vs CarveFungi):** gene-level localisation accuracy at CarveFungi's four
   categories (ER/mito/peroxisome/other) on the genes all three share, **with a McNemar exact test and a
   bootstrap 95% CI on the accuracy difference** — so the claim is honestly "matches" or "beats", not a
   bare point estimate.

The certified method is score-driven and does not use ``transport_cost``, so no proteome annotation is
needed (scalar cost). ASCII-only output.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import cobra  # noqa: E402
from cobra.flux_analysis import find_blocked_reactions  # noqa: E402
from scipy.stats import binomtest  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_replicate_yeast_gem import (  # noqa: E402
    _name,
    _norm,
    build_draft,
    curated_gene_compartments,
    curated_reaction_compartments,
    load_yeast_scores,
)
from benchmark_vs_carvefungi import (  # noqa: E402
    _CF_COARSE,
    _YEASTGEM_COARSE,
    _coarse,
    gene_compartments_coarse,
)

from raven_toolbox.localization import (  # noqa: E402
    apply_assignment,
    assign_compartments,
)


def comparison1(proposal, draft, relocate, biomass_id, curated_rxn, curated_gene):
    res_comp = {rid: _norm(cs[0]) for rid, cs in proposal.placements.items() if cs}
    common = set(res_comp) & set(curated_rxn)
    rxn = sum(res_comp[r] == curated_rxn[r] for r in common) / len(common) if common else None
    res_gene = {g: {_norm(c) for c in cs} for g, cs in proposal.gene_compartments.items() if cs}
    common_g = set(res_gene) & set(curated_gene)
    gene = sum(bool(res_gene[g] & curated_gene[g]) for g in common_g) / len(common_g) if common_g else None
    applied = apply_assignment(draft, proposal, default_compartment="c", base_metabolite=_name)
    applied.objective = biomass_id
    growth = applied.slim_optimize(error_value=0.0) or 0.0
    blocked = set(find_blocked_reactions(applied)) & set(relocate)
    return {"reaction_agreement": round(rxn, 4), "reaction_n": len(common),
            "gene_agreement": round(gene, 4), "gene_n": len(common_g),
            "transports": len(proposal.added_transports), "growth": round(growth, 4),
            "blocked_fraction": round(len(blocked) / len(relocate), 4)}


def leg_a_with_stats(cf_model, yeast, ours_gene_comps, seed):
    truth = gene_compartments_coarse(yeast, _YEASTGEM_COARSE)
    cf = gene_compartments_coarse(cf_model, _CF_COARSE)
    ours = {g: {_coarse(c, _YEASTGEM_COARSE) for c in cs} for g, cs in ours_gene_comps.items() if cs}
    common = sorted(set(truth) & set(cf) & set(ours))
    ours_ok = [bool(ours[g] & truth[g]) for g in common]
    cf_ok = [bool(cf[g] & truth[g]) for g in common]
    our_acc = sum(ours_ok) / len(common)
    cf_acc = sum(cf_ok) / len(common)

    # McNemar exact test on discordant pairs (paired: same genes)
    b = sum(o and not c for o, c in zip(ours_ok, cf_ok, strict=True))  # ours right, CF wrong
    c = sum(c and not o for o, c in zip(ours_ok, cf_ok, strict=True))  # CF right, ours wrong
    p = binomtest(min(b, c), b + c, 0.5).pvalue if (b + c) else 1.0

    # bootstrap 95% CI on (our_acc - cf_acc); deterministic via a seeded index generator
    import random
    rng = random.Random(seed)
    n = len(common)
    diffs = []
    for _ in range(10000):
        idx = [rng.randrange(n) for _ in range(n)]
        diffs.append(sum(ours_ok[i] for i in idx) / n - sum(cf_ok[i] for i in idx) / n)
    diffs.sort()
    lo, hi = diffs[249], diffs[9750]
    verdict = ("beats" if (p < 0.05 and our_acc > cf_acc)
               else "loses" if (p < 0.05 and our_acc < cf_acc) else "matches")
    return {"common_genes": n, "our_accuracy": round(our_acc, 4), "carvefungi_accuracy": round(cf_acc, 4),
            "ours_right_cf_wrong": b, "cf_right_ours_wrong": c, "mcnemar_p": round(p, 4),
            "acc_diff": round(our_acc - cf_acc, 4), "ci95": [round(lo, 4), round(hi, 4)],
            "verdict": verdict}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--yeast-gem", type=Path, default=Path("C:/Work/GitHub/yeast-GEM/model/yeast-GEM.xml"))
    ap.add_argument("--carvefungi-model", type=Path,
                    default=Path(".research_tmp/carvefungi_yeast_model.sbml"))
    ap.add_argument("--data-dir", type=Path, default=Path("data/deeploc"))
    ap.add_argument("--min-growth-fraction", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--out", type=Path, default=Path(".research_tmp/certified_yeast.json"))
    args = ap.parse_args(argv)

    yeast = cobra.io.read_sbml_model(str(args.yeast_gem))
    curated_growth = yeast.slim_optimize()
    curated_rxn = curated_reaction_compartments(yeast)
    curated_gene = curated_gene_compartments(yeast)
    draft, biomass_id = build_draft(yeast)
    scores = load_yeast_scores(args.data_dir)
    relocate = [r.id for r in draft.reactions if not r.boundary and r.id != biomass_id]
    min_growth = args.min_growth_fraction * curated_growth

    print(f"curated growth {curated_growth:.4f}; relocate {len(relocate)}; min_growth {min_growth:.4f}")
    print("running assign_compartments ...", flush=True)
    t0 = time.monotonic()
    prop = assign_compartments(
        draft, scores, relocate, default_compartment="c", base_metabolite=_name,
        biomass_reaction=biomass_id, min_growth=min_growth)
    wall = round(time.monotonic() - t0, 1)
    print(f"  {wall}s  certified={prop.certified}", flush=True)

    out = {"method": "assign_compartments", "wall_s": wall, "certified": prop.certified,
           "curated_growth": round(curated_growth, 4), "min_growth": round(min_growth, 4),
           "comparison1": comparison1(prop, draft, relocate, biomass_id, curated_rxn, curated_gene)}

    c1 = out["comparison1"]
    print("\n=== Comparison 1 (vs curated yeast-GEM) ===")
    print(f"  reaction agreement {c1['reaction_agreement']:.1%} ({c1['reaction_n']}); "
          f"gene agreement {c1['gene_agreement']:.1%} ({c1['gene_n']})")
    print(f"  transports {c1['transports']}; growth {c1['growth']}; "
          f"blocked fraction {c1['blocked_fraction']:.1%}")

    if args.carvefungi_model.is_file():
        cf_model = cobra.io.read_sbml_model(str(args.carvefungi_model))
        ours_gc = {g: cs for g, cs in prop.gene_compartments.items() if cs}
        out["leg_a"] = leg_a_with_stats(cf_model, yeast, ours_gc, args.seed)
        la = out["leg_a"]
        print(f"\n=== Comparison 2 Leg A (vs CarveFungi, 4 categories, {la['common_genes']} shared genes) ===")
        print(f"  ours {la['our_accuracy']:.1%}  CarveFungi {la['carvefungi_accuracy']:.1%}  "
              f"diff {la['acc_diff']:+.1%} (95% CI [{la['ci95'][0]:+.1%}, {la['ci95'][1]:+.1%}])")
        print(f"  discordant: ours-right/CF-wrong {la['ours_right_cf_wrong']}, "
              f"CF-right/ours-wrong {la['cf_right_ours_wrong']}; McNemar p={la['mcnemar_p']}")
        print(f"  VERDICT: certified {la['verdict'].upper()} CarveFungi on gene-level accuracy")
    else:
        print(f"\n[skip Leg A] CarveFungi model not at {args.carvefungi_model}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=str))
    (args.out.with_name("certified_yeast_placements.json")).write_text(json.dumps(
        {"gene_compartments": {g: cs for g, cs in prop.gene_compartments.items() if cs},
         "reaction_placements": {r: cs for r, cs in prop.placements.items() if cs}}, default=str))
    print(f"\nwritten -> {args.out}")


if __name__ == "__main__":
    main()
