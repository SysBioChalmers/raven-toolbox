#!/usr/bin/env python3
"""Ablation benchmarks for two `assign_compartments` features that were previously only
correctness-tested, never measured on a real model: **transport-reaction pruning** and **gap-filling**.

Run on curated *S. cerevisiae* yeast-GEM, flattened to one compartment and reassigned (the same draft the
other yeast studies use). Three parts, one JSON output:

1. **Transport pruning** (`prune_transports=True` vs `False`). How many transport reactions pruning
   removes, whether the pruned set is a strict subset of the unpruned one, whether reaction-level
   placement accuracy moves, and the runtime cost.

2. **Gap-fill, natural draft** (`universal=` passed vs not, no induced gaps). Whether gap-fill fires at
   all when the flattened draft is simply reassigned — i.e. does the feature add reactions gratuitously
   at genome scale, or does transport addition alone restore growth?

3. **Gap-fill, knockout-recovery** (ground-truthed). Remove known growth-essential reactions from the
   certified model one at a time and gap-fill from a universal that contains them (the same
   `cobra.flux_analysis.gapfill` call `assign_compartments` uses internally). Of the reactions gap-fill
   adds, how many restore growth, and how many are the *exact* removed reaction — the precision of the
   additions against ground truth. cobra's gapfill has a known numerical-validation failure mode; its
   rate is reported honestly, since `assign_compartments` degrades to "no gap-fill" when it trips.

ASCII-only output. Deterministic (seeded knockout sample).
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
from cobra.flux_analysis import gapfill as cobra_gapfill  # noqa: E402
from cobra.flux_analysis import single_reaction_deletion  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_replicate_yeast_gem import (  # noqa: E402
    _name,
    _norm,
    build_draft,
    curated_reaction_compartments,
    load_yeast_scores,
)

from raven_toolbox.localization import apply_assignment, assign_compartments  # noqa: E402


def _run(draft, biomass_id, scores, relocate, *, min_growth, time_limit, **kw):
    """One assign_compartments -> apply_assignment pass; returns (proposal, applied, growth, wall)."""
    t0 = time.monotonic()
    prop = assign_compartments(draft, scores, relocate, default_compartment="c", base_metabolite=_name,
                               biomass_reaction=biomass_id, min_growth=min_growth,
                               time_limit=time_limit, **kw)
    wall = round(time.monotonic() - t0, 1)
    applied = apply_assignment(draft, prop, default_compartment="c", base_metabolite=_name,
                               universal=kw.get("universal"))
    applied.objective = biomass_id
    growth = round(applied.slim_optimize(error_value=0.0) or 0.0, 4)
    return prop, applied, growth, wall


def _agreement(prop, curated_rxn):
    rc = {rid: _norm(cs[0]) for rid, cs in prop.placements.items() if cs}
    common = set(rc) & set(curated_rxn)
    return (round(sum(rc[r] == curated_rxn[r] for r in common) / len(common), 4) if common else None,
            len(common))


def part_prune(draft, biomass_id, scores, relocate, curated_rxn, *, min_growth, time_limit):
    """prune_transports True vs False: transport count, subset relation, accuracy, runtime."""
    on, _, g_on, w_on = _run(draft, biomass_id, scores, relocate, min_growth=min_growth,
                             time_limit=time_limit, prune_transports=True)
    off, _, g_off, w_off = _run(draft, biomass_id, scores, relocate, min_growth=min_growth,
                                time_limit=time_limit, prune_transports=False)
    t_on, t_off = set(on.added_transports), set(off.added_transports)
    acc_on, n_on = _agreement(on, curated_rxn)
    acc_off, n_off = _agreement(off, curated_rxn)
    removed = len(t_off) - len(t_on)
    return {
        "transports_pruned": len(t_on), "transports_unpruned": len(t_off),
        "removed_by_pruning": removed,
        "removed_fraction": round(removed / len(t_off), 4) if t_off else None,
        "pruned_is_strict_subset": t_on <= t_off,          # pruning only ever removes
        "only_in_pruned": len(t_on - t_off),               # expect 0
        "agreement_pruned": acc_on, "agreement_unpruned": acc_off, "agreement_n": n_on,
        "accuracy_unchanged": acc_on == acc_off,
        "growth_pruned": g_on, "growth_unpruned": g_off,
        "wall_pruned_s": w_on, "wall_unpruned_s": w_off,
    }


def part_gapfill_natural(draft, biomass_id, scores, relocate, universal, *, min_growth, time_limit):
    """Does gap-fill fire when the flattened draft is simply reassigned (no induced gaps)?"""
    no_u, _, g_no, _ = _run(draft, biomass_id, scores, relocate, min_growth=min_growth,
                            time_limit=time_limit)
    with_u, _, g_u, _ = _run(draft, biomass_id, scores, relocate, min_growth=min_growth,
                             time_limit=time_limit, universal=universal)
    return {
        "added_reactions_without_universal": len(no_u.added_reactions),
        "added_reactions_with_universal": len(with_u.added_reactions),
        "certified_without_universal": no_u.certified,
        "certified_with_universal": with_u.certified,
        "growth_without_universal": g_no, "growth_with_universal": g_u,
        "gratuitous_gapfill": len(with_u.added_reactions) > 0,   # expect False
    }


def part_gapfill_knockout(applied, *, floor_fraction, sample_size, seed):
    """Ground-truthed recovery: remove each essential reaction, gap-fill from a universal that has it.

    The universal is a copy of the certified model, so every removed reaction is a candidate and the
    ground truth is exact. Uses the same cobra.flux_analysis.gapfill call assign_compartments._gapfill
    wraps. Reports recovery rate, exact-match rate, and cobra's numerical-failure rate."""
    import random

    g0 = applied.slim_optimize()
    floor = floor_fraction * g0
    dl = single_reaction_deletion(applied, processes=1).reset_index(drop=True)
    essential = []
    for _, row in dl.iterrows():
        ids = row["ids"]
        if len(ids) != 1:
            continue
        rid = next(iter(ids))
        gr = row["growth"]
        if (gr is None or gr != gr or gr < floor) and not applied.reactions.get_by_id(rid).boundary:
            essential.append(rid)

    universal = applied.copy()
    sample = random.Random(seed).sample(essential, min(sample_size, len(essential)))
    recovered = exact = cobra_failure = added_nothing = 0
    for rid in sample:
        with applied:
            applied.remove_reactions([applied.reactions.get_by_id(rid)])
            try:
                sols = cobra_gapfill(applied, universal, lower_bound=max(floor, 1e-4),
                                     demand_reactions=False, iterations=1)
                added = [r.id for r in sols[0]] if sols else []
            except Exception:  # noqa: BLE001 — cobra tolerance/backend quirk; _gapfill swallows it too
                cobra_failure += 1
                continue
        if not added:
            added_nothing += 1
            continue
        recovered += 1
        exact += rid in added
    return {
        "internal_essential_reactions": len(essential),
        "sampled": len(sample),
        "recovered_growth": recovered,
        "recovery_rate": round(recovered / len(sample), 4) if sample else None,
        "exact_reaction_readded": exact,
        "exact_rate_of_recovered": round(exact / recovered, 4) if recovered else None,
        "cobra_gapfill_numerical_failures": cobra_failure,
        "added_nothing": added_nothing,
        "floor_fraction": floor_fraction,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--yeast-gem", type=Path, default=Path("C:/Work/GitHub/yeast-GEM/model/yeast-GEM.xml"))
    ap.add_argument("--data-dir", type=Path, default=Path("data/deeploc"))
    ap.add_argument("--min-growth-fraction", type=float, default=0.5)
    ap.add_argument("--time-limit", type=float, default=120.0)
    ap.add_argument("--knockout-sample", type=int, default=50)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--out", type=Path, default=Path(".research_tmp/assignment_ablations.json"))
    args = ap.parse_args(argv)

    yeast = cobra.io.read_sbml_model(str(args.yeast_gem))
    curated_growth = yeast.slim_optimize()
    curated_rxn = curated_reaction_compartments(yeast)
    draft, biomass_id = build_draft(yeast)
    scores = load_yeast_scores(args.data_dir)
    relocate = [r.id for r in draft.reactions if not r.boundary and r.id != biomass_id]
    min_growth = args.min_growth_fraction * curated_growth
    print(f"curated growth {curated_growth:.4f}; draft {len(draft.reactions)} reactions; "
          f"relocate {len(relocate)}; min_growth {min_growth:.4f}", flush=True)

    out: dict = {"min_growth": round(min_growth, 4), "curated_growth": round(curated_growth, 4)}

    print("\n[1/3] transport-pruning ablation ...", flush=True)
    out["pruning"] = part_prune(draft, biomass_id, scores, relocate, curated_rxn,
                                min_growth=min_growth, time_limit=args.time_limit)
    p = out["pruning"]
    print(f"  pruned {p['transports_pruned']} vs unpruned {p['transports_unpruned']} transports "
          f"({p['removed_by_pruning']} removed, {p['removed_fraction']:.1%}); "
          f"strict subset={p['pruned_is_strict_subset']}")
    print(f"  reaction agreement pruned {p['agreement_pruned']} vs unpruned {p['agreement_unpruned']} "
          f"(unchanged={p['accuracy_unchanged']}); runtime {p['wall_pruned_s']}s vs {p['wall_unpruned_s']}s")

    print("\n[2/3] gap-fill on the natural draft (no induced gaps) ...", flush=True)
    out["gapfill_natural"] = part_gapfill_natural(draft, biomass_id, scores, relocate, universal=draft,
                                                  min_growth=min_growth, time_limit=args.time_limit)
    g = out["gapfill_natural"]
    print(f"  added reactions: without universal {g['added_reactions_without_universal']}, "
          f"with universal {g['added_reactions_with_universal']} "
          f"(gratuitous gap-fill={g['gratuitous_gapfill']})")

    print("\n[3/3] gap-fill knockout-recovery (ground-truthed) ...", flush=True)
    _base, applied, _g, _w = _run(draft, biomass_id, scores, relocate, min_growth=min_growth,
                                  time_limit=args.time_limit)
    out["gapfill_knockout"] = part_gapfill_knockout(applied, floor_fraction=0.05,
                                                    sample_size=args.knockout_sample, seed=args.seed)
    k = out["gapfill_knockout"]
    print(f"  {k['sampled']} essential knockouts: recovered {k['recovered_growth']} "
          f"({k['recovery_rate']:.1%}); of those, {k['exact_reaction_readded']} re-added the exact "
          f"reaction ({k['exact_rate_of_recovered']:.1%})")
    print(f"  cobra.gapfill numerical failures: {k['cobra_gapfill_numerical_failures']}/{k['sampled']}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=str))
    print(f"\nwritten -> {args.out}")


if __name__ == "__main__":
    main()
