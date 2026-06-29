#!/usr/bin/env python3
"""Definitive CarveFungi vs ours assignment head-to-head, using CarveFungi's ACTUAL MILP (CPLEX).

This is the tight-gap companion to ``benchmark_carvefungi_milp.py`` (the Gurobi port). The Gurobi port
is faithful but the carve is a hard MILP that Gurobi leaves at ~10-14% gap, where the accuracy
comparison is unstable. CarveFungi runs its ``minmax_reduction`` in **CPLEX** to a 0.1% pool gap, so
running its own code is both maximally faithful and the only way to get a trustworthy accuracy number.

RUN THIS IN A CPLEX-ENABLED PYTHON (3.11/3.12 where ``import cplex`` works — CPLEX 22.2 has no 3.14
binding). No Gurobi needed.

Arms (only the objective's parsimony differs; same candidate set, same scores, same connectivity):
* Arm A (CarveFungi): ``minmax_reduction(model, scores)`` unmodified.
* Arm B (ours): same, with each inter-compartment transport reaction's score reduced by
  ``--transport-cost`` (our transport-minimisation term, fed through their objective).

Inputs:
* ``--carvefungi-dir`` : a clone of github.com/SandraCastilloPriego/CarveFungi (for bin/ + the
  universal model data/reactionDatabase/bigModelv2.21b.sbml + universal_v2.21.csv).
* ``--scores`` : the DeepLoc-injected CarveFungi yeast score dict (reaction_id,score) -- produced by
  running CarveFungi's unmodified scoring with a DeepLoc-derived loc file (see
  docs/studies/carvefungi_milp_benchmark.md). Provided alongside this script.
* ``--yeast-gem`` : yeast-GEM.xml, for the EC-mapped curated gold reference.

ASCII-only output.
"""
from __future__ import annotations

import argparse
import re
import sys
import warnings
from collections import defaultdict
from pathlib import Path

warnings.filterwarnings("ignore")

import cobra  # noqa: E402
import pandas as pd  # noqa: E402

UNIV2YEAST = {"c": "c", "m": "m", "x": "p", "r": "er", "n": "n", "g": "g", "e": "e", "l": "lp"}
COLLAPSE = {"erm": "er", "mm": "m", "gm": "g", "vm": "v"}


def base_id(rid: str) -> str:
    return re.sub(r"_[A-Z]+$", "", rid)


def comp_suffix(rid: str) -> str:
    m = re.search(r"_([A-Z]+)$", rid)
    return m.group(1).lower() if m else ""


def parse_solution(sol: dict[str, float]) -> tuple[dict[str, set[str]], set[str]]:
    """From a CarveFungi solution dict (yf_/yr_ binaries): (base reaction -> compartments kept,
    set of active reaction ids)."""
    kept: dict[str, set[str]] = {}
    active: set[str] = set()
    for name, val in sol.items():
        if (name.startswith("yf_") or name.startswith("yr_")) and val > 0.5:
            rid = name[3:]
            active.add(rid)
            kept.setdefault(base_id(rid), set()).add(comp_suffix(rid))
    return kept, active


def gold_reference(yeast_gem: Path, universal_csv: Path):
    y = cobra.io.read_sbml_model(str(yeast_gem))
    ec2comp: dict[str, set[str]] = defaultdict(set)
    for r in y.reactions:
        if r.boundary or not r.genes:
            continue
        comps = {m.compartment for m in r.metabolites if m.compartment}
        if len(comps) != 1:
            continue
        c = COLLAPSE.get(next(iter(comps)), next(iter(comps)))
        for key in ("ec-code", "eccodes"):
            v = r.annotation.get(key, [])
            for e in (v if isinstance(v, list) else [v]):
                ec2comp[e].add(c)
    u = pd.read_csv(universal_csv)
    base2ec: dict[str, list[str]] = {}
    for rid, ecv in zip(u["IDs"], u["ECs"].fillna(""), strict=True):
        ecs = [e.strip() for e in re.split(r"[;, ]+", str(ecv)) if e.strip() and e[0].isdigit()]
        if ecs:
            base2ec[base_id(str(rid))] = ecs
    return dict(ec2comp), base2ec


def ec_eval(kept, ec2comp, base2ec, restrict=None):
    n = recall = exact = 0
    for base, comps in kept.items():
        if restrict is not None and base not in restrict:
            continue
        ecs = base2ec.get(base, [])
        gold = set().union(*[ec2comp[e] for e in ecs if e in ec2comp]) if ecs else set()
        if not gold:
            continue
        assigned = {UNIV2YEAST.get(c, c) for c in comps}
        n += 1
        recall += bool(assigned & gold)
        exact += (assigned == gold)
    return {"n": n, "recall": recall / max(1, n), "exact": exact / max(1, n)}


def best_solution(solutions):
    if not solutions:
        raise RuntimeError("minmax_reduction returned no solutions (CPLEX populate failed?)")
    return max(solutions, key=lambda s: s.get_obj()).get_solution()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--carvefungi-dir", required=True, type=Path)
    ap.add_argument("--scores", required=True, type=Path, help="reaction_id,score CSV")
    ap.add_argument("--yeast-gem", required=True, type=Path)
    ap.add_argument("--universal", type=Path, help="bigModelv2.21b.sbml (default: under carvefungi-dir)")
    ap.add_argument("--universal-csv", type=Path, help="universal_v2.21.csv (default: under carvefungi-dir)")
    ap.add_argument("--transport-cost", type=float, default=0.3)
    args = ap.parse_args()

    sys.path.insert(0, str(args.carvefungi_dir / "bin"))
    from CarveMeFuncPool import minmax_reduction  # noqa: E402  (needs cplex)

    rdb = args.carvefungi_dir / "data" / "reactionDatabase"
    universal = args.universal or rdb / "bigModelv2.21b.sbml"
    universal_csv = args.universal_csv or rdb / "universal_v2.21.csv"

    model = cobra.io.read_sbml_model(str(universal))
    scores = dict(zip(*[pd.read_csv(args.scores)[c] for c in ("reaction_id", "score")], strict=True))
    is_transport = {r.id: len({m.compartment for m in r.metabolites}) > 1 for r in model.reactions}

    print("Arm A: CarveFungi minmax_reduction (CPLEX, unmodified) ...", flush=True)
    a, act_a = parse_solution(best_solution(minmax_reduction(model, scores)))
    print("Arm B: ours (transports penalised by transport_cost) ...", flush=True)
    scores_b = {rid: (s - args.transport_cost if is_transport.get(rid) else s) for rid, s in scores.items()}
    b, act_b = parse_solution(best_solution(minmax_reduction(model, scores_b)))

    ec2comp, base2ec = gold_reference(args.yeast_gem, universal_csv)
    common = set(a) & set(b)
    ea, eb = ec_eval(a, ec2comp, base2ec), ec_eval(b, ec2comp, base2ec)
    eac, ebc = ec_eval(a, ec2comp, base2ec, common), ec_eval(b, ec2comp, base2ec, common)
    tr_a = sum(1 for rid in act_a if is_transport.get(rid))
    tr_b = sum(1 for rid in act_b if is_transport.get(rid))

    print("\n=== Definitive CPLEX head-to-head (CarveFungi's own MILP) ===")
    print(f"Arm A (CarveFungi): {len(a)} base reactions kept, {tr_a} transports "
          f"({tr_a/max(1,len(a)):.3f}/base)")
    print(f"Arm B (ours):       {len(b)} base reactions kept, {tr_b} transports "
          f"({tr_b/max(1,len(b)):.3f}/base)")
    print(f"common base reactions: {len(common)}; "
          f"same compartment set on common: "
          f"{sum(1 for x in common if a[x] == b[x]) / max(1, len(common)):.1%}")
    print("gold reference (EC-mapped to curated yeast-GEM), each arm's own kept set:")
    print(f"  Arm A: n={ea['n']} recall={ea['recall']:.1%} exact={ea['exact']:.1%}")
    print(f"  Arm B: n={eb['n']} recall={eb['recall']:.1%} exact={eb['exact']:.1%}")
    print("gold reference on the COMMON kept set (fair, same denominator):")
    print(f"  Arm A: n={eac['n']} recall={eac['recall']:.1%} exact={eac['exact']:.1%}")
    print(f"  Arm B: n={ebc['n']} recall={ebc['recall']:.1%} exact={ebc['exact']:.1%}")


if __name__ == "__main__":
    main()
