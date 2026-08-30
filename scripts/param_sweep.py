#!/usr/bin/env python3
"""Gurobi parameter benchmark on the byte-identical yeast placement MILP.

Solves the SAME model (a fixed MPS) under a grid of parameter settings and reports, for each:
solve wall-time, objective, and how much the chosen co-optimal placement moves relative to a fixed
reference (Threads=1, Seed=0, MIPGap=0, Presolve=2). Because the model is identical every time, any
change in the placement is caused purely by the parameter -- this quantifies which parameters must be
pinned for reproducibility and what each costs in speed. ASCII-only.
"""
from __future__ import annotations

import hashlib
import time
from pathlib import Path

import gurobipy as gp

MPS = Path(".research_tmp/master_seed1.mps")
REF = dict(Threads=1, Seed=0, MIPGap=0.0, Presolve=2,
           IntFeasTol=1e-9, FeasibilityTol=1e-9, OptimalityTol=1e-9)


def placement_hash(model) -> tuple[str, float, int]:
    xs = {v.VarName: v.X for v in model.getVars() if v.VarName.startswith("x_")}
    chosen = sorted(name for name, val in xs.items() if val > 0.5)
    h = hashlib.sha256("\n".join(chosen).encode()).hexdigest()[:12]
    return h, model.ObjVal, len(chosen)


def solve(params: dict) -> dict:
    with gp.Env(empty=True) as env:
        env.setParam("OutputFlag", 0)
        env.start()
        m = gp.read(str(MPS), env)
        for k, v in params.items():
            m.setParam(k, v)
        t0 = time.perf_counter()
        m.optimize()
        wall = time.perf_counter() - t0
        h, obj, n = placement_hash(m)
        return dict(wall=wall, obj=obj, n=n, hash=h, status=m.Status)


def dump_reference_names():
    """Write the reference vertex's chosen x-variables (sorted) for MATLAB cross-check."""
    with gp.Env(empty=True) as env:
        env.setParam("OutputFlag", 0)
        env.start()
        m = gp.read(str(MPS), env)
        for k, v in REF.items():
            m.setParam(k, v)
        m.optimize()
        chosen = sorted(v.VarName for v in m.getVars()
                        if v.VarName.startswith("x_") and v.X > 0.5)
        Path(".research_tmp/ref_placement.txt").write_text("\n".join(chosen) + "\n")
        print(f"wrote .research_tmp/ref_placement.txt ({len(chosen)} vars)")


def run():
    print(f"model: {MPS}  ({MPS.stat().st_size} bytes)")
    ref = solve(REF)
    print(f"\nREFERENCE  Threads=1 Seed=0 MIPGap=0 Presolve=2  "
          f"-> obj={ref['obj']:.6f} n={ref['n']} hash={ref['hash']} wall={ref['wall']:.1f}s\n")
    dump_reference_names()

    # one-parameter-at-a-time sweeps from the reference
    sweeps = {
        "Threads":        [1, 2, 4, 8],
        "Seed":           [0, 1, 7, 12345],
        "MIPGap":         [0.0, 1e-9, 1e-6, 1e-4, 1e-2],
        "Presolve":       [-1, 0, 1, 2],
        "MIPFocus":       [0, 1, 2, 3],
        "Heuristics":     [0.0, 0.05, 0.5],
        "Method":         [-1, 0, 1, 2],
        "FeasibilityTol": [1e-9, 1e-7, 1e-6],
        "OptimalityTol":  [1e-9, 1e-7, 1e-6],
        "IntFeasTol":     [1e-9, 1e-6, 1e-5],
    }
    print(f"{'parameter':<12} {'value':>8}  {'wall(s)':>8} {'obj':>12} {'n':>5}  {'hash':>12}  agree%")
    print("-" * 74)
    for pname, values in sweeps.items():
        for val in values:
            params = dict(REF)
            params[pname] = val
            r = solve(params)
            agree = "REF" if r["hash"] == ref["hash"] else "DIFF"
            same_obj = "" if abs(r["obj"] - ref["obj"]) < 1e-6 else " obj!"
            print(f"{pname:<12} {str(val):>8}  {r['wall']:>8.2f} {r['obj']:>12.6f} "
                  f"{r['n']:>5}  {r['hash']:>12}  {agree}{same_obj}")
        print()


if __name__ == "__main__":
    run()
