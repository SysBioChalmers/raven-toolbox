#!/usr/bin/env python3
"""Does constraint ROW order change which co-optimal vertex Gurobi picks?
Rebuild the placement model with its constraints in a permuted order (variables
untouched) and compare the chosen placement to the reference. If it differs, the
MATLAB and Python live masters must present rows in the same order to agree."""
from __future__ import annotations

import hashlib
from pathlib import Path

import gurobipy as gp

MPS = Path(".research_tmp/master_seed1.mps")
PARAMS = dict(OutputFlag=0, Threads=1, Seed=0, MIPGap=0.0, Presolve=2,
              IntFeasTol=1e-9, FeasibilityTol=1e-9, OptimalityTol=1e-9)


def chosen_hash(m):
    chosen = sorted(v.VarName for v in m.getVars()
                    if v.VarName.startswith("x_") and v.X > 0.5)
    return hashlib.sha256("\n".join(chosen).encode()).hexdigest()[:12]


def solve_permuted(perm_key):
    with gp.Env(empty=True) as env:
        for k, v in PARAMS.items():
            env.setParam(k, v)
        env.start()
        src = gp.read(str(MPS), env)
        src.update()
        dst = gp.Model("perm", env)
        # copy variables in the SAME column order
        vmap = {}
        for v in src.getVars():
            vmap[v.VarName] = dst.addVar(lb=v.LB, ub=v.UB, obj=v.Obj,
                                         vtype=v.VType, name=v.VarName)
        dst.ModelSense = src.ModelSense
        dst.update()
        # copy constraints in a permuted ROW order
        cons = list(src.getConstrs())
        order = sorted(range(len(cons)), key=lambda i: perm_key(cons[i].ConstrName, i))
        for i in order:
            c = cons[i]
            row = src.getRow(c)
            expr = gp.LinExpr()
            for j in range(row.size()):
                expr.add(vmap[row.getVar(j).VarName], row.getCoeff(j))
            dst.addLConstr(expr, c.Sense, c.RHS, name=c.ConstrName)
        dst.update()
        dst.optimize()
        return chosen_hash(dst), dst.ObjVal


ref_h, ref_obj = solve_permuted(lambda name, i: i)                     # identity order
rev_h, _ = solve_permuted(lambda name, i: -i)                          # fully reversed
name_h, _ = solve_permuted(lambda name, i: name)                       # sorted by name
print(f"identity row order : hash={ref_h} obj={ref_obj:.6f}")
print(f"reversed row order : hash={rev_h}  {'SAME' if rev_h==ref_h else 'DIFFERENT'}")
print(f"name-sorted order  : hash={name_h}  {'SAME' if name_h==ref_h else 'DIFFERENT'}")
