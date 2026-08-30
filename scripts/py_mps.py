import json, csv
from pathlib import Path
import cobra, gurobipy as gp
from gurobipy import GRB
X = Path(".research_tmp/xarm")
meta = json.loads((X/"meta.json").read_text())
comps = None
scores = {}
with (X/"scores.csv").open() as fh:
    rd = csv.reader(fh); hdr = next(rd); comps = hdr[1:]
    for row in rd: scores[row[0]] = {c: float(v) for c, v in zip(comps, row[1:])}
comps = sorted(comps)
draft = cobra.io.read_sbml_model(str(X/"mat_draft.xml")) if False else cobra.io.read_sbml_model(str(X/"py_draft.xml"))
biomass = meta["biomass_id"]
relocate = set(meta["relocate"])
movable = sorted(r.id for r in draft.reactions if r.id in relocate and not r.boundary and r.id != biomass
                 and len({m.compartment for m in draft.reactions.get_by_id(r.id).metabolites}) == 1) if False else None
# movable = internal, non-boundary, non-biomass, single-compartment; sorted
mov = []
for r in draft.reactions:
    if r.id in relocate and not r.boundary and r.id != biomass:
        cset = {m.compartment for m in r.metabolites if m.compartment}
        if len(cset) == 1: mov.append(r.id)
mov = sorted(mov)
genes = sorted({g.id for rid in mov for g in draft.reactions.get_by_id(rid).genes if g.id in scores})
pen = 0.5
m = gp.Model("place"); m.Params.OutputFlag = 0
x = {(r, c): m.addVar(vtype=GRB.BINARY, name=f"x_{r}_{c}") for r in mov for c in comps}
y = {(g, c): m.addVar(vtype=GRB.BINARY, name=f"y_{g}_{c}") for g in genes for c in comps}
m.update()
gene_rxns = {g: [] for g in genes}
for r in mov:
    for g in {gg.id for gg in draft.reactions.get_by_id(r).genes} & set(genes):
        gene_rxns[g].append(r)
cons = []
for r in mov:
    cons.append((f"place_{r}", gp.quicksum(x[r, c] for c in comps), '=', 1))
    for g in sorted({gg.id for gg in draft.reactions.get_by_id(r).genes} & set(genes)):
        for c in comps:
            cons.append((f"couple_{r}_{g}_{c}", x[r, c] - y[g, c], '<', 0))
for g in genes:
    cons.append((f"gene1_{g}", gp.quicksum(y[g, c] for c in comps), '>', 1))
    for c in comps:
        cons.append((f"has_{g}_{c}", y[g, c] - gp.quicksum(x[r, c] for r in gene_rxns[g]), '<', 0))
for name, lhs, sense, rhs in sorted(cons, key=lambda t: t[0]):
    if sense == '=': m.addConstr(lhs == rhs, name=name)
    elif sense == '<': m.addConstr(lhs <= rhs, name=name)
    else: m.addConstr(lhs >= rhs, name=name)
m.setObjective(gp.quicksum((scores[g][c] - pen) * y[g, c] for g in genes for c in comps), GRB.MAXIMIZE)
m.update()
m.write(str(X/"py_model.mps"))
print(f"wrote Python MPS: {m.NumVars} vars, {m.NumConstrs} constraints")
