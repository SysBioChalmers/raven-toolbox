#!/usr/bin/env python3
"""Is less better? Dig into the transports Arm A keeps vs Arm B drops in the CarveFungi head-to-head.

Three validations of the transport-parsimony effect, all from the cached MILP solutions (no re-solve):
1. CURATED MATCH  -- do kept/dropped transports correspond to a curated yeast-GEM (literature-backed)
   transport? If the reduction were "smart", dropped transports would match LESS than kept ones.
2. FUNCTIONAL     -- remove the dropped-but-curated transports from yeast-GEM (a proper model with
   genes + validated growth) and measure the growth / individual-essentiality impact.
3. CONNECTIVITY   -- structural dead-end metabolites and isolated sub-networks (the carve guarantees
   flux-connectivity WITH exchanges; this exposes the latent gaps exchanges paper over).

Inputs mirror run_carvefungi_cplex.py: needs the CarveFungi clone (universal model), yeast-GEM, and the
per-arm solution cache that run_carvefungi_cplex.py writes (arm_A.json / arm_B.json, each with an
"active" reaction-id list). ASCII-only output. cobra + pandas required.
"""
from __future__ import annotations

import argparse
import json
import re
import warnings
from collections import Counter, defaultdict, deque
from pathlib import Path

warnings.filterwarnings("ignore")

import cobra  # noqa: E402
import pandas as pd  # noqa: E402

from raven_toolbox.localization import (  # noqa: E402
    DEFAULT_COMPARTMENT_MAP,
    SubstrateOntology,
    annotate_proteome,
    default_substrate_of,
    evidence_aware_transport_cost,
    load_deeploc,
)

UNIV2YEAST = {"c": "c", "m": "m", "x": "p", "r": "er", "n": "n", "g": "g", "e": "e", "l": "lp"}
CNAME = {"c": "cytosol", "m": "mito", "x": "peroxisome", "r": "ER", "n": "nucleus",
         "g": "Golgi", "e": "extracell", "l": "lipid"}
TRIVIAL = {"h+", "water", "dioxygen", "carbon dioxide"}  # transported ~everywhere; not discriminating


def norm(name: str) -> str:
    s = re.sub(r"\(\d*[+-]\)", "", name.lower().strip()).replace(" zwitterion", "").replace("'", "")
    return re.sub(r"\s+", " ", s).strip()


def shuttled(reaction):
    """metabolite NAMES that appear in >1 compartment within the reaction (the transported species)."""
    by = defaultdict(set)
    for m in reaction.metabolites:
        by[norm(m.name)].add(m.compartment)
    return {n for n, cs in by.items() if len(cs) > 1}


def is_transport(reaction):
    return len({m.compartment for m in reaction.metabolites}) > 1


def transport_section(U, rxnU, actA, actB):
    A = {r for r in actA if r in rxnU and is_transport(rxnU[r])}
    B = {r for r in actB if r in rxnU and is_transport(rxnU[r])}
    shared, a_only, b_only = A & B, A - B, B - A
    print("== transports ==")
    print(f"Arm A={len(A)}  Arm B={len(B)}  shared={len(shared)}  "
          f"A_only(dropped by B)={len(a_only)}  B_only={len(b_only)}")
    pair = Counter(tuple(sorted({m.compartment for m in rxnU[r].metabolites})) for r in a_only)
    print("dropped-by-B by compartment pair: "
          + ", ".join(f"{'-'.join(CNAME.get(c, c) for c in p)}:{n}" for p, n in pair.most_common()))
    cargo = Counter(n for r in a_only for n in (shuttled(rxnU[r]) - TRIVIAL))
    print("most-dropped cargo: " + ", ".join(f"{n}({c})" for n, c in cargo.most_common(8)))
    return A, B, shared, a_only, b_only


def curated_edges(Y):
    edge2rxn = defaultdict(set)
    for r in Y.reactions:
        if r.boundary:
            continue
        for n in shuttled(r):
            comps = frozenset(m.compartment for m in r.metabolites if norm(m.name) == n)
            edge2rxn[(n, comps)].add(r.id)
    pairs = set(edge2rxn)
    return edge2rxn, pairs


def rxn_edges(reaction):
    """(norm metabolite, yeast-mapped compartment pair) edges this carved transport realises."""
    ypair = frozenset(UNIV2YEAST.get(c, c) for c in {m.compartment for m in reaction.metabolites})
    return {(n, ypair) for n in (shuttled(reaction) - TRIVIAL)}


def curated_match_section(rxnU, y_pairs, sets):
    print("\n== curated-yeast-GEM match (strict = metabolite + same compartment pair) ==")
    for label, ids in sets:
        classifiable = matched = 0
        for r in ids:
            es = rxn_edges(rxnU[r])
            if not es:
                continue
            classifiable += 1
            matched += any(e in y_pairs for e in es)
        rate = matched / classifiable if classifiable else float("nan")
        print(f"  {label:26s} classifiable={classifiable:3d}  match curated={matched:3d} ({rate:.0%})")


def functional_section(Y, rxnU, edge2rxn, a_only, shared):
    print("\n== functional impact on yeast-GEM (proper model: genes + validated growth) ==")

    def y_rxns(ids):
        s = set()
        for r in ids:
            for e in rxn_edges(rxnU[r]):
                s |= edge2rxn.get(e, set())
        return s

    g0 = Y.slim_optimize()
    print(f"yeast-GEM baseline growth: {g0:.4f}")

    def growth_without(rxn_ids):
        with Y:
            Y.remove_reactions([Y.reactions.get_by_id(x) for x in rxn_ids if Y.reactions.has_id(x)],
                               remove_orphans=False)
            return Y.slim_optimize() or 0.0

    dropped_essential: list[str] = []
    for label, ids in (("dropped-by-B", a_only), ("shared", shared)):
        yr = y_rxns(ids)
        ess = [x for x in sorted(yr) if growth_without({x}) < 0.01 * g0]
        print(f"  {label:14s}: maps to {len(yr):2d} curated transports; remove all -> growth "
              f"{growth_without(yr):.4f}; individually essential: {len(ess)}")
        if label == "dropped-by-B":
            dropped_essential = ess
            for x in ess:
                print(f"      essential: {x:10s} {Y.reactions.get_by_id(x).name}")
    return dropped_essential


def connectivity_section(rxnU, E_BND, actA, actB, U):
    print("\n== structural connectivity (internal network; exchanges excluded) ==")

    def analyse(active):
        rxns = [rxnU[r] for r in active if r in rxnU and r not in E_BND]
        prod, cons, mets, adj = set(), set(), set(), defaultdict(set)
        for r in rxns:
            rev = r.lower_bound < 0
            for m, c in r.metabolites.items():
                mets.add(m.id)
                adj[("r", r.id)].add(("m", m.id))
                adj[("m", m.id)].add(("r", r.id))
                if (r.upper_bound > 0 and c > 0) or (rev and c != 0):
                    prod.add(m.id)
                if (r.upper_bound > 0 and c < 0) or (rev and c != 0):
                    cons.add(m.id)
        deadend = {m for m in mets if not (m in prod and m in cons)}
        seen, comps = set(), []
        for node in adj:
            if node in seen:
                continue
            q, n_r = deque([node]), 0
            seen.add(node)
            while q:
                x = q.popleft()
                n_r += x[0] == "r"
                for y in adj[x]:
                    if y not in seen:
                        seen.add(y)
                        q.append(y)
            comps.append(n_r)
        return mets, deadend, sorted(comps, reverse=True)

    mA, dA, cA = analyse(actA)
    mB, dB, cB = analyse(actB)
    print(f"  Arm A: {len(mA)} mets; dead-ends {len(dA)} ({len(dA)/len(mA):.1%}); "
          f"components {len(cA)} (largest {cA[0]} rxns, rest {sum(cA[1:])})")
    print(f"  Arm B: {len(mB)} mets; dead-ends {len(dB)} ({len(dB)/len(mB):.1%}); "
          f"components {len(cB)} (largest {cB[0]} rxns, rest {sum(cB[1:])})")
    both = mA & mB
    mid2c = {m.id: m.compartment for m in U.metabolites}
    stranded = (dB - dA) & both
    rev = (dA - dB) & both
    bycomp = Counter(CNAME.get(mid2c.get(m), "?") for m in stranded)
    print(f"  balanced in A but dead-end in B: {len(stranded)} (reverse: {len(rev)}); "
          f"by compartment: " + ", ".join(f"{c}:{n}" for c, n in bycomp.most_common()))


def evidence_cost_on_yeast(Y, data_dir, sibling_weight, base_cost):
    """Per-metabolite (norm-name-keyed) evidence-aware transport cost, from the yeast proteome."""
    import tempfile

    fastas = sorted(data_dir.glob("*_proteins_*.fasta"))
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as tf:
        tf.write("".join(p.read_text() for p in fastas))
        proteome = Path(tf.name)
    ann = annotate_proteome(proteome, threads=4)
    proteome.unlink(missing_ok=True)
    sdf = pd.concat([load_deeploc(c, compartment_map=DEFAULT_COMPARTMENT_MAP).df
                     for c in sorted(data_dir.glob("*_deeploc_*.csv"))])
    gene_comps = {g: {sdf.loc[g].astype(float).idxmax()} for g in sdf.index if g in ann}
    onto = SubstrateOntology.load()
    cost = evidence_aware_transport_cost(Y, ann, gene_comps, substrate_of=default_substrate_of,
                                         ontology=onto, sibling_weight=sibling_weight,
                                         base_cost=base_cost, base_metabolite=lambda m: norm(m.name))
    return cost, len(ann)


def evidence_section(Y, rxnU, y_pairs, cost, base_cost, sets, dropped_essential):
    """Would our evidence-aware cost *rescue* the transports the blanket penalty (Arm B) drops?

    For each transport, look up its cargo's evidence-aware cost (from the yeast proteome). A transport
    is 'kept' when any cargo scores below ``base_cost`` (has transporter support). The blanket penalty
    drops all of Arm A's extra transports regardless; a *selective* cost keeps the curated ones and the
    spurious ones stay expensive."""
    print("\n== evidence-aware rescue of the dropped transports (would our cost keep them?) ==")

    def summarise(label, ids):
        rows = []  # (is_curated, would_keep, min_cargo_cost)
        for r in ids:
            cargo = [c for c in (shuttled(rxnU[r]) - TRIVIAL) if c in cost]
            if not cargo:
                continue  # cargo absent from yeast-GEM -> not comparable
            mn = min(cost[c] for c in cargo)
            rows.append((any(e in y_pairs for e in rxn_edges(rxnU[r])), mn < base_cost, mn))
        for tag, sub in (("curated", [x for x in rows if x[0]]),
                         ("non-curated", [x for x in rows if not x[0]])):
            if sub:
                keep = sum(k for _, k, _ in sub) / len(sub)
                mc = sum(c for _, _, c in sub) / len(sub)
                print(f"  {label:16s} {tag:12s} n={len(sub):3d}  would-keep {keep:3.0%}  mean cost {mc:.2f}")

    for label, ids in sets:
        summarise(label, ids)

    print("  -- the individually-essential transports the blanket cut drops (cargo cost; want < "
          f"{base_cost}) --")
    for x in dropped_essential:
        rx = Y.reactions.get_by_id(x)
        cargo = {c: round(cost[c], 2) for c in sorted(shuttled(rx) - TRIVIAL) if c in cost}
        rescued = cargo and min(cargo.values()) < base_cost
        print(f"    {x:9s} {rx.name[:40]:40s} {cargo}  [{'RESCUED' if rescued else 'not evidenced'}]")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--carvefungi-dir", required=True, type=Path)
    ap.add_argument("--yeast-gem", required=True, type=Path)
    ap.add_argument("--cache-dir", required=True, type=Path,
                    help="dir with arm_A.json / arm_B.json from run_carvefungi_cplex.py")
    ap.add_argument("--universal", type=Path)
    ap.add_argument("--data-dir", type=Path, default=Path("data/deeploc"),
                    help="proteome + DeepLoc for the evidence-aware rescue (default: data/deeploc)")
    ap.add_argument("--sibling-weight", type=float, default=0.0)
    ap.add_argument("--base-cost", type=float, default=0.5)
    args = ap.parse_args()

    universal = args.universal or args.carvefungi_dir / "data" / "reactionDatabase" / "bigModelv2.21b.sbml"
    U = cobra.io.read_sbml_model(str(universal))
    Y = cobra.io.read_sbml_model(str(args.yeast_gem))
    rxnU = {r.id: r for r in U.reactions}
    E_BND = {r.id for r in U.reactions if r.id.endswith("_E") or r.boundary}
    actA = set(json.loads((args.cache_dir / "arm_A.json").read_text())["active"])
    actB = set(json.loads((args.cache_dir / "arm_B.json").read_text())["active"])

    A, B, shared, a_only, b_only = transport_section(U, rxnU, actA, actB)
    _, y_pairs = curated_edges(Y)
    curated_match_section(rxnU, y_pairs, [("shared (both keep)", shared),
                                         ("dropped by Arm B", a_only),
                                         ("kept by Arm A (all)", A)])
    edge2rxn, _ = curated_edges(Y)
    dropped_ess = functional_section(Y, rxnU, edge2rxn, a_only, shared)
    connectivity_section(rxnU, E_BND, actA, actB, U)

    cost, n_ann = evidence_cost_on_yeast(Y, args.data_dir, args.sibling_weight, args.base_cost)
    print(f"\n(evidence-aware cost from {n_ann} annotated yeast transporter genes; "
          f"sibling_weight={args.sibling_weight})")
    evidence_section(Y, rxnU, y_pairs, cost, args.base_cost,
                     [("shared (kept)", shared), ("dropped by Arm B", a_only)], dropped_ess)


if __name__ == "__main__":
    main()
