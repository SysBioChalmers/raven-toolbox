#!/usr/bin/env python3
"""Benchmark our evidence-aware transport reduction against native CarveFungi, on CarveFungi's own carve.

CarveFungi applies **no** transport minimisation: inter-compartment transport reactions get a near-zero
score (~1e-11) in its scoring and no penalty in its carve MILP, so its network keeps every transport the
reactions can support -- functional, but bloated with unsupported transports. Our evidence-aware
transport cost reduces that network *selectively*: keep a transport when a transporter gene supports its
cargo (right substrate class / ChEBI, right membrane), drop it when nothing does.

Both are scored on the SAME candidate set -- CarveFungi's own carved transports (``arm_A.json`` from
``run_carvefungi_cplex.py``, its ``minmax_reduction`` run unmodified) -- against the curated yeast-GEM
transportome (ground truth): curated transports replicated, individually-essential ones kept, and
spurious (non-curated) ones carried. ``ours`` applies our per-metabolite cost independently (keep a
transport when its best cargo cost is below ``--keep-threshold``); it is not a re-solve of the carve
MILP. The last rows sweep ``sibling_weight`` -- crediting chemical *relatives* of a curated substrate --
to show the parsimony/completeness trade-off. yeast proteome + DeepLoc ship under data/deeploc/.
ASCII-only. cobra + pandas required.
"""
from __future__ import annotations

import argparse
import json
import re
import tempfile
import warnings
from collections import defaultdict
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
TRIVIAL = {"h+", "water", "dioxygen", "carbon dioxide"}  # transported ~everywhere; not discriminating
SIBLING_SWEEP = (0.0, 0.3, 0.5, 0.7, 1.0)


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


def curated_edges(Y):
    """(norm metabolite, compartment-pair) -> yeast-GEM transport reactions realising it."""
    edge2rxn = defaultdict(set)
    for r in Y.reactions:
        if r.boundary:
            continue
        for n in shuttled(r):
            comps = frozenset(m.compartment for m in r.metabolites if norm(m.name) == n)
            edge2rxn[(n, comps)].add(r.id)
    return edge2rxn, set(edge2rxn)


def rxn_edges(reaction):
    """(norm metabolite, yeast-mapped compartment pair) edges this carved transport realises."""
    ypair = frozenset(UNIV2YEAST.get(c, c) for c in {m.compartment for m in reaction.metabolites})
    return {(n, ypair) for n in (shuttled(reaction) - TRIVIAL)}


def _carve_to_yeast(rxnU, edge2rxn, r):
    out = set()
    for e in rxn_edges(rxnU[r]):
        out |= edge2rxn.get(e, set())
    return out


def annotate_yeast(data_dir):
    """Annotate the yeast proteome once; return (annotation, gene->compartments, ontology, n_genes)."""
    fastas = sorted(data_dir.glob("*_proteins_*.fasta"))
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as tf:
        tf.write("".join(p.read_text() for p in fastas))
        proteome = Path(tf.name)
    ann = annotate_proteome(proteome, threads=4)
    proteome.unlink(missing_ok=True)
    sdf = pd.concat([load_deeploc(c, compartment_map=DEFAULT_COMPARTMENT_MAP).df
                     for c in sorted(data_dir.glob("*_deeploc_*.csv"))])
    gene_comps = {g: {sdf.loc[g].astype(float).idxmax()} for g in sdf.index if g in ann}
    return ann, gene_comps, SubstrateOntology.load(), len(ann)


def _cost(Y, ann, gene_comps, onto, base_cost, sibling_weight):
    return evidence_aware_transport_cost(Y, ann, gene_comps, substrate_of=default_substrate_of,
                                         ontology=onto, sibling_weight=sibling_weight,
                                         base_cost=base_cost, base_metabolite=lambda m: norm(m.name))


def _essential(Y, rxnU, edge2rxn, A):
    """Carve transports mapping to a yeast-GEM transport whose single deletion abolishes growth."""
    g0 = Y.slim_optimize()
    allyr = set()
    for r in A:
        allyr |= _carve_to_yeast(rxnU, edge2rxn, r)
    ess_yeast = set()
    for x in sorted(allyr):
        with Y:
            Y.remove_reactions([Y.reactions.get_by_id(x)], remove_orphans=False)
            if (Y.slim_optimize() or 0.0) < 0.01 * g0:
                ess_yeast.add(x)
    return {r for r in A if _carve_to_yeast(rxnU, edge2rxn, r) & ess_yeast}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--carvefungi-dir", required=True, type=Path)
    ap.add_argument("--yeast-gem", required=True, type=Path)
    ap.add_argument("--cache-dir", required=True, type=Path,
                    help="dir with arm_A.json (the native carve, from run_carvefungi_cplex.py)")
    ap.add_argument("--universal", type=Path)
    ap.add_argument("--data-dir", type=Path, default=Path("data/deeploc"))
    ap.add_argument("--keep-threshold", type=float, default=0.35,
                    help="ours keeps a transport when its best cargo cost is below this (< base-cost)")
    ap.add_argument("--base-cost", type=float, default=0.5)
    args = ap.parse_args()

    universal = args.universal or args.carvefungi_dir / "data" / "reactionDatabase" / "bigModelv2.21b.sbml"
    U = cobra.io.read_sbml_model(str(universal))
    Y = cobra.io.read_sbml_model(str(args.yeast_gem))
    rxnU = {r.id: r for r in U.reactions}
    active = set(json.loads((args.cache_dir / "arm_A.json").read_text())["active"])  # native carve
    A = {r for r in active if r in rxnU and is_transport(rxnU[r])}
    edge2rxn, y_pairs = curated_edges(Y)

    ann, gene_comps, onto, n_ann = annotate_yeast(args.data_dir)
    curated = {r for r in A if any(e in y_pairs for e in rxn_edges(rxnU[r]))}
    essential = _essential(Y, rxnU, edge2rxn, A)
    nc, ne = len(curated), len(essential)
    print(f"native CarveFungi carve: {len(A)} inter-compartment transports "
          f"({nc} match curated yeast-GEM, {ne} individually essential); "
          f"evidence from {n_ann} annotated yeast transporter genes; keep-threshold {args.keep_threshold}")

    def evidenced(r, cost):
        cargo = [c for c in (shuttled(rxnU[r]) - TRIVIAL) if c in cost]
        return bool(cargo) and min(cost[c] for c in cargo) < args.keep_threshold

    print("\n== benchmark: replicating the curated yeast-GEM transportome ==")
    print(f"   {'approach':28s} {'kept':>4} {'curated repl':>13} {'essential':>10} {'spurious':>9}")

    def row(name, kept):
        print(f"   {name:28s} {len(kept):4d} {len(kept & curated):8d}/{nc:<4d} "
              f"{len(kept & essential):6d}/{ne:<3d} {len(kept - curated):9d}")

    row("CarveFungi (native)", set(A))
    row("ours: coarse", {r for r in A if evidenced(r, _cost(Y, ann, gene_comps, None, args.base_cost, 0.0))})
    for sw in SIBLING_SWEEP:
        cost = _cost(Y, ann, gene_comps, onto, args.base_cost, sw)
        label = "ours: +ChEBI" if sw == 0.0 else f"ours: +ChEBI +sibling {sw:g}"
        row(label, {r for r in A if evidenced(r, cost)})


if __name__ == "__main__":
    main()
