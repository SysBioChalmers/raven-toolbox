#!/usr/bin/env python3
"""Benchmark our evidence-aware transport reduction against a native reference fungal reconstruction.

The reference tool applies **no** transport minimisation: inter-compartment transport reactions get a
near-zero score (~1e-11) in its scoring and no penalty in its carve MILP, so its network keeps every
transport its genome-scale gene content can support -- functional, but bloated with unsupported
transports. Our evidence-aware transport cost reduces that network *selectively*: keep a transport when
a transporter gene supports its cargo (right substrate class / ChEBI, right membrane) OR the network
needs it to stay feasible; drop it only when neither holds.

Takes the reference tool's OWN reconstructed, gene-annotated, functional model (built by
``build_reference_carve_model.py`` -- NOT a bare reaction-id cache, which is lossy: it drops the tool's
uptake reactions and each reaction's solved direction, so it cannot grow standalone and cannot support a
real feasibility check). On that model:

* **native** -- every inter-compartment transport the tool kept (baseline; no evidence used).
* **ours** -- a greedy, feasibility-respecting reduction: rank transports whose best cargo has no
  transporter support (cost >= ``--keep-threshold``) by cost, worst-evidenced first; tentatively drop
  each one and re-run FBA; keep it dropped if growth survives, otherwise put it back (feasibility
  overrides the missing evidence). This is one-at-a-time and greedy, not a joint MILP re-solve, so it is
  an upper bound on what feasibility alone would allow to drop -- a joint solve could only remove the
  same or fewer.

Both are scored against the curated yeast-GEM transportome (ground truth): curated transports
replicated, individually-essential ones kept, and spurious (non-curated) ones carried. The sibling sweep
shows the parsimony/completeness trade-off. yeast proteome + DeepLoc ship under data/deeploc/.
ASCII-only. cobra + pandas required.
"""
from __future__ import annotations

import argparse
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
GROWTH_FRACTION = 0.01  # a candidate transport is "feasibility-needed" if dropping it costs > 99% growth


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


def _carve_to_yeast(rxn_of, edge2rxn, rid):
    out = set()
    for e in rxn_edges(rxn_of[rid]):
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


def _essential(Y, rxn_of, edge2rxn, transport_ids):
    """Carve transports mapping to a yeast-GEM transport whose single deletion abolishes growth."""
    g0 = Y.slim_optimize()
    allyr = set()
    for rid in transport_ids:
        allyr |= _carve_to_yeast(rxn_of, edge2rxn, rid)
    ess_yeast = set()
    for x in sorted(allyr):
        with Y:
            Y.remove_reactions([Y.reactions.get_by_id(x)], remove_orphans=False)
            if (Y.slim_optimize() or 0.0) < GROWTH_FRACTION * g0:
                ess_yeast.add(x)
    return {rid for rid in transport_ids if _carve_to_yeast(rxn_of, edge2rxn, rid) & ess_yeast}


def best_cargo_cost(reaction, cost, base_cost):
    """The most-evidenced cargo cost for a transport (no cargo scored -> the full, unsupported prior)."""
    cargo = [c for c in (shuttled(reaction) - TRIVIAL) if c in cost]
    return min((cost[c] for c in cargo), default=base_cost)


def feasibility_respecting_reduction(model, transport_ids, cost, keep_threshold, base_cost):
    """Greedily drop unsupported transports (cost >= keep_threshold) worst-evidence-first, keeping any
    whose removal breaks growth. Evidence-supported transports (cost < keep_threshold) are never
    candidates -- feasibility can only ADD transports back, never remove a supported one. One-at-a-time
    on a private copy of ``model``, so this is an upper bound on what a joint MILP would drop (a joint
    solve sees the same trade-offs at once and could only be more conservative, i.e. keep >= as many).

    Returns (kept_ids, dropped_ids, feasibility_forced_ids, achieved_growth).

    Uses bounds=(0, 0) knockout rather than physical ``remove_reactions``/``add_reactions``: cobra's
    remove/re-add round trip does not reliably restore the solver's internal LP structure (the model
    can come back permanently infeasible even after "putting back" a reaction), whereas a bounds-based
    knockout -- the same technique ``single_reaction_deletion`` uses -- is exactly reversible.
    """
    work = model.copy()
    g0 = work.slim_optimize()
    # Many transports share the exact same fallback cost (unmatched cargo -> base_cost verbatim), so
    # pre-sort by id for a deterministic tie-break: otherwise the removal order among tied candidates
    # follows Python set-iteration order, which is reproducible per-run but not a meaningful signal, and
    # would make results vary across otherwise-identical runs / obscure genuine cross-variant deltas.
    candidates = sorted(rid for rid in transport_ids
                        if best_cargo_cost(work.reactions.get_by_id(rid), cost, base_cost) >= keep_threshold)
    candidates.sort(key=lambda rid: best_cargo_cost(work.reactions.get_by_id(rid), cost, base_cost),
                    reverse=True)  # worst-evidenced (highest cost) first; stable sort preserves id order for ties
    dropped, forced = [], []
    for rid in candidates:
        rxn = work.reactions.get_by_id(rid)
        lb, ub = rxn.bounds
        rxn.bounds = (0, 0)
        g = work.slim_optimize()
        if not g or g < GROWTH_FRACTION * g0:
            rxn.bounds = (lb, ub)  # feasibility requires it -- restore despite lacking evidence
            forced.append(rid)
        else:
            dropped.append(rid)  # leave blocked
    kept = set(transport_ids) - set(dropped)
    return kept, set(dropped), set(forced), work.slim_optimize()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True, type=Path,
                    help="the reference tool's own reconstructed, gene-annotated, functional SBML "
                        "(from build_reference_carve_model.py) -- NOT a bare reaction-id cache")
    ap.add_argument("--yeast-gem", required=True, type=Path)
    ap.add_argument("--data-dir", type=Path, default=Path("data/deeploc"))
    ap.add_argument("--keep-threshold", type=float, default=0.35,
                    help="a transport is a drop-candidate when its best cargo cost is >= this")
    ap.add_argument("--base-cost", type=float, default=0.5)
    args = ap.parse_args()

    CF = cobra.io.read_sbml_model(str(args.model))
    Y = cobra.io.read_sbml_model(str(args.yeast_gem))
    rxn_of = {r.id: r for r in CF.reactions}
    transport_ids = {r.id for r in CF.reactions if is_transport(r)}
    edge2rxn, y_pairs = curated_edges(Y)

    g0 = CF.slim_optimize()
    print(f"reference reconstruction: {len(CF.reactions)} reactions, {len(CF.genes)} genes, "
          f"growth = {g0:.4f}")
    print(f"  {len(transport_ids)} inter-compartment transports")

    ann, gene_comps, onto, n_ann = annotate_yeast(args.data_dir)
    curated = {rid for rid in transport_ids if any(e in y_pairs for e in rxn_edges(rxn_of[rid]))}
    essential = _essential(Y, rxn_of, edge2rxn, transport_ids)
    nc, ne = len(curated), len(essential)
    print(f"  {nc} match curated yeast-GEM transports, {ne} individually essential "
          f"(evidence from {n_ann} annotated yeast transporter genes)")

    print("\n== benchmark: replicating the curated yeast-GEM transportome (feasibility-respecting) ==")
    print(f"   {'approach':30s} {'kept':>4} {'curated repl':>13} {'essential':>10} {'spurious':>9} "
          f"{'growth':>8}")

    def row(name, kept, growth):
        print(f"   {name:30s} {len(kept):4d} {len(kept & curated):8d}/{nc:<4d} "
              f"{len(kept & essential):6d}/{ne:<3d} {len(kept - curated):9d} {growth:8.4f}")

    row("reference (native)", transport_ids, g0)

    def variant(label, cost):
        kept, dropped, forced, g = feasibility_respecting_reduction(
            CF, transport_ids, cost, args.keep_threshold, args.base_cost)
        row(label, kept, g)
        return kept, dropped, forced

    variant("ours: coarse", _cost(Y, ann, gene_comps, None, args.base_cost, 0.0))
    for sw in SIBLING_SWEEP:
        cost = _cost(Y, ann, gene_comps, onto, args.base_cost, sw)
        label = "ours: +ChEBI" if sw == 0.0 else f"ours: +ChEBI +sibling {sw:g}"
        variant(label, cost)


if __name__ == "__main__":
    main()
