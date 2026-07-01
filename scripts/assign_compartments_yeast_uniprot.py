"""End-to-end example: reproduce yeast-GEM's compartmentalisation from UniProt evidence.

This wires the whole pipeline together on a real fungal genome-scale model:

    UniProt subcellular-location annotations  (agnostic localisation evidence)
        -> raven_toolbox.localization.fetch_uniprot_localization   (gene x compartment scores)
        -> raven_toolbox.localization.assign_compartments                  (functionality-constrained MILP)
        -> a functional, compartmentalised model

It treats the curated compartment of each reaction in yeast-GEM as ground truth and asks: how
well does *agnostic* UniProt evidence, reconciled with network functionality, reproduce it — and
is the result a model that still grows? The MILP keeps biomass producible, so a reaction is placed
against its localisation score whenever functionality requires it (e.g. an essential mitochondrial
reaction stays mitochondrial even if UniProt annotates its gene as cytosolic).

Requirements
------------
* ``raven-toolbox`` with the modern localisation loaders (DeepLoc/MULocDeep/COMPARTMENTS/UniProt).
* A MILP solver (Gurobi/CPLEX/HiGHS) — GLPK is impractical at genome scale.
* Network access (the UniProt REST API) and the yeast-GEM SBML file.

Usage
-----
    python yeast_gem_uniprot.py --model /path/to/yeast-GEM.xml            # quick demo (subset)
    python yeast_gem_uniprot.py --model /path/to/yeast-GEM.xml --all      # full benchmark (slow)
"""
from __future__ import annotations

import argparse
import time
from collections import Counter

import cobra

from raven_toolbox.localization import (
    apply_assignment,
    assign_compartments,
    fetch_uniprot_localization,
)


def _base(met: cobra.Metabolite) -> str:
    # yeast-GEM keys the same species to a different id per compartment, so the
    # compartment-agnostic key is the metabolite name, not an id suffix.
    return met.name


def _sole_compartment(rxn: cobra.Reaction) -> str | None:
    comps = {m.compartment for m in rxn.metabolites}
    return next(iter(comps)) if len(comps) == 1 else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True, help="path to the yeast-GEM SBML file")
    ap.add_argument("--organism", type=int, default=559292,
                    help="UniProt organism/taxon id (default: 559292, S. cerevisiae S288C)")
    ap.add_argument("--all", action="store_true",
                    help="relocate every UniProt-covered reaction (slow); else a subset")
    ap.add_argument("--max-reactions", type=int, default=250,
                    help="cap on relocated reactions unless --all (default 250)")
    args = ap.parse_args()

    print(f"[1/4] fetching UniProt subcellular locations for organism {args.organism} ...")
    scores = fetch_uniprot_localization(args.organism)  # ordered-locus ids, mapped to model codes
    print(f"      {len(scores.genes)} genes annotated; compartments = {scores.compartments}")

    print(f"[2/4] loading model {args.model} ...")
    model = cobra.io.read_sbml_model(args.model)
    print(f"      {len(model.reactions)} reactions, baseline growth = {model.slim_optimize():.4f}")

    # Relocate internal, single-compartment reactions whose gene UniProt actually annotates, in a
    # compartment the evidence can speak to. Their curated compartment is the ground truth.
    scored, mappable = set(scores.genes), set(scores.compartments)
    movable = [r for r in model.reactions
               if not r.boundary and r.gene_reaction_rule
               and _sole_compartment(r) in mappable
               and any(g.id in scored for g in r.genes)]
    if not args.all:
        movable = movable[: args.max_reactions]
    truth = {r.id: _sole_compartment(r) for r in movable}
    print(f"[3/4] relocating {len(truth)} reactions "
          f"(curated split: {dict(Counter(truth.values()))})")

    t0 = time.time()
    # base_metabolite=name unifies species across compartments; transportable=[] keeps it tractable
    # at genome scale (no new transporters — pure relocation, biomass still enforced).
    res = assign_compartments(model, scores, list(truth),
                              base_metabolite=_base, transportable=[], time_limit=1200)
    print(f"      MILP status = {res.status}, solved in {time.time() - t0:.0f}s, "
          f"min_growth = {res.min_growth:.4f}")

    recovered = sum(1 for rid in truth if res.placements.get(rid) == [truth[rid]])
    out = apply_assignment(model, res, base_metabolite=_base)
    print("[4/4] results:")
    print(f"      agreement with curated yeast-GEM: {recovered}/{len(truth)} "
          f"= {100 * recovered / len(truth):.1f}%")
    print(f"      applied model still grows:         {out.slim_optimize():.4f} "
          f"(baseline {model.slim_optimize():.4f})")
    print("\nReactions placed against their UniProt score by the functionality constraint are the "
          "interesting cases — the network needs them where curation put them.")


if __name__ == "__main__":
    main()
