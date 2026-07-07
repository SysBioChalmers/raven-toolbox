#!/usr/bin/env python3
"""P1 scaling test: run the certified compartment assignment on Human-GEM (9 compartments, ~13k rxns).

The yeast run is 2296 reactions x 14 compartments; Human-GEM is ~6500 relocatable reactions x 9
compartments (a materially larger flux-free master), and neither its runtime nor whether it stays
certified is known past 4 compartments. This measures both, and reports Comparison-1 agreement against
Human-GEM's OWN compartment annotation (a circular truth — the P2 independent-truth work replaces it;
here the goal is scaling + certification, not the accuracy headline).

ASCII-only. Human-GEM genes are Ensembl ids, matching the committed DeepLoc predictions.
"""
from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import cobra  # noqa: E402
import pandas as pd  # noqa: E402
from cobra.flux_analysis import find_blocked_reactions  # noqa: E402

from raven_toolbox.localization import (  # noqa: E402
    LocalizationScores,
    apply_assignment,
    assign_compartments,
    load_deeploc,
)
from raven_toolbox.manipulation.compartments import merge_compartments  # noqa: E402

# DeepLoc label -> Human-GEM compartment code (note ER is 'r', lysosome 'l', peroxisome 'x'; Human-GEM
# has no plasma-membrane compartment so cell membrane maps to extracellular; Plastid is dropped).
HUMANGEM_MAP = {
    "cytoplasm": "c", "cytosol": "c", "nucleus": "n", "nucleoplasm": "n",
    "mitochondrion": "m", "mitochondria": "m", "mitochondrial": "m",
    "peroxisome": "x", "endoplasmic reticulum": "r", "golgi apparatus": "g", "golgi": "g",
    "lysosome/vacuole": "l", "lysosome": "l",
    "extracellular": "e", "extracellular space": "e", "cell membrane": "e", "plasma membrane": "e",
}
_PARENT = {"i": "m"}  # inner mitochondria -> mitochondrion (DeepLoc cannot resolve it)


def _name(m):
    return m.name


def _norm(c):
    return _PARENT.get(c, c)


def _sole(r):
    comps = {x.compartment for x in r.metabolites}
    return next(iter(comps)) if len(comps) == 1 else None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", type=Path, default=Path("C:/Work/GitHub/Human-GEM/model/Human-GEM.yml"))
    ap.add_argument("--data-dir", type=Path, default=Path("data/deeploc/humangem"))
    ap.add_argument("--min-growth-fraction", type=float, default=0.5)
    ap.add_argument("--time-limit", type=float, default=900.0, help="master solve budget (s)")
    ap.add_argument("--prune-transports", action="store_true", help="also FVA-prune (slow; off for the scaling test)")
    ap.add_argument("--out", type=Path, default=Path(".research_tmp/certified_humangem.json"))
    args = ap.parse_args(argv)

    model = cobra.io.load_yaml_model(str(args.model))
    biomass_id = next(r.id for r in model.reactions if r.objective_coefficient != 0)
    curated_growth = model.slim_optimize()
    curated_rxn = {r.id: _norm(c) for r in model.reactions
                   if not r.boundary and (c := _sole(r)) is not None}
    curated_gene = {}
    for rid, c in curated_rxn.items():
        for g in model.reactions.get_by_id(rid).genes:
            curated_gene.setdefault(g.id, set()).add(c)
    print(f"Human-GEM: {len(model.reactions)} rxns, {len(model.genes)} genes, "
          f"{len(model.compartments)} compartments, growth {curated_growth:.2f}", flush=True)

    draft, _d, _u = merge_compartments(model, merged_id="c", merged_name="cytosol",
                                       base_metabolite=_name, drop_single_metabolite_reactions=False)
    draft.objective = biomass_id
    frames = [load_deeploc(c, compartment_map=HUMANGEM_MAP).df
              for c in sorted(args.data_dir.glob("Human-GEM_deeploc_*.csv"))]
    scores = LocalizationScores(pd.concat(frames))
    relocate = [r.id for r in draft.reactions if not r.boundary and r.id != biomass_id]
    min_growth = args.min_growth_fraction * curated_growth
    print(f"flattened draft: {len(draft.reactions)} rxns; scores {len(scores.df)} genes x "
          f"{list(scores.df.columns)}; relocate {len(relocate)}; min_growth {min_growth:.2f}", flush=True)

    t0 = time.monotonic()
    prop = assign_compartments(
        draft, scores, relocate, default_compartment="c", base_metabolite=_name,
        biomass_reaction=biomass_id, min_growth=min_growth, time_limit=args.time_limit,
        prune_transports=args.prune_transports)
    wall = round(time.monotonic() - t0, 1)
    print(f"certified: {wall}s  status={prop.status}  certified={prop.certified}", flush=True)

    out = {"organism": "Human-GEM", "wall_s": wall, "status": prop.status, "certified": prop.certified,
           "compartments": len(model.compartments), "relocate": len(relocate),
           "curated_growth": round(curated_growth, 3), "min_growth": round(min_growth, 3)}
    if prop.placements:
        res_comp = {rid: _norm(cs[0]) for rid, cs in prop.placements.items() if cs}
        common = set(res_comp) & set(curated_rxn)
        out["reaction_agreement"] = round(sum(res_comp[r] == curated_rxn[r] for r in common) / len(common), 4)
        out["reaction_n"] = len(common)
        res_gene = {g: {_norm(c) for c in cs} for g, cs in prop.gene_compartments.items() if cs}
        common_g = set(res_gene) & set(curated_gene)
        out["gene_agreement"] = round(sum(bool(res_gene[g] & curated_gene[g]) for g in common_g) / len(common_g), 4)
        out["gene_n"] = len(common_g)
        out["transports"] = len(prop.added_transports)
        applied = apply_assignment(draft, prop, default_compartment="c", base_metabolite=_name)
        applied.objective = biomass_id
        out["growth"] = round(applied.slim_optimize(error_value=0.0) or 0.0, 3)
        blocked = set(find_blocked_reactions(applied)) & set(relocate)
        out["blocked_fraction"] = round(len(blocked) / len(relocate), 4)
        print(f"  reaction agreement {out['reaction_agreement']:.1%} ({out['reaction_n']}); "
              f"gene agreement {out['gene_agreement']:.1%} ({out['gene_n']})", flush=True)
        print(f"  transports {out['transports']}; growth {out['growth']}; "
              f"blocked {out['blocked_fraction']:.1%}", flush=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=str))
    print(f"written -> {args.out}  |  DONE", flush=True)


if __name__ == "__main__":
    main()
