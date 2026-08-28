#!/usr/bin/env python3
"""Extract per-compartment structure for the biology comparison (yeast-GEM / certified / CarveFungi).

Materialises the certified placement, then for each of the three models writes:
  * transporters_per_compartment  (inter-compartment reactions touching each compartment)
  * reactions_per_compartment     (single-compartment reactions per compartment)
  * gene_compartments             (gene -> sorted list of compartments of its single-compartment rxns)
  * reaction_compartments         (single-compartment reaction id -> compartment)
to .research_tmp/compartment_structure.json, for downstream pathway / dual-localisation scoring.

ASCII-only.
"""
from __future__ import annotations

import json
import sys
import time
import warnings
from collections import Counter
from pathlib import Path

warnings.filterwarnings("ignore")
import cobra  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_replicate_yeast_gem import _name, build_draft, load_yeast_scores  # noqa: E402

from raven_toolbox.localization import (  # noqa: E402
    apply_assignment,
    assign_compartments,
)


def _sole(rxn):
    comps = {m.compartment for m in rxn.metabolites}
    return next(iter(comps)) if len(comps) == 1 else None


def structure(model, biomass_id=None):
    tr_per = Counter()
    rxn_per = Counter()
    rxn_comp = {}
    n_transports = 0
    for r in model.reactions:
        if r.boundary:
            continue
        comps = {m.compartment for m in r.metabolites if m.compartment}
        if len(comps) == 1:
            c = next(iter(comps))
            rxn_per[c] += 1
            rxn_comp[r.id] = c
        elif len(comps) > 1:
            n_transports += 1
            for c in comps:
                tr_per[c] += 1
    gene_comp = {}
    for rid, c in rxn_comp.items():
        for g in model.reactions.get_by_id(rid).genes:
            gene_comp.setdefault(g.id, set()).add(c)
    return {
        "transporters_per_compartment": dict(tr_per),
        "reactions_per_compartment": dict(rxn_per),
        "total_transports": n_transports,
        "gene_compartments": {g: sorted(cs) for g, cs in gene_comp.items()},
        "reaction_compartments": rxn_comp,
        "compartments": sorted({m.compartment for m in model.metabolites if m.compartment}),
    }


def main():
    yeast = cobra.io.read_sbml_model("C:/Work/GitHub/yeast-GEM/model/yeast-GEM.xml")
    cg = yeast.slim_optimize()
    draft, biomass_id = build_draft(yeast)
    scores = load_yeast_scores(Path("data/deeploc"))
    relocate = [r.id for r in draft.reactions if not r.boundary and r.id != biomass_id]

    print("running certified (prune=True) to materialise the transports ...", flush=True)
    t0 = time.monotonic()
    prop = assign_compartments(
        draft, scores, relocate, default_compartment="c", base_metabolite=_name,
        biomass_reaction=biomass_id, min_growth=0.5 * cg)
    certified = apply_assignment(draft, prop, default_compartment="c", base_metabolite=_name)
    print(f"  certified materialised in {time.monotonic() - t0:.0f}s; "
          f"transports={len(prop.added_transports)}", flush=True)
    cobra.io.write_sbml_model(certified, ".research_tmp/certified_yeast_model.xml")

    cf = cobra.io.read_sbml_model(".research_tmp/carvefungi_yeast_model.sbml")

    out = {
        "yeast_gem": structure(yeast),
        "certified": structure(certified),
        "carvefungi": structure(cf),
        "certified_gene_compartments_from_proposal":
            {g: cs for g, cs in prop.gene_compartments.items() if cs},
    }
    Path(".research_tmp/compartment_structure.json").write_text(json.dumps(out, default=str))
    for name in ("yeast_gem", "certified", "carvefungi"):
        s = out[name]
        print(f"\n{name}: {s['total_transports']} transports; compartments {s['compartments']}")
        print("  transporters/compartment:", s["transporters_per_compartment"])
        print("  reactions/compartment:", s["reactions_per_compartment"])
    print("\nwritten -> .research_tmp/compartment_structure.json")


if __name__ == "__main__":
    main()
