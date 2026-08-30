#!/usr/bin/env python3
"""Solve ONLY the placement master once and hash the placement. Run under two
PYTHONHASHSEED values: identical hashes prove the master no longer depends on
set-iteration order."""
from __future__ import annotations

import hashlib
import json
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_replicate_yeast_gem import _name, build_draft, load_yeast_scores  # noqa: E402
import cobra  # noqa: E402
from raven_toolbox.localization import certify as C  # noqa: E402

xml = Path("C:/Work/GitHub/yeast-GEM/model/yeast-GEM.xml")
yeast = cobra.io.read_sbml_model(str(xml))
draft, biomass_id = build_draft(yeast)
scores = load_yeast_scores(Path("data/deeploc"))
relocate = [r.id for r in draft.reactions if not r.boundary and r.id != biomass_id]

compartments = sorted(set(draft.compartments) | set(scores.compartments))
sc = C._prepare_scope(draft, scores, relocate, transportable=None,
                      base_metabolite=_name, biomass_reaction=biomass_id, min_growth=0.0)

# Dump the model optlang actually hands to Gurobi, so we can diff it across seeds.
seed = __import__("os").environ.get("PYTHONHASHSEED", "?")
_orig = C._pin_deterministic
def _dump(prob, opt):  # noqa: ANN001
    _orig(prob, opt)
    try:
        opt.problem.update()
        opt.problem.write(f".research_tmp/master_seed{seed}.mps")
    except Exception as e:  # noqa: BLE001
        print("dump failed:", e)
C._pin_deterministic = _dump

t0 = time.monotonic()
status, placements, _ = C._solve_placement_master(
    draft, sc.movable, sc.genes_in_scope, sc.gene_rxns, sc.score_df, compartments,
    multi_compartment_penalty=0.5, forced={}, colocation_groups=[], time_limit=600)
wall = time.monotonic() - t0
place = {rid: sorted(cs) for rid, cs in sorted(placements.items())}
blob = json.dumps(place, sort_keys=True)
n_placed = sum(1 for cs in place.values() if cs)
print(f"status={status}  n_placed={n_placed}  "
      f"sha256={hashlib.sha256(blob.encode()).hexdigest()[:16]}  wall={wall:.1f}s")
