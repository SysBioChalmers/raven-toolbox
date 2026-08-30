#!/usr/bin/env python3
"""Run assign_compartments once and dump a hash of the placement, to check that
the result no longer depends on Python's per-process set-iteration order."""
from __future__ import annotations

import hashlib
import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_replicate_yeast_gem import (  # noqa: E402
    _name,
    build_draft,
    load_yeast_scores,
)
import cobra  # noqa: E402
from raven_toolbox.localization import assign_compartments  # noqa: E402

DATA = Path("data/deeploc")
xml = Path("C:/Work/GitHub/yeast-GEM/model/yeast-GEM.xml")
yeast = cobra.io.read_sbml_model(str(xml))
draft, biomass_id = build_draft(yeast)
scores = load_yeast_scores(DATA)
relocate = [r.id for r in draft.reactions
            if not r.boundary and r.id != biomass_id]
import time  # noqa: E402
t0 = time.monotonic()
proposal = assign_compartments(
    draft, scores, relocate, transport_cost=0.0, default_compartment="c",
    base_metabolite=_name, biomass_reaction=biomass_id, min_growth=0.0,
    time_limit=600, multi_localize=False)
wall = time.monotonic() - t0
place = {rid: sorted(cs) for rid, cs in sorted(proposal.placements.items())}
blob = json.dumps(place, sort_keys=True)
print(f"objective={proposal.objective:.6f}  n_placed={len(place)}  "
      f"sha256={hashlib.sha256(blob.encode()).hexdigest()[:16]}  wall={wall:.1f}s")
