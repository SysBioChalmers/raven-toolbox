#!/usr/bin/env python3
"""Export the exact scope behind master_seed1.mps so MATLAB can rebuild the live
placement master on the identical problem: the flattened draft (SBML), the score
matrix, and meta (relocate list, biomass id, compartments)."""
from __future__ import annotations

import json
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_replicate_yeast_gem import _name, build_draft, load_yeast_scores  # noqa: E402
import cobra  # noqa: E402
from raven_toolbox.localization import certify as C  # noqa: E402

OUT = Path(".research_tmp/live")
OUT.mkdir(parents=True, exist_ok=True)

yeast = cobra.io.read_sbml_model("C:/Work/GitHub/yeast-GEM/model/yeast-GEM.xml")
draft, biomass_id = build_draft(yeast)
scores = load_yeast_scores(Path("data/deeploc"))
relocate = [r.id for r in draft.reactions if not r.boundary and r.id != biomass_id]

cobra.io.write_sbml_model(draft, str(OUT / "draft.xml"))

comps = sorted(set(draft.compartments) | set(scores.compartments))
df = scores.df
# scores.csv: gene, <comp1>, <comp2>, ...
lines = ["gene," + ",".join(comps)]
for g in df.index:
    row = [str(g)]
    for c in comps:
        v = df.at[g, c] if c in df.columns else 0.0
        row.append("0" if v is None else repr(float(v)))  # shortest exact float64 round-trip
    lines.append(",".join(row))
(OUT / "scores.csv").write_text("\n".join(lines) + "\n")

json.dump({"relocate": relocate, "biomass_id": biomass_id,
           "compartments": comps, "penalty": 0.5},
          open(OUT / "meta.json", "w"), indent=0)
print(f"exported scope to {OUT}: {len(relocate)} relocate, biomass={biomass_id}, "
      f"comps={comps}, genes_scored={len(df.index)}")
