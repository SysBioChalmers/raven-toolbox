"""Export the exact assign_compartments inputs so a MATLAB arm uses identical
scores, relocate set and curated ground truth. Writes to .research_tmp/xarm/.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import cobra

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_replicate_yeast_gem import (  # noqa: E402
    build_draft, load_yeast_scores, curated_reaction_compartments,
    curated_gene_compartments,
)

YEAST = Path("C:/Work/GitHub/yeast-GEM/model/yeast-GEM.xml")
DATA = Path("data/deeploc")
OUT = Path(".research_tmp/xarm")
OUT.mkdir(parents=True, exist_ok=True)

yeast = cobra.io.read_sbml_model(str(YEAST))
curated_growth = float(yeast.slim_optimize())
curated_rxn = curated_reaction_compartments(yeast)      # rxn id -> compartment
curated_gene = curated_gene_compartments(yeast)         # gene -> set of compartments
draft, biomass_id = build_draft(yeast)
scores = load_yeast_scores(DATA)
relocate = [r.id for r in draft.reactions if not r.boundary and r.id != biomass_id]
min_growth = 0.5 * curated_growth

# scores matrix: gene, then one column per compartment (same ids yeast-GEM uses)
df = scores.df
comps = list(df.columns)
with (OUT / "scores.csv").open("w") as fh:
    fh.write("gene," + ",".join(comps) + "\n")
    for gene, row in df.iterrows():
        fh.write(str(gene) + "," + ",".join(f"{row[c]:.6g}" for c in comps) + "\n")

with (OUT / "curated_rxn.csv").open("w") as fh:
    fh.write("rxn,comp\n")
    for r, c in curated_rxn.items():
        fh.write(f"{r},{c}\n")

with (OUT / "curated_gene.csv").open("w") as fh:
    fh.write("gene,comps\n")
    for g, cs in curated_gene.items():
        fh.write(f"{g},{';'.join(sorted(cs))}\n")

(OUT / "meta.json").write_text(json.dumps({
    "biomass_id": biomass_id,
    "curated_growth": curated_growth,
    "min_growth": min_growth,
    "default_compartment": "c",
    "compartments": comps,
    "relocate": relocate,
}, indent=2))

print(f"genes {len(df)}; compartments {comps}; relocate {len(relocate)}; "
      f"curated_growth {curated_growth:.4f}; min_growth {min_growth:.4f}")
print(f"written -> {OUT}")
