import sys, json
from pathlib import Path
import cobra
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_replicate_yeast_gem import build_draft, load_yeast_scores
yeast = cobra.io.read_sbml_model("C:/Work/GitHub/yeast-GEM/model/yeast-GEM.xml")
draft, biomass_id = build_draft(yeast)
scores = load_yeast_scores(Path("data/deeploc"))
relocate = [r.id for r in draft.reactions if not r.boundary and r.id != biomass_id]
# genes_in_scope: scored genes on movable reactions
movable = [r for r in draft.reactions if r.id in set(relocate)]
scoped = set()
for r in movable:
    for g in r.genes:
        if g.id in scores.df.index: scoped.add(g.id)
Path(".research_tmp/xarm/py_scope.json").write_text(json.dumps({
  "movable": sorted(relocate), "genes_in_scope": sorted(scoped)}))
print(f"python: movable {len(relocate)}, genes_in_scope {len(scoped)}")
