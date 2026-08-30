import sys, json
from pathlib import Path
import cobra
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_replicate_yeast_gem import build_draft
yeast = cobra.io.read_sbml_model("C:/Work/GitHub/yeast-GEM/model/yeast-GEM.xml")
draft, biomass_id = build_draft(yeast)
out = {
 "n_rxns": len(draft.reactions),
 "n_mets": len(draft.metabolites),
 "n_genes": len(draft.genes),
 "rxn_ids": sorted(r.id for r in draft.reactions),
 "met_names": sorted(set(m.name for m in draft.metabolites)),
 "biomass": biomass_id,
}
Path(".research_tmp/xarm/py_draft.json").write_text(json.dumps(out))
print(f"python draft: {out['n_rxns']} rxns, {out['n_mets']} mets, {out['n_genes']} genes")
