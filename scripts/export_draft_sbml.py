import sys
from pathlib import Path
import cobra
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_replicate_yeast_gem import build_draft
yeast = cobra.io.read_sbml_model("C:/Work/GitHub/yeast-GEM/model/yeast-GEM.xml")
draft, biomass_id = build_draft(yeast)
cobra.io.write_sbml_model(draft, ".research_tmp/xarm/py_draft.xml")
print(f"wrote draft SBML: {len(draft.reactions)} rxns, {len(draft.genes)} genes, biomass {biomass_id}")
