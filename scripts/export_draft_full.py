import json
from pathlib import Path
import cobra, warnings
warnings.filterwarnings('ignore')
d = cobra.io.read_sbml_model('.research_tmp/xarm/py_draft.xml')
out={}
for r in d.reactions:
    stoich={m.name: round(c,9) for m,c in r.metabolites.items()}
    out[r.id]={'s':stoich,'lb':r.lower_bound,'ub':r.upper_bound,'gpr':r.gene_reaction_rule}
Path('.research_tmp/xarm/py_draft_full.json').write_text(json.dumps(out))
print('exported full py draft', len(out))
