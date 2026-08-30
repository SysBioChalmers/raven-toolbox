import json
from pathlib import Path
import cobra
X = Path(".research_tmp/xarm")
py = json.loads((X/"py_draft.json").read_text())
mat = cobra.io.read_sbml_model(str(X/"mat_draft.xml"))
pyR = set(py["rxn_ids"]); matR = set(r.id for r in mat.reactions)
pyM = set(py["met_names"]); matM = set(m.name for m in mat.metabolites)
print(f"reactions: py {len(pyR)}  mat {len(matR)}  py-only {len(pyR-matR)}  mat-only {len(matR-pyR)}")
print(f"  py-only: {sorted(pyR-matR)[:10]}")
print(f"  mat-only: {sorted(matR-pyR)[:10]}")
print(f"met NAMES: py {len(pyM)}  mat {len(matM)}  py-only {len(pyM-matM)}  mat-only {len(matM-pyM)}")
print(f"  py-only mets: {sorted(pyM-matM)[:8]}")
print(f"  mat-only mets: {sorted(matM-pyM)[:8]}")
# genes
pyG = set(py.get("gene_ids", [])) if "gene_ids" in py else None
matG = set(g.id for g in mat.genes)
print(f"genes: mat {len(matG)}")
