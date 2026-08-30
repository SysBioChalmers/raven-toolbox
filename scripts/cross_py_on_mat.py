import sys, json, csv
from pathlib import Path
import cobra
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_replicate_yeast_gem import load_yeast_scores, _name
from raven_toolbox.localization import assign_compartments

def main():
    X = Path(".research_tmp/xarm")
    meta = json.loads((X/"meta.json").read_text())
    rel = json.loads((X/"mat_relocate.json").read_text())
    relocate, biomass = rel["relocate"], rel["biomass"]
    draft = cobra.io.read_sbml_model(str(X/"mat_draft.xml"))
    print(f"MATLAB draft in cobra: {len(draft.reactions)} rxns, {len(draft.genes)} genes")
    scores = load_yeast_scores(Path("data/deeploc"))
    scoped = set(g.id for r in draft.reactions if r.id in set(relocate)
                 for g in r.genes if g.id in scores.df.index)
    print(f"Python-on-MATLAB-draft genes_in_scope: {len(scoped)}")
    prop = assign_compartments(draft, scores, relocate, default_compartment="c",
        base_metabolite=_name, biomass_reaction=biomass, min_growth=meta["min_growth"])
    cur = {}
    with (X/"curated_rxn.csv").open() as fh:
        r = csv.reader(fh); next(r)
        for row in r: cur[row[0]] = row[1]
    n = m = 0
    for rid, cs in prop.placements.items():
        if cs and rid in cur:
            n += 1; m += (cs[0] == cur[rid])
    print(f"Python-on-MATLAB-draft certified={prop.certified} rxn agreement vs curated: {100*m/n:.1f}% ({m}/{n})")

if __name__ == "__main__":
    main()
