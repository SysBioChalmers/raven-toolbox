#!/usr/bin/env python3
"""Run the FULL Python assign_compartments pipeline on the exported draft and dump
the final result (placement per movable reaction, added transports, growth) so it
can be compared byte-for-byte against the MATLAB pipeline on the same draft."""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_replicate_yeast_gem import _name, load_yeast_scores  # noqa: E402
import cobra  # noqa: E402
from raven_toolbox.localization import assign_compartments  # noqa: E402

def main():
  # (guarded so cobra's find_blocked_reactions can spawn workers on Windows)
  X = Path(".research_tmp/live")
  meta = json.loads((X / "meta.json").read_text())
  draft = cobra.io.read_sbml_model(str(X / "draft.xml"))
  scores = load_yeast_scores(Path("data/deeploc"))
  biomass_id = meta["biomass_id"]
  relocate = meta["relocate"]
  min_growth = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0
  max_rounds = int(sys.argv[2]) if len(sys.argv) > 2 else 8

  proposal = assign_compartments(
    draft, scores, relocate, transport_cost=0.5, default_compartment="c",
    base_metabolite=_name, biomass_reaction=biomass_id, min_growth=min_growth,
    time_limit=600, multi_localize=False, max_rounds=max_rounds,
    prune_transports=False, minimize_transports=False)

  place = {rid: sorted(cs) for rid, cs in sorted(proposal.placements.items()) if cs}
  out = {
    "status": proposal.status,
    "objective": round(float(proposal.objective), 6),
    "n_placed": len(place),
    "placement": place,
    "added_transports": sorted(getattr(proposal, "added_transports", []) or []),
    "added_reactions": sorted(getattr(proposal, "added_reactions", []) or []),
    "unplaced": sorted(proposal.unplaced_reactions or []),
    "min_growth": min_growth,
  }
  tag = f"_mg{min_growth:g}".replace(".", "p")
  Path(f".research_tmp/py_full{tag}.json").write_text(json.dumps(out, indent=1, sort_keys=True))
  # flat placement for easy diff: "rid<TAB>comp" per single-comp reaction
  flat = [f"{rid}\t{cs[0]}" for rid, cs in place.items() if len(cs) == 1]
  Path(f".research_tmp/py_full_placement{tag}.tsv").write_text("\n".join(sorted(flat)) + "\n")
  tflat = sorted(f"{t[0]}\t{t[1]}" for t in out["added_transports"])
  Path(f".research_tmp/py_transports{tag}.tsv").write_text("\n".join(tflat) + "\n")
  print(f"status={out['status']} obj={out['objective']} n_placed={out['n_placed']} "
        f"transports={len(out['added_transports'])} added_rxns={len(out['added_reactions'])} "
        f"unplaced={len(out['unplaced'])} tag={tag}")
  print("growths:", getattr(proposal, "growths", None), " min_growth_used:", proposal.min_growth)


if __name__ == "__main__":
    main()
