#!/usr/bin/env python3
"""Score how well each model reproduces canonical yeast localization biology.

Reads the verified marker set (pathways + dual-localized enzymes, from the biology research workflow;
default .research_tmp/biology_markers.json) and the per-model gene->compartment maps
(.research_tmp/compartment_structure.json), and reports, for yeast-GEM / certified / CarveFungi:

  * pathway localization: fraction of each pathway's marker genes placed in the literature compartment,
    scored only over marker genes the model actually contains (coverage reported separately);
  * dual localization: whether each documented dual-localized gene is placed in >=2 of its
    literature compartments.

CarveFungi compartment codes (r=ER, x=peroxisome) are mapped to yeast-GEM's (er, p); yeast-GEM membrane
sub-compartments (erm/mm/gm/vm) are collapsed to their parent organelle. ASCII-only.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

# CarveFungi universal-DB codes -> yeast-GEM codes (for scoring against the same expected compartment)
CF_TO_YG = {"r": "er", "x": "p"}
# yeast-GEM membrane sub-compartments -> parent organelle
YG_PARENT = {"erm": "er", "mm": "m", "gm": "g", "vm": "v"}


def _norm_model_comps(comps, model):
    out = set()
    for c in comps:
        if model == "carvefungi":
            c = CF_TO_YG.get(c, c)
        out.add(YG_PARENT.get(c, c))
    return out


def _expected_for_model(code, model):
    # a curated compartment that CarveFungi cannot represent (g/v/ce/lp) is scored as N/A there
    if model == "carvefungi" and code in {"g", "v", "ce", "lp"}:
        return None
    return code


def score_pathways(markers, struct):
    rows = []
    for p in markers["pathways"]:
        exp = p["expected_compartment"]
        orfs = [g["orf"] for g in p["marker_genes"]]
        rec = {"pathway": p["name"], "expected": exp, "n_markers": len(orfs),
               "confidence": p.get("confidence", "")}
        for model in ("yeast_gem", "certified", "carvefungi"):
            gc = struct[model]["gene_compartments"]
            exp_m = _expected_for_model(exp, model)
            present = [o for o in orfs if o in gc]
            if exp_m is None:
                rec[model] = {"na": True, "coverage": len(present)}
                continue
            hits = sum(1 for o in present if exp_m in _norm_model_comps(gc[o], model))
            rec[model] = {"hit_rate": round(hits / len(present), 3) if present else None,
                          "hits": hits, "present": len(present)}
        rows.append(rec)
    return rows


def score_dual(markers, struct):
    rows = []
    for d in markers["dual_localized"]:
        orf = d["orf"]
        want = set(d["compartments"])
        rec = {"orf": orf, "name": d.get("standard_name", ""), "want": sorted(want),
               "mechanism": d.get("mechanism", "")}
        for model in ("yeast_gem", "certified", "carvefungi"):
            gc = struct[model]["gene_compartments"]
            if orf not in gc:
                rec[model] = "absent"
                continue
            got = _norm_model_comps(gc[orf], model)
            want_m = {_expected_for_model(c, model) for c in want}
            want_m.discard(None)
            overlap = got & want_m
            rec[model] = {"placed": sorted(got), "overlap": sorted(overlap),
                          "captured_dual": len(overlap) >= 2}
        rows.append(rec)
    return rows


def _agg_pathways(rows, model):
    vals = [r[model]["hit_rate"] for r in rows
            if isinstance(r[model], dict) and r[model].get("hit_rate") is not None]
    cov = [r[model].get("present", 0) for r in rows if isinstance(r[model], dict)]
    tot = sum(r["n_markers"] for r in rows)
    return (sum(vals) / len(vals) if vals else float("nan"), sum(cov), tot)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--markers", type=Path, default=Path(".research_tmp/biology_markers.json"))
    ap.add_argument("--structure", type=Path, default=Path(".research_tmp/compartment_structure.json"))
    ap.add_argument("--out", type=Path, default=Path(".research_tmp/biology_score.json"))
    args = ap.parse_args()

    markers = json.loads(args.markers.read_text())
    struct = json.loads(args.structure.read_text())
    pth = score_pathways(markers, struct)
    dual = score_dual(markers, struct)

    print("=== pathway localization (fraction of markers in the literature compartment) ===")
    print(f"{'pathway':32}{'exp':>4}{'yeast-GEM':>11}{'certified':>11}{'CarveFungi':>12}")
    for r in pth:
        def cell(m, r=r):
            v = r[m]
            if not isinstance(v, dict) or v.get("na"):
                return "n/a"
            hr = v.get("hit_rate")
            return f"{hr:.0%}({v['present']})" if hr is not None else "-(0)"
        print(f"{r['pathway'][:31]:32}{r['expected']:>4}{cell('yeast_gem'):>11}"
              f"{cell('certified'):>11}{cell('carvefungi'):>12}")
    for m, lab in [("yeast_gem", "yeast-GEM"), ("certified", "certified"), ("carvefungi", "CarveFungi")]:
        acc, cov, tot = _agg_pathways(pth, m)
        print(f"  {lab}: mean pathway hit-rate {acc:.1%}; marker coverage {cov}/{tot}")

    print("\n=== dual localization (captured >=2 of the literature compartments?) ===")
    print(f"{'gene':10}{'want':>10}  yeast-GEM / certified / CarveFungi")
    for r in dual:
        def dcell(m, r=r):
            v = r[m]
            if v == "absent":
                return "absent"
            return ("DUAL:" + "+".join(v["overlap"])) if v["captured_dual"] else \
                   ("mono:" + "+".join(v["placed"]) if v["placed"] else "none")
        print(f"{(r['name'] or r['orf'])[:9]:10}{'+'.join(r['want']):>10}  "
              f"{dcell('yeast_gem')} / {dcell('certified')} / {dcell('carvefungi')}")
    for m, lab in [("yeast_gem", "yeast-GEM"), ("certified", "certified"), ("carvefungi", "CarveFungi")]:
        present = [r for r in dual if r[m] != "absent"]
        capt = sum(1 for r in present if r[m]["captured_dual"])
        print(f"  {lab}: dual captured {capt}/{len(present)} present ({len(dual)} total)")

    args.out.write_text(json.dumps({"pathways": pth, "dual": dual}, indent=2, default=str))
    print(f"\nwritten -> {args.out}")


if __name__ == "__main__":
    main()
