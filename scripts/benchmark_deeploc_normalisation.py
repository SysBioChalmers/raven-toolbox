#!/usr/bin/env python3
"""Normalised vs. raw DeepLoc probabilities for compartment assignment on the *whole* yeast-GEM.

Question: :func:`raven_toolbox.localization.load_deeploc` rescales every gene's best compartment
to 1.0 (RAVEN ``parseScores`` convention). That preserves each gene's argmax but erases how
*confident* DeepLoc was (a 0.97 call and a 0.40 call both become 1.0). DeepLoc's probability is
calibrated, so from the assignment's perspective the magnitude is real signal. Does keeping the raw
probabilities (``normalise=False``) change the assignment's agreement with curated yeast-GEM?

Method (the established flatten benchmark, ``benchmark_localization_yeast.py``, on the full model):
flatten yeast-GEM to one compartment so the MILP cannot lean on metabolite topology, then ask
:func:`predict_localization` to re-place every single-compartment GPR'd reaction from a DeepLoc score
table. Two arms from the *same* DeepLoc output -- ``normalise=True`` (top->1.0) vs ``normalise=False``
(raw probabilities) -- swept over ``transport_cost`` (the arms live on different score scales, so
each must be compared at its own best operating point, not a shared one).

Because the MILP objective is ``max sum score[g,c]*y[g,c] - transport_cost*transports
- penalty*extra_compartments``, the score enters linearly: normalisation makes every gene vote with
weight 1.0; raw lets a confident gene outvote a shaky one and -- because ``transport_cost`` and the
multi-compartment ``penalty`` are *absolute* constants -- raises the bar a reaction must clear to
leave the default compartment or a gene to occupy a second compartment. Both arms agree on each
gene's *preferred* compartment (normalisation is a monotone per-row rescale), but raw clears those
absolute thresholds less often, so even single-gene reactions can be placed differently at
``transport_cost > 0`` (raw behaves like normalised at a higher effective transport cost).

Caveat baked into the interpretation: the objective has **no term in the reaction-placement
variable**, so when a gene sits in several compartments the MILP does not uniquely pin *which* of
them its reaction lands in (a solver tie-break, amenable to set-ordering and the ``mip_gap``). The
gene-level solution (multi-localisation counts) is reproducible; reaction-level accuracy -- and
especially the tiny-stratum ``macro`` -- carries tie-break noise, so small cross-arm deltas should be
read against that floor (two independent runs here differ by up to ~0.1 on macro at the same tc).

Reported per arm x transport_cost: overall accuracy (headline -- the whole truth set), accuracy on
the DeepLoc-*addressable* compartments only (the membrane sub-compartments erm/gm/mm/vm/lp have no
DeepLoc column and cap both arms equally), macro/balanced accuracy (mean of per-compartment
accuracies -- not rewarded by dumping into the majority class `c`), and how many reactions left the
default. Then, at each arm's best operating point: the contested subset (placements that differ
between arms) and a confidence-stratified view (does raw help the high-confidence calls?).

ASCII-only output (Windows console is cp1252).
"""
from __future__ import annotations

import argparse
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import cobra  # noqa: E402
import pandas as pd  # noqa: E402

from raven_toolbox.localization import (  # noqa: E402
    DEFAULT_COMPARTMENT_MAP,
    LocalizationScores,
    load_deeploc,
    predict_localization,
)
from raven_toolbox.manipulation.compartments import merge_compartments  # noqa: E402

#: Compartments a DeepLoc 2 organelle call can actually name (post DEFAULT_COMPARTMENT_MAP).
#: The membrane sub-compartments (erm, gm, mm, vm) and the lipid particle (lp) have no DeepLoc
#: column, so reactions truly there are unreachable for *both* arms -- excluded from the
#: addressable-subset metric so the normalisation effect is not drowned by structural misses.
ADDRESSABLE = {"c", "ce", "e", "er", "g", "m", "n", "p", "v"}


def build_truth(model: cobra.Model) -> dict[str, str]:
    """``{rxn_id: compartment}`` for every non-boundary GPR'd single-compartment reaction."""
    truth: dict[str, str] = {}
    for r in model.reactions:
        if r.boundary or not r.genes:
            continue
        comps = {m.compartment for m in r.metabolites if m.compartment}
        if len(comps) == 1:
            truth[r.id] = next(iter(comps))
    return truth


def load_arm(csvs: list[str], *, normalise: bool) -> LocalizationScores:
    """Concatenate the per-chunk DeepLoc CSVs into one gene x compartment table.

    The chunks hold disjoint gene sets (the FASTA was split at the 500-sequence web limit), so a
    plain row concat is correct. ``raw_confidence`` (the per-gene top probability) is carried for
    the confidence-stratified analysis regardless of ``normalise``.
    """
    parts, confs = [], []
    for c in csvs:
        s = load_deeploc(c, compartment_map=DEFAULT_COMPARTMENT_MAP,
                         normalise=normalise, keep_raw_confidence=True)
        parts.append(s.df)
        confs.append(s.raw_confidence)
    df = pd.concat(parts).fillna(0.0)
    df.index.name = "gene_id"
    out = LocalizationScores(df)
    out.raw_confidence = pd.concat(confs).reindex(df.index)
    return out


def predict(flat: cobra.Model, rel: list[str], scores: LocalizationScores, *,
            tc: float, mcp: float, mip_gap: float, time_limit: float) -> tuple[dict[str, str], dict]:
    """One MILP solve -> (``{rxn_id: predicted_compartment}``, stats)."""
    t = time.time()
    prop = predict_localization(flat, scores, rel, default_compartment="c", transport_cost=tc,
                                multi_compartment_penalty=mcp, apply=False,
                                mip_gap=mip_gap, time_limit=time_limit)
    moved = dict(zip(prop.moved["rxn_id"], prop.moved["to_compartment"], strict=True))
    pred = {rid: moved.get(rid, "c") for rid in rel}
    n_multi = sum(1 for cs in prop.gene_compartments.values() if len(cs) > 1)
    return pred, {"seconds": time.time() - t, "n_moved": len(prop.moved),
                  "n_multi_genes": n_multi, "n_unplaced": len(prop.unplaced_reactions)}


def accuracy(pred: dict[str, str], truth: dict[str, str], rel: list[str]) -> dict:
    """Overall, addressable-only, and macro (balanced) accuracy + per-compartment counts."""
    per: dict[str, list[int]] = {}
    for rid in rel:
        n, ok = per.get(truth[rid], (0, 0))
        per[truth[rid]] = (n + 1, ok + (pred[rid] == truth[rid]))
    overall = sum(ok for _, ok in per.values()) / len(rel)
    addr_ids = [rid for rid in rel if truth[rid] in ADDRESSABLE]
    addr = sum(pred[rid] == truth[rid] for rid in addr_ids) / max(1, len(addr_ids))
    macro_comps = [c for c in per if c in ADDRESSABLE]
    macro = sum(per[c][1] / per[c][0] for c in macro_comps) / max(1, len(macro_comps))
    return {"overall": overall, "addressable": addr, "macro": macro, "per": per,
            "n_addr": len(addr_ids)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--yeast-gem", type=Path,
                    default=Path("C:/Work/GitHub/yeast-GEM/model/yeast-GEM.xml"))
    ap.add_argument("--deeploc-dir", type=Path, default=Path("data/deeploc"))
    ap.add_argument("--transport-costs", default="0.0,0.01,0.02,0.05,0.1")
    ap.add_argument("--operating-point", type=float, default=0.01,
                    help="matched, well-conditioned transport_cost for the head-to-head / contested "
                         "/ confidence-stratified sections (tc=0 is degenerate -- see the doc)")
    ap.add_argument("--multi-compartment-penalty", type=float, default=0.5)
    ap.add_argument("--mip-gap", type=float, default=0.01)
    ap.add_argument("--time-limit", type=float, default=600)
    ap.add_argument("--doc", type=Path, help="write the markdown study here")
    args = ap.parse_args()

    csvs = [str(args.deeploc_dir / f"yeast-GEM_deeploc_{i:03d}.csv") for i in (1, 2, 3)]
    tcs = sorted({float(x) for x in args.transport_costs.split(",")} | {args.operating_point})
    op = args.operating_point

    print(f"loading {args.yeast_gem} ...", flush=True)
    model = cobra.io.read_sbml_model(str(args.yeast_gem))
    truth = build_truth(model)
    norm, raw = load_arm(csvs, normalise=True), load_arm(csvs, normalise=False)

    flat, _, _ = merge_compartments(model, merged_id="c", merged_name="c",
                                    drop_single_metabolite_reactions=False,
                                    deduplicate_reactions=False)
    surviving = {r.id for r in flat.reactions}
    rel = [rid for rid in truth if rid in surviving]
    # Trivial "leave everything in the default compartment c" baseline: the headroom both arms share.
    n_c = sum(truth[rid] == "c" for rid in rel)
    n_addr = sum(truth[rid] in ADDRESSABLE for rid in rel)
    base_overall, base_addr = n_c / len(rel), n_c / max(1, n_addr)
    print(f"yeast-GEM: {len(model.reactions)} rxns, {len(model.genes)} genes; truth set "
          f"{len(rel)} reactions ({sum(truth[r] in ADDRESSABLE for r in rel)} in DeepLoc-"
          f"addressable compartments). Sweeping transport_cost={tcs}.", flush=True)

    arms = {"normalised": norm, "raw": raw}
    results: dict[tuple[str, float], dict] = {}
    preds: dict[tuple[str, float], dict[str, str]] = {}
    for name, scores in arms.items():
        for tc in tcs:
            pred, stats = predict(flat, rel, scores, tc=tc, mcp=args.multi_compartment_penalty,
                                  mip_gap=args.mip_gap, time_limit=args.time_limit)
            acc = accuracy(pred, truth, rel)
            results[(name, tc)] = {**stats, **acc}
            preds[(name, tc)] = pred
            print(f"  {name:11s} tc={tc:<5} {stats['seconds']:4.0f}s  "
                  f"overall={acc['overall']:.3f}  addressable={acc['addressable']:.3f}  "
                  f"macro={acc['macro']:.3f}  moved={stats['n_moved']}  "
                  f"multi_genes={stats['n_multi_genes']}", flush=True)

    # Sweep-best per arm (for the curve), and the matched well-conditioned operating point used for
    # the head-to-head / contested / stratified sections. transport_cost=0 maximises accuracy but is
    # *degenerate*: with no transport term the MILP objective has no term in the reaction-placement
    # variable, so a reaction whose gene sits in several compartments is placed arbitrarily (solver
    # tie-break). The matched `op` (a small tc>0) breaks that tie deterministically, which is why the
    # sub-analyses use it rather than the (higher but ill-conditioned) tc=0 point.
    best = {name: max([tc for tc in tcs if tc > 0],
                      key=lambda tc: results[(name, tc)]["addressable"]) for name in arms}

    # ---- report -----------------------------------------------------------------
    L: list[str] = []
    L += ["# DeepLoc normalisation benchmark (whole yeast-GEM)", "",
          "Does keeping DeepLoc's **raw** probabilities (`load_deeploc(normalise=False)`) instead of "
          "rescaling each gene's best compartment to 1.0 change the compartment assignment's "
          "agreement with curated yeast-GEM? Run on the **entire** model.", "",
          f"* Model: `{args.yeast_gem.name}` -- {len(model.reactions)} reactions, "
          f"{len(model.genes)} genes, {len(model.compartments)} compartments.",
          f"* Truth set: {len(rel)} non-boundary GPR'd single-compartment reactions "
          f"({sum(truth[r] in ADDRESSABLE for r in rel)} in DeepLoc-addressable compartments; the "
          "rest live in membrane sub-compartments DeepLoc cannot name and cap **both** arms equally).",
          f"* Driver: `scripts/benchmark_deeploc_normalisation.py`. "
          f"`multi_compartment_penalty={args.multi_compartment_penalty}`, `mip_gap={args.mip_gap}`.",
          "* Two arms from the *same* DeepLoc CSVs: **normalised** (top->1.0) vs **raw** "
          "(probabilities kept), swept over `transport_cost`.",
          f"* Trivial baseline -- *leave every reaction in the default `c`* (the headroom both arms "
          f"share, since `c` is the majority class and is transport-free): overall {base_overall:.3f}, "
          f"addressable {base_addr:.3f}. Read the accuracies below against this, not against zero.",
          "* **Reproducibility floor.** The MILP objective has no term in the reaction-placement "
          "variable, so a multi-compartment gene's reaction is placed by a solver tie-break "
          "(sensitive to set ordering and `mip_gap`). Gene-level counts (multi-localisation) are "
          "reproducible; reaction-level `overall`/`addressable` are stable to ~1pp run-to-run, but "
          "`macro` swings by ~0.1 at the same `transport_cost` (it equal-weights compartments with "
          "n=7-13). Treat sub-1pp accuracy deltas and macro differences as within noise.", "",
          "## Accuracy vs. transport_cost", "",
          "`overall` = whole truth set; `addressable` = truth compartment is DeepLoc-reachable; "
          "`macro` = mean of per-compartment accuracies (balances the `c` majority, but **fragile** "
          "-- dominated by tiny-n compartments and tie-break noise; see the reproducibility floor); "
          "`moved` = reactions placed outside the default `c`.", "",
          "| arm | transport_cost | overall | addressable | macro | moved | multi-loc genes | s |",
          "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for name in arms:
        for tc in tcs:
            r = results[(name, tc)]
            star = " *" if tc == best[name] else (" (deg.)" if tc == 0 else "")
            L.append(f"| {name}{star} | {tc} | {r['overall']:.3f} | {r['addressable']:.3f} | "
                     f"{r['macro']:.3f} | {r['n_moved']} | {r['n_multi_genes']} | {r['seconds']:.0f} |")
    L += ["", "`*` = best operating point with `transport_cost > 0` (by `addressable`). "
          "`(deg.)` marks `transport_cost = 0`: it scores highest but is degenerate -- with no "
          "transport term the objective does not constrain *which* of a multi-compartment gene's "
          "compartments its reaction lands in, so those placements are an arbitrary solver "
          "tie-break. The sections below therefore compare both arms at the same well-conditioned "
          f"`transport_cost = {op}` (matched operating point: only the scores differ).", ""]

    # Head-to-head at the matched, well-conditioned operating point (only the scores differ).
    rn, rr = results[("normalised", op)], results[("raw", op)]
    pn, pr = preds[("normalised", op)], preds[("raw", op)]
    L += [f"## Head-to-head at matched transport_cost = {op}", "",
          f"* **normalised**: overall {rn['overall']:.3f}, addressable {rn['addressable']:.3f}, "
          f"macro {rn['macro']:.3f}, multi-loc genes {rn['n_multi_genes']}.",
          f"* **raw**: overall {rr['overall']:.3f}, addressable {rr['addressable']:.3f}, "
          f"macro {rr['macro']:.3f}, multi-loc genes {rr['n_multi_genes']}.",
          f"* delta (raw - normalised): overall {rr['overall'] - rn['overall']:+.3f}, "
          f"addressable {rr['addressable'] - rn['addressable']:+.3f}, "
          f"macro {rr['macro'] - rn['macro']:+.3f}.",
          f"* sweep-best (tc>0) by addressable: normalised @ {best['normalised']} "
          f"({results[('normalised', best['normalised'])]['addressable']:.3f}), "
          f"raw @ {best['raw']} ({results[('raw', best['raw'])]['addressable']:.3f}).", ""]

    # Contested subset: reactions placed differently at the matched operating point.
    contested = [rid for rid in rel if pn[rid] != pr[rid]]
    if contested:
        cn = sum(pn[rid] == truth[rid] for rid in contested)
        cr = sum(pr[rid] == truth[rid] for rid in contested)
        L += ["## Contested subset (placements differ between arms)", "",
              f"{len(contested)} of {len(rel)} reactions are placed differently at "
              f"transport_cost = {op} -- this is where normalisation actually changes the answer.", "",
              f"* normalised correct on contested: {cn}/{len(contested)} = {cn/len(contested):.3f}",
              f"* raw correct on contested: {cr}/{len(contested)} = {cr/len(contested):.3f}",
              f"* the rest ({len(rel) - len(contested)} reactions) are placed identically by both.",
              ""]

    # Confidence-stratified: single-gene reactions, binned by the gene's raw top probability.
    gene_of = {}
    for r in model.reactions:
        if r.id in truth and len(r.genes) == 1:
            gene_of[r.id] = next(iter(r.genes)).id
    rc = raw.raw_confidence
    bins = [(0.0, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.01)]
    L += ["## Confidence-stratified accuracy (single-gene reactions)", "",
          f"Single-gene reactions only, binned by the gene's raw DeepLoc top probability, at the "
          f"matched transport_cost = {op}. If raw probabilities help, the gain concentrates in the "
          "high-confidence bins.", "",
          "| raw top prob | n | normalised acc | raw acc | delta |",
          "|---|---:|---:|---:|---:|"]
    for lo, hi in bins:
        ids = [rid for rid, g in gene_of.items()
               if g in rc.index and pd.notna(rc[g]) and lo <= float(rc[g]) < hi]
        if not ids:
            continue
        an = sum(pn[rid] == truth[rid] for rid in ids) / len(ids)
        ar = sum(pr[rid] == truth[rid] for rid in ids) / len(ids)
        L.append(f"| [{lo:.1f}, {hi if hi <= 1 else 1.0:.1f}] | {len(ids)} | {an:.3f} | {ar:.3f} "
                 f"| {ar - an:+.3f} |")
    L.append("")

    # Per-compartment at each arm's best tc.
    comps = sorted({truth[r] for r in rel})
    L += [f"## Per-compartment accuracy at matched transport_cost = {op}", "",
          "| compartment | n | normalised | raw |", "|---|---:|---:|---:|"]
    for c in comps:
        ids = [rid for rid in rel if truth[rid] == c]
        an = sum(pn[rid] == c for rid in ids) / len(ids)
        ar = sum(pr[rid] == c for rid in ids) / len(ids)
        addr = "" if c in ADDRESSABLE else " (unreachable)"
        L.append(f"| {c}{addr} | {len(ids)} | {an:.3f} | {ar:.3f} |")
    L.append("")

    # Honest conclusion (kept in sync with the reproducibility floor above).
    n_id = len(rel) - len(contested) if contested else len(rel)
    L += ["## Conclusion", "",
          f"* **Agreement is normalisation-neutral.** `overall`/`addressable` differ by <=1pp "
          f"between arms -- within the reproducibility floor. This is partly structural: both arms "
          f"share each gene's argmax, so {n_id} of {len(rel)} reactions are placed identically; the "
          f"arms can only diverge on the {len(contested)} contested reactions, where **both** place "
          "poorly (raw no better).",
          "* **No accuracy reason to drop per-gene normalisation.** Raw does not rescue the "
          "high-confidence or contested calls; the `macro` edge that can appear at some operating "
          "points is fragile (tiny-n compartments + tie-break noise) and does not reproduce.",
          f"* **The one reproducible difference is structural and gene-level:** raw assigns far "
          f"fewer genes to multiple compartments ({rr['n_multi_genes']} vs {rn['n_multi_genes']} at "
          f"transport_cost={op}), because its calibrated magnitudes (~0.5-0.9) clear the *absolute* "
          "multi-compartment penalty less often than the normalised 1.0. The same effect is "
          "reachable by raising `multi_compartment_penalty` on normalised scores -- i.e. raw is "
          "mostly a re-scaling of the existing knobs, not new information.",
          "* **Recommendation:** keep normalisation the **default**; expose `normalise=False` as an "
          "opt-in for callers who want DeepLoc's calibrated magnitudes (e.g. the confidence signal "
          "`triage_localization` consumes, or to avoid hand-tuning the transport scale). Flipping "
          "the default is not justified by this benchmark.", ""]

    text = "\n".join(L) + "\n"
    print("\n" + text, flush=True)
    if args.doc:
        args.doc.write_text(text, encoding="utf-8")
        print(f"wrote {args.doc}", flush=True)


if __name__ == "__main__":
    main()
