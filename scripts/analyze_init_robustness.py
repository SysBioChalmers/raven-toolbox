#!/usr/bin/env python3
"""Robustness of (f)tINIT to degraded expression input (Phase 4d.7).

The clean-data calibration (``analyze_init_params.py``) asks "which parameter value is
best on good data". This asks the complementary question: **as the input expression data
gets noisier or sparser, does the pipeline still produce a *functional* model, and which
parameters keep it functional?** "Functional" is measured rigorously — the fraction of the
essential metabolic tasks (``metabolicTasks_Essential.txt``) the extracted model can still
perform, via :func:`ravengem.tasks.check_tasks`.

Three independent degradations of the gene-expression vector (severity = higher is worse):

* ``dropout``    — set a random fraction of genes to 0 (→ gene score -5, a strong *remove*
                   signal). Simulates shallow sequencing / single-cell dropout.
* ``noise``      — multiply each level by ``exp(N(0, sigma))`` (sigma = severity).
                   Simulates measurement noise.
* ``downsample`` — drop a random fraction of genes from the dataset entirely (→ their
                   reactions fall back to ``no_gene_score``). Simulates a smaller panel.

Two phases:

* **gradient** — for each (degradation, severity) run ftINIT *with* and *without* the task
  layer (task-prep vs no-task prep) and record the functional-task pass-rate, reaction
  count, and reaction-set Jaccard vs the clean-data model. Shows where expression-only
  ftINIT loses functionality and whether tasks+gap-fill protect it.
* **rescue** — at a fixed severe degradation, sweep the robustness levers on the *no-task*
  pipeline (``no_gene_score``, ``force_on``) and contrast with the task pipeline, to show
  what restores functionality.

``--algo ftinit`` (default) or ``tinit`` (then ``prod_weight``/``eps`` are the extra
levers). Resumable: each config is pickled and a re-run skips finished ones. Reuses the
cached preps from the Human-GEM validation run (``rg_prep.pkl`` no-task, ``rg_prep_tasks.pkl``
task). Robustness runs use a loose ``mip_gap``/``time_limit`` (functionality, not the exact
optimum, is what matters) so the grid is affordable.

Usage
-----
    python scripts/analyze_init_robustness.py --algo ftinit --phase gradient,rescue \
        --work ~/hgem_compare --cell HCT116
"""
from __future__ import annotations

import argparse
import pickle
import time
from dataclasses import dataclass, field
from pathlib import Path

import cobra
import numpy as np

from ravengem.init import (
    ftinit,
    gene_scores_from_expression,
    get_init_model,
    score_reactions_from_genes,
)
from ravengem.tasks import check_tasks, parse_task_list

# Degradation grid (severity per kind). 0.0 is the shared clean baseline.
GRADIENT = {
    "dropout": (0.25, 0.5, 0.75, 0.9),
    "noise": (0.5, 1.0, 2.0),
    "downsample": (0.3, 0.6, 0.8),
}
RESCUE_KIND, RESCUE_LEVEL = "dropout", 0.75   # a severe-but-not-degenerate point
NO_GENE_SCORES = (-2.0, -1.0, -0.5, -4.0)     # default first
FORCE_ONS = (0.1, 0.02, 0.05, 0.2)
PROD_WEIGHTS = (0.5, 0.0, 1.0, 2.0)           # tINIT only
EPS_VALS = (1.0, 0.5, 0.1)                    # tINIT only

# Loose solver tolerances for the robustness grid (speed; functionality is the metric).
MIP_GAP, TIME_LIMIT = 0.01, 300.0


@dataclass
class Result:
    label: str
    seconds: float
    status: str
    n_rxns: int
    n_pass: int
    n_tasks: int
    frac_pass: float
    reactions: list[str] = field(default_factory=list)
    jaccard_clean: float | None = None


def _jaccard(a: set[str], b: set[str]) -> float:
    return len(a & b) / len(a | b) if (a or b) else 1.0


def degrade(expr: dict[str, float], kind: str, level: float, seed: int) -> dict[str, float]:
    """Return a degraded copy of the expression dict (severity ``level``)."""
    if level <= 0:
        return dict(expr)
    rng = np.random.default_rng(seed)
    genes = list(expr)
    if kind == "dropout":
        out = dict(expr)
        for g in rng.choice(genes, size=int(level * len(genes)), replace=False):
            out[g] = 0.0
        return out
    if kind == "noise":
        return {g: max(v * float(np.exp(rng.normal(0.0, level))), 0.0) for g, v in expr.items()}
    if kind == "downsample":
        keep = set(rng.choice(genes, size=int((1 - level) * len(genes)), replace=False))
        return {g: v for g, v in expr.items() if g in keep}
    raise ValueError(f"unknown degradation kind {kind!r}")


def functionality(model: cobra.Model, tasks) -> tuple[int, int]:
    """(passed, total) essential tasks the extracted model can perform."""
    results = check_tasks(model, tasks)
    return sum(t.passed for t in results), len(results)


def _build_ftinit(prep, ref, expr, *, no_gene_score, force_on):
    g = gene_scores_from_expression(expr, 1.0)
    r = score_reactions_from_genes(ref, g, no_gene_score=no_gene_score)
    return ftinit(prep, r, gene_scores=g, series="1+1", force_on=force_on,
                  mip_gap=MIP_GAP, time_limit=TIME_LIMIT)


def _build_tinit(ref, expr, essential, *, no_gene_score, prod_weight, eps):
    g = gene_scores_from_expression(expr, 1.0)
    r = score_reactions_from_genes(ref, g, no_gene_score=no_gene_score)
    return get_init_model(ref, rxn_scores=r, essential_rxns=essential, prod_weight=prod_weight,
                          eps=eps, mip_gap=MIP_GAP, time_limit=TIME_LIMIT).model


def _measure(label, builder, tasks, clean_set=None) -> Result:
    t = time.time()
    try:
        model = builder()
        n_pass, n_tasks = functionality(model, tasks)
        rset = sorted(x.id for x in model.reactions)
        r = Result(label, time.time() - t, "ok", len(rset), n_pass, n_tasks,
                   n_pass / n_tasks if n_tasks else 0.0, rset)
        if clean_set is not None:
            r.jaccard_clean = _jaccard(set(rset), clean_set)
    except Exception as ex:  # noqa: BLE001  (infeasible/failed build is the headline finding)
        r = Result(label, time.time() - t, f"FAIL:{type(ex).__name__}", 0, 0, len(tasks), 0.0)
    return r


def _table(title, results, note="") -> list[str]:
    lines = [f"### {title}", ""]
    if note:
        lines += [note, ""]
    lines.append("| config | time (s) | status | n_rxns | tasks passed | frac | Jaccard vs clean |")
    lines.append("|---|--:|---|--:|--:|--:|--:|")
    for r in results:
        jac = f"{r.jaccard_clean:.3f}" if r.jaccard_clean is not None else "-"
        lines.append(f"| {r.label} | {r.seconds:.0f} | {r.status} | {r.n_rxns} | "
                     f"{r.n_pass}/{r.n_tasks} | {r.frac_pass:.3f} | {jac} |")
    lines.append("")
    for ln in lines:
        print(ln)
    return lines


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", type=Path, default=Path.home() / "hgem_compare")
    ap.add_argument("--human-gem", type=Path, default=Path.home() / "github" / "Human-GEM")
    ap.add_argument("--cell", default="HCT116")
    ap.add_argument("--algo", choices=("ftinit", "tinit"), default="ftinit")
    ap.add_argument("--phase", default="gradient,rescue", help="gradient,rescue")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--doc", type=Path, default=None)
    args = ap.parse_args()

    out = args.out or args.work / f"init_robustness_{args.algo}_{args.cell}.pkl"
    store: dict = pickle.load(open(out, "rb")) if out.exists() else {}

    def save():
        tmp = Path(f"{out}.part")
        pickle.dump(store, open(tmp, "wb"))
        tmp.replace(out)

    def cached(key, fn):
        if key not in store:
            print(f"[{args.algo}] {key[1]} ...", flush=True)
            store[key] = fn()
            save()
        return store[key]

    t0 = time.time()
    ref = cobra.io.read_sbml_model(str(args.work / "raven_refModel.xml"))
    ref.solver = cobra.Configuration().solver
    expr: dict[str, float] = {}
    with open(args.human_gem / "data" / "datasets" / "Hart2015_RNAseq.txt") as f:
        h = f.readline().rstrip("\n").split("\t")
        c = h.index(args.cell)
        for line in f:
            p = line.rstrip("\n").split("\t")
            expr[p[0]] = float(p[c])
    tasks = parse_task_list(str(args.human_gem / "data" / "metabolicTasks" /
                                "metabolicTasks_Essential.txt"))
    prep_nt = pickle.load(open(args.work / "rg_prep.pkl", "rb"))
    prep_tk = pickle.load(open(args.work / "rg_prep_tasks.pkl", "rb"))
    essential = list(prep_nt.essential_rxns)
    print(f"[{time.time()-t0:.0f}s] ref {len(ref.reactions)} rxns, {len(tasks)} tasks, "
          f"cell={args.cell}, algo={args.algo}", flush=True)

    def build(prep, e, **kw):
        if args.algo == "ftinit":
            return lambda: _build_ftinit(prep, ref, e, no_gene_score=kw.get("no_gene_score", -2.0),
                                         force_on=kw.get("force_on", 0.1))
        # tINIT ignores prep; its analogue of the task layer is forcing task-essential
        # reactions kept. The task pipeline (prep_tk) → pass essential; no-task → none.
        ess = essential if prep is prep_tk else []
        return lambda: _build_tinit(ref, e, ess, no_gene_score=kw.get("no_gene_score", -2.0),
                                    prod_weight=kw.get("prod_weight", 0.5), eps=kw.get("eps", 1.0))

    phases = set(args.phase.split(","))
    doc = [f"# (f)tINIT robustness to degraded input — Human-GEM / {args.cell} / {args.algo}", "",
           "Functional = fraction of essential metabolic tasks the extracted model can perform "
           "(check_tasks). Configs use a loose MIP gap (speed). Generated by "
           "`scripts/analyze_init_robustness.py`.", ""]

    # Clean baselines (level 0) per pipeline, for the Jaccard reference.
    clean_nt = cached(("clean", "no-task clean"), lambda: _measure(
        "clean (no-task)", build(prep_nt, expr), tasks))
    clean_tk = cached(("clean", "task clean"), lambda: _measure(
        "clean (task)", build(prep_tk, expr), tasks))
    clean_set_nt, clean_set_tk = set(clean_nt.reactions), set(clean_tk.reactions)
    doc += _table("Clean-data baseline (no degradation)", [clean_nt, clean_tk])

    if "gradient" in phases:
        for kind, levels in GRADIENT.items():
            rows = [clean_nt, clean_tk]
            for lvl in levels:
                e = degrade(expr, kind, lvl, args.seed)
                rows.append(cached((f"grad_{kind}", f"no-task {kind}={lvl}"), lambda e=e, lvl=lvl, kind=kind:
                            _measure(f"no-task {kind}={lvl}", build(prep_nt, e), tasks, clean_set_nt)))
                rows.append(cached((f"grad_{kind}", f"task {kind}={lvl}"), lambda e=e, lvl=lvl, kind=kind:
                            _measure(f"task {kind}={lvl}", build(prep_tk, e), tasks, clean_set_tk)))
            doc += _table(f"Gradient: {kind} (no-task vs task pipeline)", rows,
                          "Higher severity = noisier/sparser input. Watch frac (functional "
                          "task pass-rate): the gap between no-task and task rows is what the "
                          "task+gap-fill layer buys.")

    if "rescue" in phases:
        e = degrade(expr, RESCUE_KIND, RESCUE_LEVEL, args.seed)
        rows = []
        tag = f"{RESCUE_KIND}={RESCUE_LEVEL}"
        if args.algo == "ftinit":
            for ngs in NO_GENE_SCORES:
                rows.append(cached(("rescue", f"no-task no_gene_score={ngs}"), lambda ngs=ngs:
                            _measure(f"no-task no_gene_score={ngs}", build(prep_nt, e, no_gene_score=ngs),
                                     tasks, clean_set_nt)))
            for fo in FORCE_ONS:
                rows.append(cached(("rescue", f"no-task force_on={fo}"), lambda fo=fo:
                            _measure(f"no-task force_on={fo}", build(prep_nt, e, force_on=fo),
                                     tasks, clean_set_nt)))
        else:
            for pw in PROD_WEIGHTS:
                rows.append(cached(("rescue", f"tinit prod_weight={pw}"), lambda pw=pw:
                            _measure(f"prod_weight={pw}", build(ref, e, prod_weight=pw), tasks)))
            for ev in EPS_VALS:
                rows.append(cached(("rescue", f"tinit eps={ev}"), lambda ev=ev:
                            _measure(f"eps={ev}", build(ref, e, eps=ev), tasks)))
        rows.append(cached(("rescue", "task pipeline"), lambda:
                    _measure("task pipeline (gap-fill)", build(prep_tk, e), tasks, clean_set_tk)))
        doc += _table(f"Rescue at {tag}: which lever restores functionality?", rows,
                      "Baseline no-task at this severity is the first row of the matching gradient "
                      "table; the task pipeline (last row) is the reference 'robust' result.")

    if args.doc:
        args.doc.write_text("\n".join(doc) + "\n")
        print(f"\nwrote {args.doc}", flush=True)


if __name__ == "__main__":
    main()
