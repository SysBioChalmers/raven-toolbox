#!/usr/bin/env python3
"""Regenerate a genuine, functional, gene-annotated fungal reconstruction carve (yeast) for benchmarking.

Drives the reference reconstruction tool's OWN pipeline (its EggNOG-based scoring, its ``carve_model``
MILP, its gene-annotation step) end to end, so the output is what the tool actually produces --
unlike a bare reaction-id cache (which drops the tool's uptake reactions, loses each reaction's solved
direction, and carries no genes, so it cannot grow standalone).

Steps (mirrors the tool's own ``carveFungi()`` orchestration, minus its ensemble/ensemble-averaging and
its Unix-``sed`` SBML-annotation post-step -- both cosmetic, not needed for a functional model):

1. Score every universal-database reaction from the tool's shipped EggNOG functional annotations
   (``get_gene_to_score_dict`` + ``assign_scores``), then refine by localisation
   (``scoring_compartments``). The tool's shipped localisation file is RefSeq-keyed while its EggNOG
   file is ORF-keyed (zero overlap -- confirmed dead on arrival), so this substitutes an ORF-keyed
   compartment file built from this project's own yeast DeepLoc predictions (disclosed below).
2. ``carve_model`` -- the tool's own MILP (its ``minmax_reduction``) on a defined minimal medium
   (glucose/ammonium/phosphate/sulfate/Fe2+/O2/water -- the tool's own hardcoded constrained-medium
   exchange set), returning up to 5 pool solutions as functional models; keep the best-objective one.
3. ``add_gene_annotation`` -- GPRs from the EC-to-gene map (the tool's own function, pure cobra).
4. ``add_missing_reactions`` -- the tool's own biomass/exchange completion step.
5. Verify the model grows and has genes; write it to ``--out``.

Needs: the tool's clone (``--carvefungi-dir``), a CPLEX-capable Python (this project's ``.research_tmp/
cpxenv``), and this project's yeast DeepLoc predictions (``data/deeploc/``). ASCII-only output.

Usage::

    .research_tmp/cpxenv/Scripts/python.exe scripts/build_reference_carve_model.py \\
        --carvefungi-dir C:/Work/GitHub/CarveFungi --out .research_tmp/cf_cache/reference_carve_model.sbml
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# DeepLoc column -> the tool's 4-category compartment-scoring label (E=ER, M=mito, P=peroxisome,
# O=other/catch-all for cytosol+nucleus+Golgi). Cutoffs below match the tool's own per-category
# thresholds in EggNogScoring.scoring_compartments (get_final_score calls): M/O 0.45, P 0.1, E 0.2.
_DEEPLOC_DIRECT = {"E": "Endoplasmic reticulum", "M": "Mitochondrion", "P": "Peroxisome"}
_DEEPLOC_OTHER = ("Cytoplasm", "Nucleus", "Golgi apparatus")
_NEUTRAL_ROW = {"E": 0.2, "M": 0.45, "P": 0.1, "O": 0.45}  # == the category cutoffs -> compartment_score
# = 0 -> scoring_compartments takes its `else` branch and keeps the EC-only score unmodified: the
# correct "no localisation evidence" behaviour for a gene outside DeepLoc's coverage, rather than
# fabricating a confident-looking prediction.


def build_compartment_file(deeploc_dir: Path, eggnog_genes: set[str], out_csv: Path) -> Path:
    """ORF-keyed compartment file in the tool's expected schema (columns Ids, E, M, P, O)."""
    frames = [pd.read_csv(f) for f in sorted(deeploc_dir.glob("*_deeploc_*.csv"))]
    dl = pd.concat(frames, ignore_index=True)
    rows = []
    covered = set()
    for _, r in dl.iterrows():
        gene = r["Protein_ID"]
        covered.add(gene)
        row = {"Ids": gene}
        for label, col in _DEEPLOC_DIRECT.items():
            row[label] = float(r[col])
        row["O"] = float(max(r[c] for c in _DEEPLOC_OTHER))
        rows.append(row)
    missing = eggnog_genes - covered
    for gene in sorted(missing):
        rows.append({"Ids": gene, **_NEUTRAL_ROW})
    df = pd.DataFrame(rows)[["E", "M", "P", "O", "Ids"]]
    df.to_csv(out_csv, index=False)
    print(f"compartment file: {len(covered)} genes from DeepLoc, {len(missing)} filled neutral "
          f"(no DeepLoc coverage) -> {out_csv}")
    return out_csv


def set_minimal_medium(model):
    """Reproduce the tool's own constrained-medium bound-setting (its ``carveFungi(open_bounds=False)``):
    cap all reversible bounds to +-100, then open only its 9 hardcoded minimal-medium exchanges."""
    for r in model.reactions:
        if r.lower_bound < 0:
            r.lower_bound = -100.0
        if r.upper_bound > 0:
            r.upper_bound = 100.0
    medium = ["UF03376_E", "UF02549_E", "UF03382_E", "UF03474_E", "UF02765_E",
              "UF03268_E", "UF03456_E", "UF03314_E"]
    for rid in medium:
        model.reactions.get_by_id(rid).lower_bound = -100.0
    model.reactions.get_by_id("UF03288_E").lower_bound = -10.0  # glucose, flux-capped carbon source


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--carvefungi-dir", required=True, type=Path)
    ap.add_argument("--fungi-id", default="yeast_reconstruction")
    ap.add_argument("--deeploc-dir", type=Path, default=Path("data/deeploc"))
    ap.add_argument("--eggnog-file", type=Path,
                    default=None, help="default: <carvefungi-dir>/data/annotations/functional/*.annotations")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--time-limit", type=float, default=1200.0, help="CPLEX pool-solve seconds")
    args = ap.parse_args(argv)

    cf_dir = args.carvefungi_dir.resolve()
    deeploc_dir = args.deeploc_dir.resolve()
    out_path = args.out.resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    eggnog_file = args.eggnog_file
    if eggnog_file is None:
        cands = sorted((cf_dir / "data" / "annotations" / "functional").glob("*.emapper.annotations"))
        if not cands:
            sys.exit("no eggnog annotation file found; pass --eggnog-file")
        eggnog_file = cands[0]
    eggnog_file = eggnog_file.resolve()

    import os
    os.chdir(cf_dir / "bin")
    sys.path.insert(0, str(cf_dir / "bin"))
    import CarveMeFuncPool  # noqa: E402  (needs cplex)
    import cobra  # noqa: E402
    import CreateModelEggNogPool  # noqa: E402
    import EggNogScoring  # noqa: E402

    print(f"[1/5] scoring reactions from {eggnog_file.name} ...")
    ec_to_score, ecs_genes = EggNogScoring.get_gene_to_score_dict(args.fungi_id, str(eggnog_file))
    eggnog_genes = {g for genes in ecs_genes.values() for g in genes}
    print(f"      {len(ec_to_score)} scored ECs, {len(eggnog_genes)} genes referenced")

    universal_data = EggNogScoring.assign_scores(ec_to_score)

    print("[2/5] building compartment file (this project's yeast DeepLoc, ORF-keyed) ...")
    compartment_csv = out_path.parent / f"{args.fungi_id}_compartments.csv"
    build_compartment_file(deeploc_dir, eggnog_genes, compartment_csv)
    final_scores = EggNogScoring.scoring_compartments(
        str(compartment_csv), universal_data, ecs_genes, args.fungi_id)

    print("[3/5] carving (CPLEX minmax_reduction, minimal medium, "
          f"time-limit {args.time_limit:.0f}s) ...")
    universal_model_path = cf_dir / "data" / "reactionDatabase" / "bigModelv2.21b.sbml"
    universal_model = cobra.io.read_sbml_model(str(universal_model_path))
    set_minimal_medium(universal_model)

    gap_holder: dict[str, float] = {}

    def _quiet_time_limited_pool(solver):
        cpx = solver
        cpx.parameters.mip.pool.relgap.set(0.001)  # the tool's own setting, unchanged
        cpx.parameters.timelimit.set(float(args.time_limit))
        try:
            cpx.populate_solution_pool()
        except Exception as exc:  # noqa: BLE001
            print("  populate raised:", exc, flush=True)
            return []
        try:
            gap_holder["gap"] = cpx.solution.MIP.get_mip_relative_gap()
            print(f"  CPLEX stopped at {gap_holder['gap']:.3%} gap", flush=True)
        except Exception:  # noqa: BLE001
            pass
        names = cpx.variables.get_names()
        numsol = cpx.solution.pool.get_num()
        print(f"  pool: {numsol} solutions", flush=True)
        return [CarveMeFuncPool.solution(cpx.solution.pool.get_objective_value(i),
                                         dict(zip(names, cpx.solution.pool.get_values(i), strict=True)))
                for i in range(numsol)]

    CarveMeFuncPool.generate_soln_pool = _quiet_time_limited_pool
    objective, models = CarveMeFuncPool.carve_model(
        universal_model, final_scores, eps=1e-3, min_growth=0.1, min_atpm=0.1, feast=1e-7, opti=1e-7)
    if not models:
        sys.exit("carve_model returned no solutions (CPLEX populate failed?)")
    best_i = max(range(len(objective)), key=lambda i: objective[i])
    model = models[best_i]
    print(f"      {len(models)} pool solutions; best objective {objective[best_i]:.3f} "
          f"(gap {gap_holder.get('gap', float('nan')):.1%}); "
          f"{len(model.reactions)} reactions before gene/gap-fill steps")

    print("[4/5] gene annotation + biomass/exchange completion ...")
    model = CreateModelEggNogPool.add_gene_annotation(model, ecs_genes, universal_data)
    model = CreateModelEggNogPool.add_missing_reactions(args.fungi_id, model, universal_model)

    growth = model.slim_optimize()
    n_genes = len(model.genes)
    n_gpr = sum(1 for r in model.reactions if r.gene_reaction_rule)
    print(f"[5/5] result: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites, "
          f"{n_genes} genes, {n_gpr} reactions with a GPR; growth = {growth}")
    if not growth or growth < 1e-6:
        print("WARNING: model does not grow -- inspect before using it for the benchmark.")

    cobra.io.write_sbml_model(model, str(out_path))
    print(f"written -> {out_path}")


if __name__ == "__main__":
    main()
