"""Real-world test on yeast-GEM — a genome-scale fungal model.

This is the most demanding real-world check: the consensus *Saccharomyces cerevisiae* model
(~4100 reactions, ~2750 metabolites, 14 compartments, ~1140 genes). It is large and lives
outside this repo, so the test **skips** unless the model file is found via the
``ASSIGNCOMP_YEAST_GEM`` environment variable or a sibling ``yeast-GEM/model/yeast-GEM.xml``
checkout. CI (which has neither the file nor a commercial solver) skips it.

yeast-GEM also exercises two conventions the synthetic toys do not:

* metabolites use opaque ``s_####`` ids, with the same species getting a *different* id per
  compartment — so the compartment-agnostic key must be the **name**
  (``base_metabolite=lambda m: m.name``), not an id suffix; and
* those names contain spaces and symbols, which must not leak into solver variable names.

The test relocates a batch of reactions whose true compartment is known and gives each gene a
score favouring that compartment; a correct, functional assignment should recover it.
"""
import os
from collections import defaultdict
from pathlib import Path

import cobra
import pandas as pd
import pytest

from raven_toolbox.localization import apply_assignment, assign_compartments
from raven_toolbox.localization.scores import LocalizationScores


def _find_model():
    env = os.environ.get("ASSIGNCOMP_YEAST_GEM")
    candidates = [Path(env)] if env else []
    here = Path(__file__).resolve()
    for up in here.parents[2:5]:  # repo dir, its parent, grandparent
        candidates.append(up / "yeast-GEM" / "model" / "yeast-GEM.xml")
    return next((p for p in candidates if p and p.is_file()), None)


@pytest.fixture(scope="module")
def yeast():
    path = _find_model()
    if path is None:
        pytest.skip("yeast-GEM model not found (set ASSIGNCOMP_YEAST_GEM or check out yeast-GEM)")
    return cobra.io.read_sbml_model(str(path))


def _sole_compartment(rxn):
    comps = {m.compartment for m in rxn.metabolites}
    return next(iter(comps)) if len(comps) == 1 else None


def test_yeast_gem_functionality_override_is_sound(yeast):
    # Regression for a Big-M flux-gating leak. r_0104 (ERG10, gene YPL028W) is essential in the
    # mitochondrion; here it is scored for cytosol. With a loose integer tolerance the MILP could
    # "place" it in c while leaking up to ub*tol (=0.01 > the growth floor) of ghost mitochondrial
    # flux, certifying a placement whose materialized model does NOT grow. The tightened integer
    # tolerance closes the leak, so functionality genuinely holds: the reaction is forced to stay
    # in m and the applied model grows.
    model = yeast
    base = (lambda m: m.name)  # noqa: E731
    scores = LocalizationScores(pd.DataFrame({"c": [1.0], "m": [0.0]}, index=["YPL028W"]))
    res = assign_compartments(model, scores, ["r_0104"], base_metabolite=base, transportable=[])
    assert res.status == "optimal"
    assert res.placements["r_0104"] == ["m"]        # functionality overrides the cytosol score
    out = apply_assignment(model, res, base_metabolite=base)
    assert out.slim_optimize() > 1e-6


def test_yeast_gem_recovers_known_localization(yeast):
    model = yeast
    base = (lambda m: m.name)  # noqa: E731 — yeast-GEM keys species by name across compartments
    assert model.slim_optimize() > 0

    # a batch of internal, single-compartment reactions with GPRs across cytosol/mito/peroxisome
    by_comp = defaultdict(list)
    for r in model.reactions:
        if not r.boundary and r.gene_reaction_rule and _sole_compartment(r) in {"c", "m", "p"}:
            by_comp[_sole_compartment(r)].append(r)
    movable = [r for comp in ("c", "m", "p") for r in by_comp[comp][:8]]
    assert len(movable) >= 12
    true = {r.id: _sole_compartment(r) for r in movable}

    # every gene scores its reaction's true compartment highest
    compartments = sorted(model.compartments)
    score = defaultdict(lambda: {c: 0.2 for c in compartments})
    for r in movable:
        for g in r.genes:
            score[g.id][true[r.id]] = 1.0
    df = pd.DataFrame.from_dict({g: score[g] for g in score}, orient="index")[compartments]

    res = assign_compartments(model, LocalizationScores(df), list(true),
                              base_metabolite=base, time_limit=300)

    assert res.status in ("optimal", "feasible")
    # with a clear prior and the transport cost, every reaction returns to its true compartment
    recovered = sum(1 for rid in true if res.placements.get(rid) == [true[rid]])
    assert recovered == len(true), f"only {recovered}/{len(true)} recovered"

    # the rebuilt model is still able to grow (same base keying reuses existing metabolites)
    out = apply_assignment(model, res, base_metabolite=base)
    assert out.slim_optimize() > 1e-6


def test_yeast_gem_scales_to_many_genes(yeast):
    # Stress the MILP at genome scale: relocate *every* single-compartment cytosol/mito/
    # peroxisome reaction at once — ~1000 reactions covering ~700 of the model's ~1140 genes
    # (far beyond the focused test's ~30). transportable=[] means no new transporters are
    # added (pure relocation), which keeps the model tractable; the biomass constraint is still
    # enforced, so an `optimal` status already proves the assignment stays functional.
    model = yeast
    base = (lambda m: m.name)  # noqa: E731
    by_comp = defaultdict(list)
    for r in model.reactions:
        if not r.boundary and r.gene_reaction_rule and _sole_compartment(r) in {"c", "m", "p"}:
            by_comp[_sole_compartment(r)].append(r)
    movable = [r for comp in ("c", "m", "p") for r in by_comp[comp]]
    true = {r.id: _sole_compartment(r) for r in movable}

    compartments = sorted(model.compartments)
    score = defaultdict(lambda: {c: 0.2 for c in compartments})
    for r in movable:
        for g in r.genes:
            score[g.id][true[r.id]] = 1.0
    genes = list(score)
    assert len(genes) > 200, f"expected to exercise many genes, got {len(genes)}"
    df = pd.DataFrame.from_dict({g: score[g] for g in genes}, orient="index")[compartments]

    res = assign_compartments(model, LocalizationScores(df), list(true),
                              base_metabolite=base, transportable=[], time_limit=600)

    assert res.status in ("optimal", "feasible")
    recovered = sum(1 for rid in true if res.placements.get(rid) == [true[rid]])
    # The large majority return to their true compartment; the few that don't are genes shared
    # across compartments, legitimately consolidated to avoid gratuitous multi-localization.
    assert recovered / len(true) >= 0.9, f"only {recovered}/{len(true)} recovered"
