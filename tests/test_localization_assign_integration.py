"""Integration test on a real metabolic network.

The unit tests in ``test_assign.py`` use tiny hand-built toys to pin down specific behaviours.
This module instead exercises the MILP on a real published model — cobra's bundled E. coli core
(``textbook``: ~95 reactions, ~137 genes, real GPRs) — to confirm it builds, solves at real
scale, and returns a model that still grows.
"""
import cobra
import pandas as pd
import pytest

from raven_toolbox.localization import apply_assignment, assign_compartments
from raven_toolbox.localization.scores import LocalizationScores


@pytest.fixture(scope="module")
def textbook():
    try:
        return cobra.io.load_model("textbook")
    except Exception as exc:  # noqa: BLE001 — bundled model unavailable in this cobra build
        pytest.skip(f"cobra textbook model not available: {exc}")


def test_real_network_solves_and_stays_functional(textbook):
    model = textbook
    base = model.slim_optimize()
    assert base and base > 0

    # relocate a batch of internal, single-compartment (cytosol) reactions that carry GPRs
    movable = [
        r.id for r in model.reactions
        if not r.boundary and r.gene_reaction_rule
        and {met.compartment for met in r.metabolites} == {"c"}
    ][:12]
    assert movable, "expected some relocatable cytosolic reactions"

    genes = sorted({g.id for rid in movable for g in model.reactions.get_by_id(rid).genes})
    # synthetic localization evidence: every gene clearly prefers cytosol over a hypothetical
    # organelle 'm'. With that prior and the transport cost, the functional optimum keeps the
    # whole batch in the cytosol — and crucially the MILP must still solve on a real network.
    scores = LocalizationScores(pd.DataFrame(
        {"c": [0.9] * len(genes), "m": [0.3] * len(genes)}, index=genes))

    res = assign_compartments(model, scores, movable)

    assert res.status == "optimal"
    # score + transport cost ⇒ all relocatable reactions stay in cytosol
    assert all(res.placements[rid] == ["c"] for rid in movable)
    # the compartmentalised model still grows (essentially unchanged, nothing moved)
    out = apply_assignment(model, res)
    assert out.slim_optimize() >= 0.99 * base
