"""The same input must give the same answer twice.

Not a cross-language check, but it protects the property the parity tiers rest
on: a comparison against MATLAB means nothing if raven-toolbox itself returns a
different model each run.

Compartment placement and gap-filling were made deterministic in #76, #83 and
the assignment gap-fill rework, all without a regression test -- the property
was verified once, by hand, and nothing has guarded it since. Tie-breaking in a
MILP is easy to reintroduce accidentally (an unordered set, a dict iteration
order, an unseeded solver), and the symptom is a benchmark that quietly moves.
"""
from __future__ import annotations

import cobra
import pytest

from raven_toolbox.gapfilling import connect_blocked_reactions
from raven_toolbox.manipulation import remove_duplicate_reactions

pytestmark = pytest.mark.parity

RUNS = 3


def _toy_model() -> cobra.Model:
    """A small network with a blocked branch and a symmetric fill choice.

    B is produced but not consumed unless one of two equally-good reactions is
    added, so any tie-breaking that is not deterministic shows up here.
    """
    model = cobra.Model("toy")
    mets = {
        mid: cobra.Metabolite(mid, compartment="c") for mid in ("a", "b", "c", "d")
    }
    model.add_metabolites(list(mets.values()))

    def rxn(rid, stoich, lb=0.0, ub=1000.0):
        r = cobra.Reaction(rid, lower_bound=lb, upper_bound=ub)
        r.add_metabolites({mets[m]: v for m, v in stoich.items()})
        return r

    model.add_reactions(
        [
            rxn("EX_a", {"a": 1}, lb=-10.0),
            rxn("R1", {"a": -1, "b": 1}),
            rxn("EX_d", {"d": -1}),
        ]
    )
    model.objective = "EX_d"
    return model


def _template() -> cobra.Model:
    """Two symmetric routes from b to d, plus a decoy."""
    template = cobra.Model("template")
    mets = {
        mid: cobra.Metabolite(mid, compartment="c") for mid in ("a", "b", "c", "d")
    }
    template.add_metabolites(list(mets.values()))

    def rxn(rid, stoich):
        r = cobra.Reaction(rid, lower_bound=0.0, upper_bound=1000.0)
        r.add_metabolites({mets[m]: v for m, v in stoich.items()})
        return r

    template.add_reactions(
        [
            rxn("T_bd_1", {"b": -1, "d": 1}),
            rxn("T_bd_2", {"b": -1, "d": 1}),
            rxn("T_bc", {"b": -1, "c": 1}),
            rxn("T_cd", {"c": -1, "d": 1}),
        ]
    )
    return template


def test_gap_fill_returns_the_same_reactions_every_run():
    results = []
    for _ in range(RUNS):
        result = connect_blocked_reactions(_toy_model(), [_template()])
        results.append(tuple(sorted(result.added_reactions)))

    assert results[0], (
        "the fixture no longer requires any fill, so this test would pass "
        "whatever the tie-breaking does -- fix the fixture, not the assertion"
    )
    assert len(set(results)) == 1, (
        f"connect_blocked_reactions returned different reaction sets across "
        f"{RUNS} runs: {set(results)}"
    )


def test_duplicate_removal_keeps_the_same_reaction_every_run():
    """Which duplicate survives must not depend on iteration order."""
    kept = []
    for _ in range(RUNS):
        model = _toy_model()
        for suffix in ("dup1", "dup2"):
            clone = cobra.Reaction(f"R1_{suffix}", lower_bound=0.0, upper_bound=1000.0)
            clone.add_metabolites(
                {m: c for m, c in model.reactions.R1.metabolites.items()}
            )
            model.add_reactions([clone])
        remove_duplicate_reactions(model)
        kept.append(tuple(sorted(r.id for r in model.reactions)))

    assert all("R1_dup1" not in run or "R1_dup2" not in run for run in kept), (
        "no duplicates were removed, so this test would pass whatever the "
        "selection does -- fix the fixture, not the assertion"
    )
    assert len(set(kept)) == 1, (
        f"remove_duplicate_reactions kept different reactions across {RUNS} "
        f"runs: {set(kept)}"
    )
