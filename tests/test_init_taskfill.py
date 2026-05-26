"""Phase 4d.4: task gap-filling (fill_tasks).

Oracle: RAVEN tinitTests T0003. Remove the exchange reactions and create a gap by
deleting R7 (e[c] -> e[s]); gap-filling against the full reference must add R7 back so
the task 'make e[s] from a[s]' becomes feasible again.
"""
from tinit_oracles import make_test_model, make_test_task

from ravengem.init import TaskFillResult, fill_tasks


def _reference_without_exchanges():
    """testModel with the exchange reactions (R1, R8) removed — the gap-fill template."""
    ref = make_test_model()
    ref.remove_reactions(["R1", "R8"], remove_orphans=False)
    return ref


def test_fills_the_gap_with_r7():
    ref = _reference_without_exchanges()
    gapped = ref.copy()
    gapped.remove_reactions(["R7"], remove_orphans=False)  # the gap
    res = fill_tasks(gapped, ref, [make_test_task()])
    assert isinstance(res, TaskFillResult)
    assert res.added_reactions == ["R7"]
    assert "R7" in {r.id for r in res.model.reactions}
    assert not res.failed_tasks


def test_no_fill_when_already_feasible():
    """A model that can already do the task gets no additions."""
    ref = _reference_without_exchanges()
    res = fill_tasks(ref.copy(), ref, [make_test_task()])
    assert res.added_reactions == []


def test_should_fail_tasks_ignored():
    from ravengem.tasks import Task

    ref = _reference_without_exchanges()
    gapped = ref.copy()
    gapped.remove_reactions(["R7"], remove_orphans=False)
    sf = Task(id="sf", should_fail=True, outputs=[("e[s]", 1.0, 1.0)])
    res = fill_tasks(gapped, ref, [sf])
    assert res.added_reactions == []  # should_fail task drives no gap-filling


def test_open_exchange_does_not_short_circuit_gapfill():
    """Boundaries are closed during gap-filling, so an open exchange can't fake feasibility.

    Give the gapped model an open exchange on e[s]; without closing boundaries the task
    'produce e[s]' would look feasible (free secretion) and R7 would never be added.
    """
    import cobra

    ref = _reference_without_exchanges()
    gapped = ref.copy()
    gapped.remove_reactions(["R7"], remove_orphans=False)
    ex_es = cobra.Reaction("EX_es", lower_bound=-1000, upper_bound=1000)
    ex_es.add_metabolites({gapped.metabolites.es: -1})
    gapped.add_reactions([ex_es])  # open exchange that must be ignored
    res = fill_tasks(gapped, ref, [make_test_task()])
    assert "R7" in res.added_reactions  # gap still detected and filled


def test_prefers_cheaper_reactions_by_score():
    """When two candidates can fill a gap, the higher-scored (cheaper) one is chosen.

    Build a gap that R7 (e[c]->e[s]) OR an alternative ALT (e[c]->e[s]) can fill; give
    ALT a much better score so it is preferred.
    """
    import cobra

    ref = _reference_without_exchanges()
    alt = cobra.Reaction("ALT", lower_bound=0, upper_bound=1000)
    alt.add_metabolites({ref.metabolites.ec: -1, ref.metabolites.es: 1})  # same as R7
    alt.gene_reaction_rule = "gALT"
    ref.add_reactions([alt])
    gapped = ref.copy()
    gapped.remove_reactions(["R7", "ALT"], remove_orphans=False)
    # ALT scored high (cost low), R7 scored low (cost high) → ALT chosen.
    res = fill_tasks(gapped, ref, [make_test_task()], rxn_scores={"ALT": 5.0, "R7": -3.0})
    assert res.added_reactions == ["ALT"]
