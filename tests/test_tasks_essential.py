"""Phase 4d.1: essential-reaction discovery for tasks (find_task_essential_reactions).

Oracle: RAVEN tinitTests T0002 — for testModel + the "make e[s] from a[s]" task, the
pre-merge essential reactions are R2 (the only a[s]<->a[c] link) and R7 (the only
e[c]->e[s] producer); the alternative internal paths make nothing else essential.
"""
import cobra
from tinit_oracles import (
    TEST_MODEL_TASK_ESSENTIAL_PREMERGE,
    make_test_model,
    make_test_task,
)

from ravengem.tasks import (
    EssentialReactionsResult,
    Task,
    find_task_essential_reactions,
)


def test_essential_reactions_match_oracle():
    res = find_task_essential_reactions(make_test_model(), [make_test_task()])
    assert isinstance(res, EssentialReactionsResult)
    assert sorted(res.reactions) == TEST_MODEL_TASK_ESSENTIAL_PREMERGE  # ['R2', 'R7']
    assert not res.failed_tasks


def test_essential_directions_are_forward():
    """R2 (a[s]->a[c]) and R7 (e[c]->e[s]) both carry positive flux for this task."""
    res = find_task_essential_reactions(make_test_model(), [make_test_task()])
    assert res.reactions == {"R2": 1, "R7": 1}


def test_task_metabolites_collected():
    """a[s] and e[s] are referenced by the task and must be protected from removal."""
    res = find_task_essential_reactions(make_test_model(), [make_test_task()])
    m = make_test_model()
    names = {res_id: f"{m.metabolites.get_by_id(res_id).name}"
             f"[{m.metabolites.get_by_id(res_id).compartment}]" for res_id in res.task_metabolites}
    assert set(names.values()) == {"a[s]", "e[s]"}


def test_no_task_no_essentials():
    res = find_task_essential_reactions(make_test_model(), [])
    assert res.reactions == {} and res.per_task == {}


def test_equation_metabolites_are_protected():
    """A task equation's metabolites count as task metabolites (protected from removal)."""
    m = make_test_model()
    task = Task(
        id="equ",
        inputs=[("a[s]", 0.0, 1000.0)],
        outputs=[("e[c]", 1.0, 1.0)],
        equations=[("a[c] => e[c]", 0.0, 1000.0)],  # references a[c], which is not an I/O met
    )
    res = find_task_essential_reactions(m, [task])
    names = {f"{m.metabolites.get_by_id(i).name}[{m.metabolites.get_by_id(i).compartment}]"
             for i in res.task_metabolites}
    assert {"a[c]", "e[c]"} <= names and "equ" not in res.failed_tasks


def test_infeasible_task_is_reported_failed():
    """A task requiring an impossible output is dropped, not crashed."""
    impossible = Task(id="bad", outputs=[("z[s]", 1.0, 1.0)])
    # z[s] doesn't exist -> unknown metabolite -> failed.
    res = find_task_essential_reactions(make_test_model(), [impossible])
    assert res.failed_tasks == ["bad"] and res.reactions == {}


def test_should_fail_task_defines_no_essentials():
    res = find_task_essential_reactions(
        make_test_model(), [Task(id="sf", should_fail=True, outputs=[("e[s]", 1.0, 1.0)])]
    )
    assert res.reactions == {} and "sf" not in res.per_task


def test_direction_majority_across_tasks():
    """A reaction essential reverse in two tasks and forward in one is recorded reverse."""
    # Build a tiny model where a single reaction must run in a chosen direction.
    m = cobra.Model("dir")
    a, b = (cobra.Metabolite(x, name=x, compartment="s") for x in "ab")
    m.add_metabolites([a, b])
    r = cobra.Reaction("REV", lower_bound=-1000, upper_bound=1000)
    r.add_metabolites({a: -1, b: 1})  # a <=> b
    m.add_reactions([r])
    m.objective = "REV"
    # Task forcing net production of b from a -> REV forward (+1).
    fwd = Task(id="fwd", inputs=[("a[s]", 0.0, 1000.0)], outputs=[("b[s]", 1.0, 1.0)])
    # Task forcing net production of a from b -> REV reverse (-1).
    rev = Task(id="rev", inputs=[("b[s]", 0.0, 1000.0)], outputs=[("a[s]", 1.0, 1.0)])
    res = find_task_essential_reactions(m, [rev, rev, fwd])
    assert res.reactions["REV"] == -1  # two reverse votes beat one forward
