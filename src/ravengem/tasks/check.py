"""Check whether a model performs a set of metabolic tasks (port of ``checkTasks``).

For each task the model is constrained by the task's allowed inputs/outputs (and
any extra reactions / bound changes), then tested for feasibility: a task *passes*
if a steady-state flux exists, unless it is marked ``should_fail`` (then it passes
iff infeasible). No cobra equivalent.

RAVEN defines inputs/outputs via a two-column metabolite RHS (``model.b``): the
net production of a metabolite, ``Sv_m``, is constrained to ``[b1, b2]`` instead of
the usual ``0``. We do the same directly on cobra's mass-balance constraints
(``model.constraints[met.id]``): an input allows net consumption (``Sv ∈ [-UB, -LB]``)
and an output allows/requires net production (``Sv ≤ UB``, and ``≥ LB`` if ``LB>0``).
Existing boundary reactions are closed first, so inputs/outputs are defined solely
by the task (RAVEN's closed-model assumption).
"""
from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import cobra
from optlang.symbolics import Zero

from ravengem.manipulation.add import add_reactions_from_equations
from ravengem.tasks.tasklist import Task, parse_task_list

_ALLMETS = "ALLMETS"
_ALLMETSIN = "ALLMETSIN"


@dataclass
class TaskResult:
    """Result of one task: ``passed`` is the verdict (accounts for ``should_fail``)."""

    id: str
    description: str
    passed: bool
    feasible: bool
    error: str | None = None


def _set_constraint_bounds(constraint, lb: float, ub: float) -> None:
    """Set an optlang constraint's bounds without a transient lb > ub."""
    if lb > constraint.ub:
        constraint.ub = ub
        constraint.lb = lb
    else:
        constraint.lb = lb
        constraint.ub = ub


def _classify(token: str) -> tuple[str, str | None]:
    """Return ``("all", None)``, ``("comp", COMP)``, or ``("met", token_upper)``."""
    upper = token.upper()
    if upper == _ALLMETS:
        return "all", None
    if upper.startswith(_ALLMETSIN + "["):
        return "comp", upper[len(_ALLMETSIN) + 1: upper.rfind("]")]
    return "met", upper


def _metabolite_bounds(
    task: Task, name_to_id: dict[str, str], comp_to_ids: dict[str, list[str]]
) -> tuple[dict[str, list[float]], list[str]]:
    """Compute ``{met_id: [lb, ub]}`` from a task's inputs/outputs (RAVEN ``b``).

    Bulk tokens (ALLMETS / ALLMETSIN) are applied before specific metabolites, as
    RAVEN does. Returns the bounds and a list of unresolved tokens (→ task error).
    """
    bounds: dict[str, list[float]] = {}
    missing: list[str] = []

    def touch(mid: str) -> list[float]:
        return bounds.setdefault(mid, [0.0, 0.0])

    for entries, is_input in ((task.inputs, True), (task.outputs, False)):
        bulk = [(t, lb, ub) for (t, lb, ub) in entries if _classify(t)[0] != "met"]
        specific = [(t, lb, ub) for (t, lb, ub) in entries if _classify(t)[0] == "met"]
        for token, lb, ub in bulk + specific:
            kind, arg = _classify(token)
            if kind == "all":
                ids = list(name_to_id.values())
            elif kind == "comp":
                ids = comp_to_ids.get(arg, [])
            else:
                mid = name_to_id.get(arg)
                if mid is None:
                    missing.append(token)
                    continue
                ids = [mid]
            for mid in ids:
                b = touch(mid)
                if is_input:
                    b[0] = -ub  # allow net consumption up to UB (RAVEN b1 = -UBin)
                    if kind == "met":
                        b[1] = -lb
                else:
                    b[1] = ub  # allow net production up to UB
                    if kind == "met" and lb > 0:
                        b[0] = lb  # require at least LB produced
    return bounds, missing


def _run_task(base: cobra.Model, task: Task, name_to_id, comp_to_ids) -> TaskResult:
    model = base.copy()
    bounds, missing = _metabolite_bounds(task, name_to_id, comp_to_ids)
    if missing:
        return TaskResult(task.id, task.description, False, False,
                          f"unknown metabolite(s): {sorted(set(missing))}")
    for mid, (lb, ub) in bounds.items():
        if (lb, ub) != (0.0, 0.0):
            _set_constraint_bounds(model.constraints[mid], lb, ub)

    if task.equations:
        specs = [
            {"id": f"TASK_TMP_{i}", "equation": equ, "bounds": (lb, ub)}
            for i, (equ, lb, ub) in enumerate(task.equations)
        ]
        add_reactions_from_equations(model, specs, mets_by="name", allow_new_mets=True)

    for rxn_id, lb, ub in task.changed:
        if rxn_id not in model.reactions:
            return TaskResult(task.id, task.description, False, False,
                              f"CHANGED RXN not in model: {rxn_id!r}")
        model.reactions.get_by_id(rxn_id).bounds = (lb, ub)

    model.objective = model.problem.Objective(Zero, direction="max")  # feasibility only
    model.slim_optimize()
    feasible = model.solver.status == "optimal"
    return TaskResult(task.id, task.description, feasible != task.should_fail, feasible)


def check_tasks(
    model: cobra.Model,
    tasks: str | Iterable[Task],
    *,
    close_boundaries: bool = True,
) -> list[TaskResult]:
    """Run a task list against ``model`` and return a :class:`TaskResult` per task.

    ``tasks`` is a parsed list of :class:`Task` or a path to a task-list file. With
    ``close_boundaries`` (default), existing exchange/sink/demand reactions are
    closed so inputs/outputs are defined purely by the tasks (as RAVEN assumes).
    """
    if isinstance(tasks, (str, bytes)) or hasattr(tasks, "__fspath__"):
        tasks = parse_task_list(tasks)
    else:
        tasks = list(tasks)

    base = model.copy()
    if close_boundaries:
        for rxn in base.boundary:
            rxn.bounds = (0.0, 0.0)
    name_to_id = {f"{m.name}[{m.compartment}]".upper(): m.id for m in base.metabolites}
    comp_to_ids: dict[str, list[str]] = {}
    for m in base.metabolites:
        comp_to_ids.setdefault((m.compartment or "").upper(), []).append(m.id)

    return [_run_task(base, task, name_to_id, comp_to_ids) for task in tasks]
