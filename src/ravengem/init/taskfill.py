"""Task gap-filling for ftINIT (port of ``ftINITFillGapsForAllTasks`` + the fill MILP).

After ftINIT extracts a context-specific model, some metabolic tasks may no longer be
feasible (the scoring removed reactions a task needs). :func:`fill_tasks` restores
feasibility by adding back the **minimum-cost** set of reactions from the reference
(template) model — cost = ``−score``, so high-scoring reactions are preferred — one
task at a time, only for tasks that are actually infeasible (a cheap LP check gates
the expensive MILP), accumulating additions across tasks.

This is a different MILP from ftINIT's main extraction: it *adds* reactions to satisfy
the task's ranged metabolite bounds (RAVEN's two-column ``b``), rather than selecting
which to keep by expression score. Exchange reactions are not used to fill gaps (task
inputs/outputs come from the task's ``b``), so they are excluded as candidates.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import cobra
from optlang.symbolics import Real, add, mul

from ravengem.tasks import Task
from ravengem.tasks.check import (
    _metabolite_bounds,
    _set_constraint_bounds,
    apply_task_constraints,
    task_name_maps,
)

_DEFAULT_SCORE = -1.0   # RAVEN: missing scores default to -1 (cost 1)
_MAX_SCORE = -0.1       # RAVEN min(score, -0.1): every added reaction costs ≥ 0.1


@dataclass
class TaskFillResult:
    """Result of :func:`fill_tasks`: the gap-filled model and what was added."""

    model: cobra.Model
    added_reactions: list[str]
    failed_tasks: list[str]


def _closed_copy(model: cobra.Model) -> cobra.Model:
    """A copy with boundary reactions closed: task I/O comes only from the task's b."""
    out = model.copy()
    for rxn in out.boundary:
        rxn.bounds = (0.0, 0.0)
    return out


def _feasible(model: cobra.Model, task: Task, name_to_id, comp_to_ids) -> bool:
    """Is ``task`` feasible in ``model`` (boundaries closed)? Tested in place, then reverted.

    Avoids copying the (genome-scale) model for each of the task list's feasibility checks
    — the copy dominated gap-fill runtime. ``with model:`` reverts the closed boundaries and
    everything ``apply_task_constraints`` does through cobra's API; the untracked direct
    metabolite mass-balance bound edits are snapshotted and restored (as in check_tasks).
    """
    bounds, missing = _metabolite_bounds(task, name_to_id, comp_to_ids)
    if missing:
        return False
    saved = {mid: (model.constraints[mid].lb, model.constraints[mid].ub) for mid in bounds}
    try:
        with model:
            for rxn in model.boundary:
                rxn.bounds = (0.0, 0.0)
            _, error = apply_task_constraints(model, task, name_to_id, comp_to_ids)
            if error is not None:
                return False
            model.slim_optimize()
            return model.solver.status == "optimal"
    finally:
        for mid, (lb, ub) in saved.items():
            _set_constraint_bounds(model.constraints[mid], lb, ub)


def _fill_one_task(
    model: cobra.Model, candidates: list[cobra.Reaction], task: Task,
    costs: dict[str, float],
) -> list[str]:
    """Min-cost set of ``candidates`` to make ``task`` feasible in ``model`` (the MILP)."""
    combined = _closed_copy(model)  # task I/O via the task's b, not the model's exchanges
    combined.add_reactions([r.copy() for r in candidates])
    name_to_id, comp_to_ids = task_name_maps(combined)
    _, error = apply_task_constraints(combined, task, name_to_id, comp_to_ids)
    if error is not None:
        raise RuntimeError(f"task {task.id!r} could not be applied to the reference: {error}")

    prob = combined.problem
    extras = []
    objective_terms = []
    for cand in candidates:
        rxn = combined.reactions.get_by_id(cand.id)
        y = prob.Variable(f"_fill_{cand.id}", type="binary")
        # off ⇒ no flux; on ⇒ the reaction's own bounds apply.
        extras += [
            y,
            prob.Constraint(rxn.flux_expression - rxn.upper_bound * y, ub=0.0,
                            name=f"_fillub_{cand.id}"),
            prob.Constraint(rxn.flux_expression - rxn.lower_bound * y, lb=0.0,
                            name=f"_filllb_{cand.id}"),
        ]
        objective_terms.append(mul([Real(costs[cand.id]), y]))
    combined.add_cons_vars(extras)
    # add() over a flat list, not Python sum() — the latter is O(n²) in sympy and with
    # thousands of candidates dominates gap-fill runtime (see ftINIT/tINIT, same fix).
    combined.objective = prob.Objective(add(objective_terms), direction="min")
    combined.slim_optimize()
    if combined.solver.status != "optimal":
        raise RuntimeError(f"gap-filling found no way to make task {task.id!r} feasible.")
    return [c.id for c in candidates
            if (combined.variables[f"_fill_{c.id}"].primal or 0.0) > 0.5]


def fill_tasks(
    model: cobra.Model,
    reference_model: cobra.Model,
    tasks: Iterable[Task],
    *,
    rxn_scores: Mapping[str, float] | None = None,
) -> TaskFillResult:
    """Add minimum-cost reference reactions so every task is feasible in ``model``.

    ``reference_model`` supplies the candidate reactions (those not already in
    ``model``, excluding exchange/boundary reactions). ``rxn_scores`` (original
    reaction id → score) sets the cost of adding each candidate as ``−min(score,
    −0.1)`` (missing → cost 1). Tasks already feasible are skipped; ``should_fail``
    tasks are ignored. The model is carried forward, so later tasks see earlier
    additions. Returns the gap-filled model and the reactions added.

    Boundary reactions are closed while testing/solving each task, so task inputs and
    outputs come solely from the task's ranged metabolite bounds (RAVEN gap-fills the
    exchange-free model). The returned model keeps its boundary reactions.
    """
    scores = dict(rxn_scores or {})
    tasks = list(tasks)
    in_model = {r.id for r in model.reactions}
    candidates = [r for r in reference_model.reactions
                  if r.id not in in_model and not r.boundary]
    costs = {r.id: -min(scores.get(r.id, _DEFAULT_SCORE), _MAX_SCORE) for r in candidates}

    out = model.copy()
    added: list[str] = []
    failed: list[str] = []
    for task in tasks:
        if task.should_fail:
            continue
        name_to_id, comp_to_ids = task_name_maps(out)
        if _feasible(out, task, name_to_id, comp_to_ids):
            continue
        # Only offer reactions not yet in the (growing) model.
        present = {r.id for r in out.reactions}
        avail = [r for r in candidates if r.id not in present]
        try:
            chosen = _fill_one_task(out, avail, task, costs)
        except RuntimeError:
            failed.append(task.id)
            continue
        if chosen:
            out.add_reactions([reference_model.reactions.get_by_id(c).copy() for c in chosen])
            added.extend(chosen)
    return TaskFillResult(out, added, failed)
