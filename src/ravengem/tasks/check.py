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
from cobra.exceptions import OptimizationError
from cobra.flux_analysis import flux_variability_analysis, pfba
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


def _build_task_model(
    base: cobra.Model, task: Task, name_to_id, comp_to_ids
) -> tuple[cobra.Model | None, set[str], str | None]:
    """Apply a task's constraints to a copy of ``base`` (feasibility objective set).

    Returns ``(model, task_metabolite_ids, error)``. ``task_metabolite_ids`` are the
    model metabolites the task references (inputs/outputs + equation mets present in
    the model) — RAVEN's ``essentialMetsForTasks``, to be protected from removal.
    ``model``/``error`` are mutually exclusive.
    """
    model = base.copy()
    bounds, missing = _metabolite_bounds(task, name_to_id, comp_to_ids)
    if missing:
        return None, set(), f"unknown metabolite(s): {sorted(set(missing))}"
    task_mets = {mid for mid in bounds}
    for mid, (lb, ub) in bounds.items():
        if (lb, ub) != (0.0, 0.0):
            _set_constraint_bounds(model.constraints[mid], lb, ub)

    if task.equations:
        existing = {m.id for m in model.metabolites}
        specs = [
            {"id": f"TASK_TMP_{i}", "equation": equ, "bounds": (lb, ub)}
            for i, (equ, lb, ub) in enumerate(task.equations)
        ]
        add_reactions_from_equations(model, specs, mets_by="name", allow_new_mets=True)
        for i in range(len(specs)):
            tmp = model.reactions.get_by_id(f"TASK_TMP_{i}")
            task_mets |= {m.id for m in tmp.metabolites if m.id in existing}

    for rxn_id, lb, ub in task.changed:
        if rxn_id not in model.reactions:
            return None, set(), f"CHANGED RXN not in model: {rxn_id!r}"
        model.reactions.get_by_id(rxn_id).bounds = (lb, ub)

    model.objective = model.problem.Objective(Zero, direction="max")  # feasibility only
    return model, task_mets, None


def _run_task(base: cobra.Model, task: Task, name_to_id, comp_to_ids) -> TaskResult:
    model, _, error = _build_task_model(base, task, name_to_id, comp_to_ids)
    if error is not None:
        return TaskResult(task.id, task.description, False, False, error)
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
    tasks = _as_tasks(tasks)
    base, name_to_id, comp_to_ids = _prepare_base(model, close_boundaries)
    return [_run_task(base, task, name_to_id, comp_to_ids) for task in tasks]


def _as_tasks(tasks: str | Iterable[Task]) -> list[Task]:
    if isinstance(tasks, (str, bytes)) or hasattr(tasks, "__fspath__"):
        return parse_task_list(tasks)
    return list(tasks)


def _prepare_base(model: cobra.Model, close_boundaries: bool):
    base = model.copy()
    if close_boundaries:
        for rxn in base.boundary:
            rxn.bounds = (0.0, 0.0)
    name_to_id = {f"{m.name}[{m.compartment}]".upper(): m.id for m in base.metabolites}
    comp_to_ids: dict[str, list[str]] = {}
    for m in base.metabolites:
        comp_to_ids.setdefault((m.compartment or "").upper(), []).append(m.id)
    return base, name_to_id, comp_to_ids


@dataclass
class EssentialReactionsResult:
    """Reactions a model *must* use to perform a task list (RAVEN ``essentialRxns``).

    ``reactions`` maps reaction id → forced flux direction (``+1`` forward, ``-1``
    reverse): the reaction must carry flux of that sign in every feasible solution of
    at least one task. ``per_task`` is the same, split by task id. ``task_metabolites``
    are the model metabolites the tasks reference (RAVEN ``essentialMetsForTasks``,
    protected from removal). ``failed_tasks`` are tasks that were infeasible or
    malformed and thus skipped (RAVEN drops these from the task list).
    """

    reactions: dict[str, int]
    per_task: dict[str, dict[str, int]]
    task_metabolites: set[str]
    failed_tasks: list[str]


def _task_essential_reactions(
    task_model: cobra.Model, candidates: list[str], tol: float
) -> dict[str, int]:
    """Reactions in ``candidates`` forced to carry flux, with direction, via FVA.

    A reaction is *essential* for the task iff zero is not attainable in any feasible
    solution — i.e. its FVA range excludes 0. This is exactly RAVEN's
    "constrain to 0 → infeasible" definition, but obtained from FVA ranges (no
    per-reaction knockout loop). The nonzero side of the range gives the forced
    direction. FVA is restricted to ``candidates`` — the reactions carrying flux in a
    minimal feasible solution, the only ones that *can* be essential (an essential
    reaction is nonzero in every feasible solution, so also in that one) — which keeps
    this cheap on genome-scale templates instead of ranging all reactions.
    """
    if not candidates:
        return {}
    fva = flux_variability_analysis(task_model, reaction_list=candidates, fraction_of_optimum=0.0)
    essential: dict[str, int] = {}
    for rxn_id, lo, hi in zip(fva.index, fva["minimum"], fva["maximum"], strict=True):
        if lo > tol:
            essential[rxn_id] = 1
        elif hi < -tol:
            essential[rxn_id] = -1
    return essential


def find_task_essential_reactions(
    model: cobra.Model,
    tasks: str | Iterable[Task],
    *,
    close_boundaries: bool = True,
    tol: float = 1e-8,
) -> EssentialReactionsResult:
    """Find the reactions a model must use to satisfy a task list.

    For each task the model is constrained as in :func:`check_tasks`, then FVA
    identifies reactions whose flux can never be zero (essential) and their forced
    direction. This is the ``prepINITModel`` step that feeds (ft)INIT: essential
    reactions are kept regardless of expression score and made irreversible in their
    forced direction. When a reaction is essential in several tasks with conflicting
    directions, the majority wins (ties → forward), matching RAVEN's ``pos < neg``.
    """
    tasks = _as_tasks(tasks)
    base, name_to_id, comp_to_ids = _prepare_base(model, close_boundaries)
    original_ids = {r.id for r in base.reactions}

    per_task: dict[str, dict[str, int]] = {}
    task_metabolites: set[str] = set()
    failed: list[str] = []
    direction_votes: dict[str, int] = {}

    for task in tasks:
        if task.should_fail:
            continue  # a task meant to fail defines no essential reactions
        task_model, task_mets, error = _build_task_model(base, task, name_to_id, comp_to_ids)
        if error is not None:
            failed.append(task.id)
            continue
        # One min-flux solve both proves feasibility and yields the essential-reaction
        # candidates (the original reactions carrying flux in a sparse solution).
        try:
            fluxes = pfba(task_model).fluxes
        except OptimizationError:
            failed.append(task.id)
            continue
        candidates = [rid for rid in original_ids if abs(fluxes.get(rid, 0.0)) > tol]
        task_metabolites |= task_mets
        essential = _task_essential_reactions(task_model, candidates, tol)
        per_task[task.id] = essential
        for rxn_id, direction in essential.items():
            direction_votes[rxn_id] = direction_votes.get(rxn_id, 0) + direction

    # Majority direction; tie (sum == 0) → forward, as RAVEN's `pos < neg`.
    reactions = {rid: (-1 if votes < 0 else 1) for rid, votes in direction_votes.items()}
    return EssentialReactionsResult(reactions, per_task, task_metabolites, failed)
