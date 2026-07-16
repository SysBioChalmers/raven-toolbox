"""Task gap-filling for ftINIT.

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

import warnings
from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import cobra
from cobra.exceptions import OptimizationError
from optlang.symbolics import Real, add, mul

from raven_toolbox.manipulation.boundary import close_model
from raven_toolbox.tasks import Task
from raven_toolbox.tasks.check import (
    _metabolite_bounds,
    _set_constraint_bounds,
    apply_task_constraints,
    task_name_maps,
)

_FILL_TIME_LIMIT = 300.0  # RAVEN ftINITFillGapsMILP TimeLimit (per task)
_FILL_SEED = 26           # RAVEN ftINITFillGapsMILP Seed

_DEFAULT_SCORE = -1.0   # RAVEN: missing scores default to -1 (cost 1)
_MAX_SCORE = -0.1       # RAVEN min(score, -0.1): every added reaction costs ≥ 0.1


@dataclass
class TaskFillResult:
    """Result of :func:`fill_tasks`: the gap-filled model and what was added."""

    model: cobra.Model
    added_reactions: list[str]
    failed_tasks: list[str]


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


def _add_reference_reactions(
    target: cobra.Model, reference_model: cobra.Model, rxn_ids: Iterable[str]
) -> None:
    """Add ``rxn_ids`` from ``reference_model`` to ``target`` by reconstruction.

    ``cobra.Reaction.copy()`` deep-copies through the reaction↔metabolite graph and hits the
    recursion limit on large (and, on newer Python, even small) models, so the reactions are
    rebuilt from their id/bounds/stoichiometry/GPR instead — matching RAVEN's ``addRxns``
    (which adds reactions from equations, not by copying objects). Metabolites the reactions
    need that are missing from ``target`` are created first.
    """
    have = {m.id for m in target.metabolites}
    new_mets: dict[str, cobra.Metabolite] = {}
    for rid in rxn_ids:
        for met in reference_model.reactions.get_by_id(rid).metabolites:
            if met.id not in have and met.id not in new_mets:
                new_mets[met.id] = cobra.Metabolite(
                    met.id, name=met.name, compartment=met.compartment,
                    formula=met.formula, charge=met.charge,
                )
    if new_mets:
        target.add_metabolites(list(new_mets.values()))
    built = []
    for rid in rxn_ids:
        src = reference_model.reactions.get_by_id(rid)
        new = cobra.Reaction(src.id, name=src.name,
                             lower_bound=src.lower_bound, upper_bound=src.upper_bound)
        built.append((new, src))
    target.add_reactions([new for new, _ in built])
    for new, src in built:
        new.add_metabolites({target.metabolites.get_by_id(m.id): c for m, c in src.metabolites.items()})
        new.gene_reaction_rule = src.gene_reaction_rule


def _set_fill_solver(model: cobra.Model, time_limit: float | None, seed: int) -> None:
    """Match RAVEN ``ftINITFillGapsMILP``: single-threaded, fixed seed, tight integrality.

    The seed matters: the MILP is degenerate (many equal-cost fills), and a fixed seed
    makes the chosen reactions reproducible run to run.
    """
    try:  # Gurobi-specific; harmless on other backends
        params = model.solver.problem.Params
        params.Threads = 1
        params.Seed = seed
        params.IntFeasTol = 1e-9
        if time_limit is not None:
            params.TimeLimit = float(time_limit)
    except Exception:  # noqa: BLE001
        if time_limit is not None:
            model.solver.configuration.timeout = int(time_limit)


def _gap_fill_task(
    reference_model: cobra.Model, present_ids: set[str], task: Task,
    costs: dict[str, float], *, time_limit: float | None, seed: int,
) -> list[str]:
    """Min-cost reference reactions that make ``task`` feasible (RAVEN ``ftINITFillGaps``).

    The MILP runs on a **closed copy of the reference model**, which already contains every
    candidate reaction, rather than copying candidates into the target model one at a time.
    This mirrors RAVEN's ``fullModel = tRefModel`` and is the crucial difference from a naive
    port: copying the (thousands of) candidate reactions into a genome-scale model per task
    both dominated the runtime and overflowed the recursion limit, so the growth task's fill
    never returned an incumbent within the time limit and growth was silently left broken.

    Every reference reaction not already ``present`` in the target model is gated by a binary
    (off ⇒ no flux) and its cost minimised subject to the task's ranged metabolite bounds.
    Returns the ids of the reactions the MILP turns on.
    """
    work = close_model(reference_model)  # task I/O via the task's b; holds all candidates
    candidates = [r.id for r in work.reactions if r.id not in present_ids and not r.boundary]
    if not candidates:  # nothing to add → task cannot be made feasible
        raise OptimizationError(f"gap-filling found no candidates for task {task.id!r}.")
    name_to_id, comp_to_ids = task_name_maps(work)
    _, error = apply_task_constraints(work, task, name_to_id, comp_to_ids)
    if error is not None:
        raise OptimizationError(
            f"task {task.id!r} could not be applied to the reference: {error}"
        )

    prob = work.problem
    extras = []
    objective_terms = []
    for cid in candidates:
        rxn = work.reactions.get_by_id(cid)
        y = prob.Variable(f"_fill_{cid}", type="binary")
        # off ⇒ no flux; on ⇒ the reaction's own bounds apply.
        extras += [
            y,
            prob.Constraint(rxn.flux_expression - rxn.upper_bound * y, ub=0.0, name=f"_fillub_{cid}"),
            prob.Constraint(rxn.flux_expression - rxn.lower_bound * y, lb=0.0, name=f"_filllb_{cid}"),
        ]
        objective_terms.append(mul([Real(costs[cid]), y]))
    work.add_cons_vars(extras)
    # add() over a flat list, not Python sum() — the latter is O(n²) in sympy and with
    # thousands of candidates dominates gap-fill runtime (see ftINIT/tINIT, same fix).
    work.objective = prob.Objective(add(objective_terms), direction="min")
    _set_fill_solver(work, time_limit, seed)
    work.slim_optimize()
    # Accept a near-optimal incumbent (time_limit); only a truly infeasible fill (no
    # incumbent) means the task cannot be satisfied from the reference.
    if work.solver.status not in ("optimal", "feasible", "suboptimal", "time_limit") or \
            work.variables[f"_fill_{candidates[0]}"].primal is None:
        raise OptimizationError(f"gap-filling found no way to make task {task.id!r} feasible.")
    return [cid for cid in candidates
            if (work.variables[f"_fill_{cid}"].primal or 0.0) > 0.5]


def fill_tasks(
    model: cobra.Model,
    reference_model: cobra.Model,
    tasks: Iterable[Task],
    *,
    rxn_scores: Mapping[str, float] | None = None,
    time_limit: float | None = _FILL_TIME_LIMIT,
    seed: int = _FILL_SEED,
) -> TaskFillResult:
    """Add minimum-cost reference reactions so every task is feasible in ``model``.

    Port of RAVEN ``ftINITFillGapsForAllTasks``: task by task, if a task is infeasible in
    the (growing) model, :func:`_gap_fill_task` finds the minimum-cost set of reference
    reactions that restores it and they are added, carrying forward so later tasks see the
    earlier additions. ``reference_model`` supplies the candidates (its reactions not yet in
    the model, excluding exchange/boundary reactions); ``rxn_scores`` (original reaction id →
    score) sets each candidate's cost as ``−min(score, −0.1)`` (missing → cost 1).
    ``should_fail`` tasks are ignored. Each gap-fill MILP is single-threaded with a fixed
    ``seed`` and bounded by ``time_limit`` (RAVEN's 300 s).

    Boundary reactions are closed while testing/solving each task, so task inputs and outputs
    come solely from the task's ranged metabolite bounds (RAVEN gap-fills the exchange-free
    model). The returned model keeps its boundary reactions. Tasks that could not be filled
    are returned in ``failed_tasks`` **and** raised as a warning — a non-empty list means the
    context model cannot perform those tasks, which callers should not ignore silently.
    """
    scores = dict(rxn_scores or {})
    tasks = list(tasks)

    out = model.copy()
    added: list[str] = []
    failed: list[str] = []
    for task in tasks:
        if task.should_fail:
            continue
        name_to_id, comp_to_ids = task_name_maps(out)
        if _feasible(out, task, name_to_id, comp_to_ids):
            continue
        # Candidates are the reference reactions not yet in the (growing) model.
        present = {r.id for r in out.reactions}
        costs = {r.id: -min(scores.get(r.id, _DEFAULT_SCORE), _MAX_SCORE)
                 for r in reference_model.reactions if r.id not in present and not r.boundary}
        try:
            chosen = _gap_fill_task(reference_model, present, task, costs,
                                    time_limit=time_limit, seed=seed)
        except OptimizationError:
            failed.append(task.id)
            continue
        if chosen:
            _add_reference_reactions(out, reference_model, chosen)
            added.extend(chosen)
    if failed:
        warnings.warn(
            f"fill_tasks: {len(failed)} task(s) could not be gap-filled and remain "
            f"infeasible: {sorted(set(failed))}",
            stacklevel=2,
        )
    return TaskFillResult(out, added, failed)
