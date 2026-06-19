"""LP-based gap-filling: fastGapFill and swiftGapFill.

For each blocked draft reaction, solves a single LP that identifies the
minimum-L1-flux (fastLP) or maximum-non-core-flux (swiftLP) extension of
the merged draft+template model that activates the blocked reaction. The
union of all template reactions activated across all LP solves is returned.

* **fastLP** minimises the sum of |template fluxes| while forcing the
  blocked reaction to carry >= epsilon flux (FASTCORE, Vlassis 2014).
* **swiftLP** maximises the sum of template fluxes in a single LP solve
  per blocked reaction (SWIFTCORE, Tefagh & Boyd 2020). Because the
  maximisation objective is degenerate, results can differ between runs
  or solvers.

Unlike :func:`~raven_toolbox.gapfilling.fill.connect_blocked_reactions`,
this is LP-only (no MILP), so it scales better to large models. However,
the returned template set is a union over individual LP solutions and may
be larger than strictly necessary.

Metabolite identifiers in *templates* must match those in *model* (matched
by reaction id via :func:`~raven_toolbox.manipulation.transfer.add_reactions_from_model`).
"""
from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Literal

import cobra
from cobra.flux_analysis import find_blocked_reactions, flux_variability_analysis

from raven_toolbox.manipulation.transfer import add_reactions_from_model


@dataclass
class FastLPResult:
    """Outcome of an LP-based gap-fill.

    Parameters
    ----------
    added_reactions:
        Template reaction IDs selected by the LP solver (union across all
        per-blocked-reaction LP solves).
    newly_connected:
        Draft reaction IDs that were blocked in the draft and appear in the
        returned *model* (not guaranteed to carry flux until verified by FVA).
    cannot_connect:
        Blocked draft reaction IDs that remain blocked even with the full
        template database added.
    candidates_per_reaction:
        Per blocked-reaction template reaction candidates (as returned by
        each individual LP solve). Useful for greedy post-processing.
    model:
        Draft model with *added_reactions* incorporated.
    """

    added_reactions: list[str] = field(default_factory=list)
    newly_connected: list[str] = field(default_factory=list)
    cannot_connect: list[str] = field(default_factory=list)
    candidates_per_reaction: dict[str, list[str]] = field(default_factory=dict)
    model: cobra.Model = field(default_factory=cobra.Model)


def _as_models(templates: cobra.Model | Iterable[cobra.Model]) -> list[cobra.Model]:
    return [templates] if isinstance(templates, cobra.Model) else list(templates)


def _merge_templates(
    model: cobra.Model, templates: list[cobra.Model]
) -> tuple[cobra.Model, list[str]]:
    """Merge all template reactions (new ones only) into a working copy.

    Returns the merged model and the ids of reactions from templates.
    """
    working = model.copy()
    template_ids: list[str] = []
    for t in templates:
        new = [r.id for r in t.reactions if r.id not in working.reactions]
        if new:
            added = add_reactions_from_model(working, t, new, genes=False, note=None)
            template_ids += [r.id for r in added]
    return working, template_ids


def _solve_fastlp(
    working: cobra.Model,
    blocked_rid: str,
    template_ids: list[str],
    *,
    epsilon: float,
) -> list[str]:
    """Run one FASTCORE LP for a blocked reaction.

    Minimises the L1 norm (sum of |flux|) over template reactions while
    forcing *blocked_rid* to carry >= *epsilon* flux. For irreversible
    template reactions the L1 norm is simply the flux value. For reversible
    template reactions, auxiliary variables are added to linearise |v|.

    Returns the list of template reaction ids with |flux| > epsilon/2.
    """
    prob = working.problem
    aux_vars: list = []
    aux_cons: list = []
    obj_terms: list = []

    with working as m:
        blocked_rxn = m.reactions.get_by_id(blocked_rid)
        blocked_rxn.lower_bound = epsilon

        for rid in template_ids:
            rxn = m.reactions.get_by_id(rid)
            if rxn.lower_bound >= 0:
                # Irreversible — L1 = flux (already non-negative)
                obj_terms.append(rxn.flux_expression)
            else:
                # Reversible — |v| = t, with t >= v and t >= -v
                t = prob.Variable(f"_l1_{rid}", lb=0)
                c1 = prob.Constraint(
                    t - rxn.flux_expression, lb=0, name=f"_l1c1_{rid}"
                )
                c2 = prob.Constraint(
                    t + rxn.flux_expression, lb=0, name=f"_l1c2_{rid}"
                )
                aux_vars.append(t)
                aux_cons.extend([c1, c2])
                obj_terms.append(t)

        if aux_vars or aux_cons:
            m.add_cons_vars(aux_vars + aux_cons)

        if obj_terms:
            m.objective = prob.Objective(sum(obj_terms), direction="min")
        sol = m.optimize()

        if sol.status != "optimal":
            return []

        return [
            rid
            for rid in template_ids
            if abs(sol.fluxes.get(rid, 0.0)) > epsilon / 2
        ]


def _solve_swiftlp(
    working: cobra.Model,
    blocked_rid: str,
    template_ids: list[str],
    *,
    epsilon: float,
) -> list[str]:
    """Run one SWIFTCORE LP for a blocked reaction.

    Maximises the sum of template reaction fluxes (for irreversible
    templates) while forcing *blocked_rid* to carry >= *epsilon* flux.
    The LP is degenerate (many optimal solutions), so results are
    stochastic across solvers and runs.

    Returns the list of template reaction ids with flux > epsilon/2.
    """
    prob = working.problem

    with working as m:
        m.reactions.get_by_id(blocked_rid).lower_bound = epsilon

        # Maximise sum of (irreversible) template fluxes
        obj_terms = [
            m.reactions.get_by_id(rid).flux_expression
            for rid in template_ids
            if m.reactions.get_by_id(rid).lower_bound >= 0
        ]
        if obj_terms:
            m.objective = prob.Objective(sum(obj_terms), direction="max")
        sol = m.optimize()

        if sol.status != "optimal":
            return []

        return [
            rid
            for rid in template_ids
            if sol.fluxes.get(rid, 0.0) > epsilon / 2
        ]


def fill_gaps_fast_lp(
    model: cobra.Model,
    templates: cobra.Model | Iterable[cobra.Model],
    *,
    epsilon: float = 1e-4,
    variant: Literal["fast", "swift"] = "fast",
    verbose: bool = True,
) -> FastLPResult:
    """LP-based gap-filling (fastGapFill / swiftGapFill).

    For each blocked draft reaction, runs a single LP that activates that
    reaction using a minimal (or maximal) set of template reactions. The
    union of all activated template reactions is returned.

    Parameters
    ----------
    model:
        Draft model to gap-fill.
    templates:
        Universal reaction database model(s).
    epsilon:
        Minimum flux threshold for a reaction to be considered active.
    variant:
        ``"fast"`` — FASTCORE L1-norm LP (minimises L1 norm of template fluxes).
        ``"swift"`` — SWIFTCORE single LP (maximises template flux sum; faster
        but stochastic).
    verbose:
        Print progress messages.

    Returns
    -------
    FastLPResult
    """
    templates = _as_models(templates)

    # ---- Find blocked reactions in draft ----
    blocked = set(find_blocked_reactions(model))
    if not blocked:
        if verbose:
            print("fill_gaps_fast_lp: no blocked reactions found.")
        return FastLPResult(model=model.copy())

    if verbose:
        print(f"fill_gaps_fast_lp: {len(blocked)} blocked reaction(s) found.")

    # ---- Merge draft + templates ----
    working, template_ids = _merge_templates(model, templates)
    if not template_ids:
        if verbose:
            print("fill_gaps_fast_lp: no new template reactions to add.")
        return FastLPResult(cannot_connect=sorted(blocked), model=model.copy())

    # ---- Determine which blocked reactions are rescuable in merged model ----
    blocked_in_model = [r for r in blocked if r in working.reactions]
    if blocked_in_model:
        fva = flux_variability_analysis(
            working, reaction_list=blocked_in_model, fraction_of_optimum=0.0
        )
        rescuable = [
            r
            for r in blocked_in_model
            if r in fva.index and fva.at[r, "maximum"] > epsilon
        ]
    else:
        rescuable = []

    cannot_connect = sorted(blocked - set(rescuable))

    if not rescuable:
        if verbose:
            print("fill_gaps_fast_lp: no blocked reactions can be rescued.")
        return FastLPResult(cannot_connect=cannot_connect, model=model.copy())

    if verbose:
        print(
            f"fill_gaps_fast_lp: {len(rescuable)}/{len(blocked)} rescuable; "
            f"running {variant}LP..."
        )

    # ---- Run LP for each rescuable blocked reaction ----
    solve_fn = _solve_fastlp if variant == "fast" else _solve_swiftlp
    candidates_per_rxn: dict[str, list[str]] = {}
    all_added: set[str] = set()
    n_skipped = 0

    for blocked_rid in rescuable:
        active = solve_fn(working, blocked_rid, template_ids, epsilon=epsilon)
        if active:
            candidates_per_rxn[blocked_rid] = active
            all_added.update(active)
        else:
            n_skipped += 1
            candidates_per_rxn[blocked_rid] = []

    if verbose and n_skipped:
        print(f"fill_gaps_fast_lp: {n_skipped} LP(s) returned no solution.")

    if verbose:
        print(f"fill_gaps_fast_lp: adding {len(all_added)} template reaction(s).")

    # ---- Build output model ----
    filled = model.copy()
    remaining = set(all_added)
    for t in templates:
        ids_from_t = [r.id for r in t.reactions if r.id in remaining]
        if ids_from_t:
            add_reactions_from_model(filled, t, ids_from_t, genes=False, note=None)
            remaining -= set(ids_from_t)

    return FastLPResult(
        added_reactions=sorted(all_added),
        newly_connected=sorted(rescuable),
        cannot_connect=cannot_connect,
        candidates_per_reaction=candidates_per_rxn,
        model=filled,
    )
