"""Reduce a model by removing/merging reactions that cannot carry flux.

Ports the cobra-absent reduction modes of RAVEN ``simplifyModel.m`` as focused
functions. The other modes are cobra one-liners and live in the migration
cheatsheet (PLAN.md §1):

* ``deleteMinMax`` (no-flux) → ``cobra.flux_analysis.find_blocked_reactions``
* ``deleteZeroInterval`` → remove reactions with ``bounds == (0, 0)`` then prune
* ``deleteUnconstrained`` → RAVEN-only ``unconstrained`` field; moot on cobra

Implemented here: ``remove_dead_end_reactions`` (``deleteInaccessible``),
``remove_duplicate_reactions`` (``deleteDuplicates``),
``constrain_reversible_reactions`` (``constrainReversible``), and
``group_linear_reactions`` (``groupLinear``).
"""
from __future__ import annotations

from typing import Iterable

import cobra
from cobra.flux_analysis import flux_variability_analysis

from ravengem.manipulation.irreversible import convert_to_irreversible


def _prune_orphan_metabolites(model: "cobra.Model") -> list[str]:
    orphans = [m for m in model.metabolites if not m.reactions]
    if orphans:
        model.remove_metabolites(orphans)
    return [m.id for m in orphans]


def _can_produce_and_consume(met) -> tuple[bool, bool]:
    """Whether the network can both produce and consume ``met`` (given directions)."""
    produce = consume = False
    for rxn in met.reactions:
        coef = rxn.get_coefficient(met)
        if coef > 0:
            produce |= rxn.upper_bound > 0
            consume |= rxn.lower_bound < 0
        elif coef < 0:
            consume |= rxn.upper_bound > 0
            produce |= rxn.lower_bound < 0
    return produce, consume


def remove_dead_end_reactions(
    model: "cobra.Model", *, reserved: Iterable[str] | None = None
) -> tuple[list[str], list[str]]:
    """Iteratively remove dead-end reactions and metabolites.

    Port of RAVEN ``simplifyModel(..., deleteInaccessible=true)``. A metabolite
    is a dead end if it participates in only one reaction, or if (accounting for
    reaction directionality) it can only be produced or only consumed — such
    metabolites cannot carry steady-state flux, so the reactions touching them
    are removed. Repeats until stable.

    Returns ``(removed_reaction_ids, removed_metabolite_ids)``.
    """
    reserved = set(reserved or [])
    removed_rxns: list[str] = []
    removed_mets: list[str] = []
    while True:
        removed_mets += _prune_orphan_metabolites(model)
        dead = [
            m
            for m in model.metabolites
            if len(m.reactions) <= 1 or not all(_can_produce_and_consume(m))
        ]
        if not dead:
            break
        rxns = {r for m in dead for r in m.reactions}
        to_delete = [r for r in rxns if r.id not in reserved]
        if not to_delete:
            break
        removed_rxns += [r.id for r in to_delete]
        model.remove_reactions(to_delete)
    return removed_rxns, removed_mets


def _signature(rxn):
    mets = frozenset((m.id, c) for m, c in rxn.metabolites.items())
    return (mets, rxn.lower_bound, rxn.upper_bound, rxn.objective_coefficient)


def remove_duplicate_reactions(
    model: "cobra.Model", *, reserved: Iterable[str] | None = None
) -> list[str]:
    """Remove all-but-one of each set of duplicate reactions.

    Port of RAVEN ``simplifyModel(..., deleteDuplicates=true)`` (``contractModel``).
    Reactions are duplicates when they have identical stoichiometry, bounds, and
    objective coefficient. One of each set is kept (reserved reactions are never
    removed). Returns the removed reaction IDs.
    """
    reserved = set(reserved or [])
    groups: dict = {}
    for rxn in model.reactions:
        groups.setdefault(_signature(rxn), []).append(rxn)

    removed: list[str] = []
    for rxns in groups.values():
        if len(rxns) <= 1:
            continue
        keep = rxns[-1]
        to_remove = [r for r in rxns if r is not keep and r.id not in reserved]
        if to_remove:
            removed += [r.id for r in to_remove]
            model.remove_reactions(to_remove)
    return removed


def constrain_reversible_reactions(
    model: "cobra.Model", *, eps: float = 1e-9
) -> list[str]:
    """Constrain reversible reactions that can only carry flux one way.

    Port of RAVEN ``simplifyModel(..., constrainReversible=true)``. Runs FVA on
    each reversible reaction; if it can only carry forward flux its lower bound
    is set to 0, and if it can only carry reverse flux it is flipped to a forward
    reaction (stoichiometry, bounds, and objective negated). Returns the changed
    reaction IDs.
    """
    revs = [r for r in model.reactions if r.lower_bound < 0 < r.upper_bound]
    if not revs:
        return []
    fva = flux_variability_analysis(model, reaction_list=revs, fraction_of_optimum=0.0)

    changed: list[str] = []
    for rxn in revs:
        lo = fva.at[rxn.id, "minimum"]
        hi = fva.at[rxn.id, "maximum"]
        min_zero, max_zero = abs(lo) < eps, abs(hi) < eps
        if min_zero == max_zero:  # both ~0 (blocked) or both nonzero (truly reversible)
            continue
        if max_zero:  # only reverse flux → flip to a forward reaction
            old_lb = rxn.lower_bound
            rxn.add_metabolites({m: -2 * c for m, c in rxn.metabolites.items()})
            rxn.bounds = (0.0, -old_lb)
            rxn.objective_coefficient = -rxn.objective_coefficient
        else:  # only forward flux
            rxn.lower_bound = 0.0
        changed.append(rxn.id)
    return changed


def group_linear_reactions(
    model: "cobra.Model", *, reserved: Iterable[str] | None = None
) -> None:
    """Merge linear (single-producer, single-consumer) reaction chains.

    Port of RAVEN ``simplifyModel(..., groupLinear=true)``. **Lossy**: gene-reaction
    associations are discarded (RAVEN does the same), since merged reactions have
    no meaningful combined GPR. The model is first made irreversible, then any
    metabolite that is produced by exactly one reaction and consumed by exactly
    one reaction is eliminated by merging the two reactions. Mutates in place.
    """
    reserved = set(reserved or [])

    # Lossy: drop all gene information.
    for rxn in model.reactions:
        rxn.gene_reaction_rule = ""
    for gene in list(model.genes):
        model.genes.remove(gene)

    convert_to_irreversible(model)

    while True:
        merged = False
        for met in list(model.metabolites):
            rxns = list(met.reactions)
            if len(rxns) != 2 or any(r.id in reserved for r in rxns):
                continue
            r1, r2 = rxns
            c1, c2 = r1.get_coefficient(met), r2.get_coefficient(met)
            if (c1 > 0) == (c2 > 0):  # need one producer and one consumer
                continue
            ratio = abs(c1 / c2)
            new_lb = max(r1.lower_bound, r2.lower_bound / ratio)
            new_ub = min(r1.upper_bound, r2.upper_bound / ratio)
            new_obj = r1.objective_coefficient + r2.objective_coefficient * ratio
            # Merge r2*ratio into r1; the shared metabolite cancels and is dropped.
            r1.add_metabolites({m: c * ratio for m, c in r2.metabolites.items()})
            model.remove_reactions([r2])
            r1.bounds = (new_lb, new_ub)
            r1.objective_coefficient = new_obj
            merged = True
            break
        if not merged:
            break
        empty = [r for r in model.reactions if not r.metabolites]
        if empty:
            model.remove_reactions(empty)
        _prune_orphan_metabolites(model)
