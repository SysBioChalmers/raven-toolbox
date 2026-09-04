"""Reactions whose flux changed from a reference case — port of RAVEN's ``followChanged``.

Compares two flux distributions and reports reactions whose flux both (a)
clears an absolute-value floor in either case, (b) differs from the
reference by at least a fixed absolute amount, and (c) differs by at least
a relative percentage — all three cutoffs must pass; none alone is
sufficient. Optionally restricted to reactions touching a given list of
metabolite *names* (matching ``followChanged.m``'s own name-based, not
id-based, lookup — RAVEN's own source flags this as questionable itself:
``%Should use id maybe``).
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import cobra
import pandas as pd

__all__ = ["ChangedReaction", "FollowChangedResult", "follow_changed", "print_changed_fluxes"]


@dataclass
class ChangedReaction:
    """One reaction whose flux changed, in the order it appears in the model."""

    reaction: str
    name: str
    equation: str
    flux: float
    reference_flux: float
    difference: float


@dataclass
class FollowChangedResult:
    """Outcome of a followChanged comparison.

    Parameters
    ----------
    changed:
        Reactions that passed all three cutoffs, in model order.
    missing_metabolites:
        Names from ``metabolite_list`` that matched no metabolite in the model.
    """

    changed: list[ChangedReaction]
    missing_metabolites: list[str] = field(default_factory=list)


def _flux(fluxes: pd.Series | Mapping[str, float], rxn_id: str) -> float:
    return float(fluxes[rxn_id])


def follow_changed(
    model: cobra.Model,
    fluxes_a: pd.Series | Mapping[str, float],
    fluxes_b: pd.Series | Mapping[str, float],
    *,
    cutoff_change: float = 1e-8,
    cutoff_flux: float = 1e-8,
    cutoff_diff: float = 1e-8,
    metabolite_list: Sequence[str] | None = None,
) -> FollowChangedResult:
    """Find reactions whose flux differs meaningfully between two solutions.

    Parameters
    ----------
    model:
        The model both flux vectors are for.
    fluxes_a:
        Flux for the test case, indexed by reaction id.
    fluxes_b:
        Flux for the reference case, indexed by reaction id.
    cutoff_change:
        Reactions whose flux differs by less than this many *percent*
        (relative to ``fluxes_a``) are excluded. Only applied when
        ``fluxes_a`` is nonzero for that reaction; a reaction with zero flux
        in ``fluxes_a`` is judged solely by ``cutoff_flux`` instead, since a
        percent change from zero is undefined.
    cutoff_flux:
        Reactions where ``|fluxes_a|`` and ``|fluxes_b|`` are *both* below
        this are excluded, regardless of the other cutoffs.
    cutoff_diff:
        Reactions where ``|fluxes_a - fluxes_b|`` is below this are excluded,
        regardless of the other cutoffs.
    metabolite_list:
        Metabolite *names* (case-insensitive); if given, only reactions
        touching at least one of them are considered. A name matching no
        metabolite is reported in the result rather than raising.

    Returns
    -------
    FollowChangedResult
    """
    missing: list[str] = []
    if metabolite_list is not None:
        touched: set[str] = set()
        for name in metabolite_list:
            mets = [m for m in model.metabolites if m.name.strip().lower() == name.strip().lower()]
            if not mets:
                missing.append(name)
                continue
            for met in mets:
                touched.update(r.id for r in met.reactions)
        candidates = [r for r in model.reactions if r.id in touched]
    else:
        candidates = list(model.reactions)

    quota = 1 + cutoff_change / 100
    changed: list[ChangedReaction] = []
    for rxn in candidates:
        a = _flux(fluxes_a, rxn.id)
        b = _flux(fluxes_b, rxn.id)
        if abs(a) < cutoff_flux and abs(b) < cutoff_flux:
            continue
        diff = a - b
        if abs(diff) < cutoff_diff:
            continue

        if a != 0:
            ratio = b / a
            is_changed = ratio >= quota or ratio < 1.0 / quota
        else:
            is_changed = abs(b) >= cutoff_flux

        if is_changed:
            changed.append(
                ChangedReaction(
                    reaction=rxn.id,
                    name=rxn.name,
                    equation=rxn.build_reaction_string(use_metabolite_names=True),
                    flux=a,
                    reference_flux=b,
                    difference=diff,
                )
            )

    return FollowChangedResult(changed=changed, missing_metabolites=missing)


def print_changed_fluxes(
    model: cobra.Model,
    fluxes_a: pd.Series | Mapping[str, float],
    fluxes_b: pd.Series | Mapping[str, float],
    *,
    print_fn=print,
    **kwargs,
) -> FollowChangedResult:
    """Print :func:`follow_changed`'s result and return it.

    A simpler, Pythonic rendering of ``followChanged.m``'s console report —
    not a literal reproduction of its text (see the module docstring's
    parity note on text formatting in general); the structured result is
    what's compared for parity.
    """
    result = follow_changed(model, fluxes_a, fluxes_b, **kwargs)
    for name in result.missing_metabolites:
        print_fn(f"Could not find any reactions with the metabolite {name}")
    print_fn(f"{len(result.changed)} reaction(s) changed:\n")
    for c in result.changed:
        print_fn(f"{c.reaction}: {c.equation}\n\t{c.name}")
        print_fn(f"\tFlux: {c.flux:g}  Reference flux: {c.reference_flux:g}  Difference: {c.difference:g}\n")
    return result
