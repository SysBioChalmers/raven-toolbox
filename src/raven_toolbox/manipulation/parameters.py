"""Set reaction bounds to a sign-aware ±% variance band around measured values, and
set exchange-flux bounds for a chosen set of metabolites.

Cobra has no idiom for the *variance band* case (e.g. "5 ± 20 %"); the other common
bound-setting cases are cobra one-liners:

* fixed lb / ub  → ``reaction.lower_bound`` / ``upper_bound`` / ``reaction.bounds``
* equality       → ``reaction.bounds = (v, v)``
* objective      → ``model.objective = {reaction: coeff}``
* unconstrained  → ``reaction.bounds = cobra.Configuration().bounds``
* a medium       → ``model.medium = {ex_id: uptake}`` (reaction-id-keyed, one bound only)

``model.medium`` doesn't cover ``set_exchange_bounds``'s case: independent lb and ub
per metabolite, a media-only compartment filter, and a check that every exchange
reaction agrees on which sign of flux is import before touching any of them.
"""
from __future__ import annotations

import warnings
from collections.abc import Iterable, Sequence

import cobra
from cobra import Metabolite, Reaction

Number = int | float


def _resolve(model: cobra.Model, reactions) -> list[Reaction]:
    if isinstance(reactions, (str, Reaction)):
        reactions = [reactions]
    out: list[Reaction] = []
    for r in reactions:
        if isinstance(r, Reaction):
            out.append(r)
        elif r in model.reactions:
            out.append(model.reactions.get_by_id(r))
        else:
            raise ValueError(f"Reaction {r!r} not found in the model.")
    return out


def _broadcast(value, n: int) -> list[float]:
    if isinstance(value, (int, float)):
        return [float(value)] * n
    vals = [float(v) for v in value]
    if len(vals) != n:
        raise ValueError(
            f"Expected 1 or {n} values to match the reactions, got {len(vals)}."
        )
    return vals


def set_variance_bounds(
    model: cobra.Model,
    reactions: str | Reaction | Iterable,
    values: Number | Sequence[Number],
    percent: Number,
) -> list[Reaction]:
    """Constrain reactions to a ``±percent/2`` band around measured values.

    For a measured value ``v`` and ``percent`` ``p``, the bounds become
    ``v * (1 - p/200) .. v * (1 + p/200)`` — i.e. ``percent`` is the *total*
    width, split half above and half below. For a negative ``v`` the two are
    swapped so that ``lb <= ub``. E.g. ``percent=5`` gives 97.5 %..102.5 % of ``v``.

    Parameters
    ----------
    reactions
        Reaction IDs or objects.
    values
        Measured value per reaction; a scalar is broadcast to all reactions.
    percent
        Total band width as a percentage.

    Returns
    -------
    list of cobra.Reaction
        The reactions affected.
    """
    rxns = _resolve(model, reactions)
    half = percent / 200.0
    for rxn, v in zip(rxns, _broadcast(values, len(rxns)), strict=True):
        lo, hi = v * (1 - half), v * (1 + half)
        rxn.bounds = (hi, lo) if v < 0 else (lo, hi)
    return rxns


def _exchange_reactions(model: cobra.Model, media_only: bool) -> list[Reaction]:
    rxns = [r for r in model.reactions if r.boundary]
    if not media_only:
        return rxns
    extracellular_ids = {
        cid for cid, name in model.compartments.items() if name.lower() == "extracellular"
    }
    if not extracellular_ids:
        raise ValueError('Could not find any compartments named "extracellular".')
    return [r for r in rxns if next(iter(r.metabolites)).compartment in extracellular_ids]


def _resolve_met(model: cobra.Model, exchanged: list[Metabolite], token) -> Metabolite | None:
    """Match ``token`` (a Metabolite, id, or name) against the *exchanged* set only —
    a metabolite that exists in the model but isn't exchanged is unused here, same as
    a metabolite that doesn't exist at all."""
    if isinstance(token, Metabolite):
        return token if token in exchanged else None
    by_id = {m.id: m for m in exchanged}
    if token in by_id:
        return by_id[token]
    by_name = {m.name.lower(): m for m in exchanged}
    return by_name.get(str(token).lower())


def set_exchange_bounds(
    model: cobra.Model,
    mets: str | Metabolite | Iterable[str | Metabolite] | None = None,
    lb: Number | Sequence[Number] | None = None,
    ub: Number | Sequence[Number] | None = None,
    *,
    close_others: bool = True,
    media_only: bool = False,
) -> list[Metabolite]:
    """Set exchange-flux bounds for a chosen set of metabolites.

    Port of RAVEN ``setExchangeBounds``. Exchange reactions are cobra's own boundary
    reactions (single metabolite); a reaction where the same metabolite is exchanged
    more than once is warned about but still handled, matching RAVEN.

    Parameters
    ----------
    mets
        Metabolite objects, ids, or names (case-insensitive), matched against the
        model's exchanged metabolites specifically — a real metabolite that isn't
        exchanged counts as unused, same as one that doesn't exist. ``None`` (default)
        means every exchanged metabolite, in which case ``lb``/``ub`` must each be a
        single value, not a sequence.
    lb, ub
        A single value applied to every metabolite in ``mets``, or one value per
        metabolite. Default to ``cobra.Configuration().bounds`` when not given.
    close_others
        Close import (only import — export is left alone) on every other exchange
        reaction not covered by ``mets``. Forced off, with a warning, if the model's
        exchange reactions don't agree on which sign of flux is import: correctly
        closing "import" would then require a different bound (lb vs ub) reaction by
        reaction, which this function does not attempt.
    media_only
        Only consider exchange reactions whose metabolite is in the compartment named
        "extracellular" (by ``model.compartments``, case-insensitive) — e.g. skips
        sink/demand reactions on an intracellular metabolite. Raises if the model has
        no compartment with that name.

    Returns
    -------
    list of cobra.Metabolite
        Metabolites in ``mets`` that were not found among the model's exchanged
        metabolites.
    """
    exch_rxns = _exchange_reactions(model, media_only)
    exchanged = list(dict.fromkeys(next(iter(r.metabolites)) for r in exch_rxns))

    signs = {next(iter(r.metabolites.values())) >= 0 for r in exch_rxns}
    if len(signs) > 1:
        warnings.warn(
            "Some exchange reactions differ in direction, and therefore have opposite "
            'meanings of lb and ub; forcing close_others=False.',
            stacklevel=2,
        )
        close_others = False
        direction = None
    else:
        direction = "forward" if True in signs else "backward"

    default_lb, default_ub = cobra.Configuration().bounds
    if lb is None:
        lb = default_lb
    if ub is None:
        ub = default_ub

    if mets is None:
        if isinstance(lb, Sequence) or isinstance(ub, Sequence):
            raise ValueError(
                "Only one upper and one lower bound may be provided if metabolites "
                "are not specified."
            )
        targets = exchanged
        lbs = _broadcast(lb, len(targets))
        ubs = _broadcast(ub, len(targets))
        unused: list[Metabolite] = []
    else:
        tokens = [mets] if isinstance(mets, (str, Metabolite)) else list(mets)
        lbs_in = _broadcast(lb, len(tokens))
        ubs_in = _broadcast(ub, len(tokens))
        targets = []
        lbs, ubs, unused = [], [], []
        for token, lo, hi in zip(tokens, lbs_in, ubs_in, strict=True):
            met = _resolve_met(model, exchanged, token)
            if met is None:
                unused.append(token)
                continue
            targets.append(met)
            lbs.append(lo)
            ubs.append(hi)

    if not targets:
        return unused

    met_to_rxns: dict[Metabolite, list[Reaction]] = {}
    for r in exch_rxns:
        met_to_rxns.setdefault(next(iter(r.metabolites)), []).append(r)
    repeated = [m for m, rs in met_to_rxns.items() if len(rs) > 1]
    if repeated:
        warnings.warn(
            "The following metabolites are involved in more than one exchange "
            f"reaction: {', '.join(m.name for m in repeated[:10])}"
            + (f", and {len(repeated) - 10} more." if len(repeated) > 10 else ""),
            stacklevel=2,
        )

    touched: set[Reaction] = set()
    for met, lo, hi in zip(targets, lbs, ubs, strict=True):
        for rxn in met_to_rxns[met]:
            rxn.lower_bound, rxn.upper_bound = lo, hi
            touched.add(rxn)

    if close_others:
        for rxn in exch_rxns:
            if rxn in touched:
                continue
            if direction == "backward":
                rxn.lower_bound = 0
            else:
                rxn.upper_bound = 0

    return unused
