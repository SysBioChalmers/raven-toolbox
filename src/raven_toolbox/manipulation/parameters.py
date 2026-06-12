"""Set reaction bounds to a sign-aware ±% variance band around measured values.

Cobra has no idiom for the *variance band* case (e.g. "5 ± 20 %"); the other common
bound-setting cases are cobra one-liners:

* fixed lb / ub  → ``reaction.lower_bound`` / ``upper_bound`` / ``reaction.bounds``
* equality       → ``reaction.bounds = (v, v)``
* objective      → ``model.objective = {reaction: coeff}``
* unconstrained  → ``reaction.bounds = cobra.Configuration().bounds``
"""
from __future__ import annotations

from collections.abc import Iterable, Sequence

import cobra
from cobra import Reaction

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
