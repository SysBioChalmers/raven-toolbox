"""Set reaction bounds to a variance band around measured values.

Narrow port of the one part of RAVEN ``setParam.m`` that cobra has no idiom for:
the ``'var'`` mode. Everything else ``setParam`` did is a cobra one-liner and is
left to cobra (see the migration cheatsheet in PLAN.md §1):

* ``lb``/``ub`` → ``reaction.lower_bound`` / ``upper_bound`` / ``reaction.bounds``
* ``eq``       → ``reaction.bounds = (v, v)``
* ``obj``      → ``model.objective = {reaction: coeff}``
* ``unc``      → ``reaction.bounds = cobra.Configuration().bounds``
* batch        → a loop over ``model.reactions.get_by_any(...)``
"""
from __future__ import annotations

from typing import Iterable, Sequence, Union

import cobra
from cobra import Reaction

Number = Union[int, float]


def _resolve(model: "cobra.Model", reactions) -> list[Reaction]:
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
    model: "cobra.Model",
    reactions: Union[str, Reaction, Iterable],
    values: Union[Number, Sequence[Number]],
    percent: Number,
) -> list[Reaction]:
    """Constrain reactions to a ``±percent/2`` band around measured values.

    Port of RAVEN ``setParam(model, 'var', rxns, values, percent)``.

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
    for rxn, v in zip(rxns, _broadcast(values, len(rxns))):
        lo, hi = v * (1 - half), v * (1 + half)
        rxn.bounds = (hi, lo) if v < 0 else (lo, hi)
    return rxns
