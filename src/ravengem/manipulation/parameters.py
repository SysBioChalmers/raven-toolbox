"""Batch-set reaction bounds and objective coefficients.

Port of RAVEN ``setParam.m``.

cobra exposes ``reaction.bounds`` / ``lower_bound`` / ``upper_bound`` and
``model.objective`` per reaction, so the value here is (a) batching over many
reactions in one call and (b) the two modes cobra has no shorthand for: a
``var`` percentage band around a measured value, and ``reset`` to the default
bounds. RAVEN's ``paramType`` string (``'lb'``/``'ub'``/``'eq'``/``'obj'``/
``'var'``/``'unc'``) becomes readable keywords; RAVEN's ``'rev'`` is dropped
(it set a RAVEN-only field that has no cobra meaning).
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


def set_parameters(
    model: "cobra.Model",
    reactions: Union[str, Reaction, Iterable],
    *,
    lb: Union[Number, Sequence[Number], None] = None,
    ub: Union[Number, Sequence[Number], None] = None,
    eq: Union[Number, Sequence[Number], None] = None,
    objective: Union[Number, Sequence[Number], None] = None,
    var: tuple[Union[Number, Sequence[Number]], Number] | None = None,
    reset: bool = False,
) -> list[Reaction]:
    """Set bounds / objective for a batch of reactions.

    Port of RAVEN ``setParam.m``. Each value may be a scalar (broadcast to all
    reactions) or a sequence matching ``reactions``.

    Parameters
    ----------
    reactions
        Reaction IDs or objects.
    lb, ub
        Lower / upper bound.
    eq
        Set lower **and** upper bound to this value (equality constraint).
    objective
        Objective coefficient. Setting this **replaces** the model objective
        (all other coefficients are zeroed first), matching RAVEN.
    var
        ``(values, percent)``: set bounds to a band around each measured value,
        i.e. ``value * (1 ± percent/200)`` (swapped for negative values). E.g.
        ``var=(v, 5)`` gives 97.5 %..102.5 % of ``v``.
    reset
        Reset bounds to the cobra configuration defaults (RAVEN ``'unc'``).

    Returns
    -------
    list of cobra.Reaction
        The reactions affected.
    """
    rxns = _resolve(model, reactions)
    n = len(rxns)

    if reset:
        cfg = cobra.Configuration()
        for r in rxns:
            r.bounds = (cfg.lower_bound, cfg.upper_bound)

    if var is not None:
        values, percent = var
        half = percent / 200.0
        for r, v in zip(rxns, _broadcast(values, n)):
            r.bounds = (v * (1 + half), v * (1 - half)) if v < 0 else (
                v * (1 - half),
                v * (1 + half),
            )

    if eq is not None:
        for r, v in zip(rxns, _broadcast(eq, n)):
            r.bounds = (v, v)

    if lb is not None and ub is not None:
        for r, lo, hi in zip(rxns, _broadcast(lb, n), _broadcast(ub, n)):
            r.bounds = (lo, hi)
    elif lb is not None:
        for r, lo in zip(rxns, _broadcast(lb, n)):
            r.lower_bound = lo
    elif ub is not None:
        for r, hi in zip(rxns, _broadcast(ub, n)):
            r.upper_bound = hi

    if objective is not None:
        for r in model.reactions:
            r.objective_coefficient = 0.0
        for r, c in zip(rxns, _broadcast(objective, n)):
            r.objective_coefficient = c

    invalid = [r.id for r in rxns if r.lower_bound > r.upper_bound]
    if invalid:
        raise ValueError(f"Invalid bounds (lb > ub) for reaction(s): {invalid}")

    return rxns
