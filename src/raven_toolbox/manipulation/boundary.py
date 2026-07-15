"""Close a model's exchange reactions for metabolic-task checking (RAVEN ``closeModel``)."""
from __future__ import annotations

import cobra

_UNIT_TOL = 1e-9


def _is_unit_exchange(rxn: cobra.Reaction) -> bool:
    """Whether ``rxn`` is a unit exchange/sink/demand: ``sum(|coeff|) == 1`` (RAVEN's rule)."""
    total = sum(abs(c) for c in rxn.metabolites.values())
    return abs(total - 1.0) < _UNIT_TOL


def close_model(model: cobra.Model) -> cobra.Model:
    """Return a copy with the exchange reactions closed, RAVEN ``closeModel``-style.

    RAVEN's ``closeModel`` adds a boundary metabolite to every reaction whose absolute
    stoichiometric coefficients sum to 1 — a single-metabolite exchange / sink / demand with
    unit coefficient — balancing it against a boundary metabolite so it can no longer carry
    flux. After it, a metabolic task's inputs and outputs are defined *solely* by the task
    constraints (as ``checkTasks`` assumes), not by leftover open exchanges.

    This closes exactly that set of reactions. The rule is **coefficient-aware**
    (``sum(|coeff|) == 1``) and so differs from cobra's ``Reaction.boundary`` (any
    single-metabolite reaction, regardless of coefficient): a single-metabolite reaction
    with ``|coeff| != 1`` is *not* an exchange in RAVEN's sense and is left open, while a
    reaction whose coefficients happen to sum to 1 is closed. Reproducing this rule is
    necessary for the task-essential-reaction discovery to match RAVEN — a mismatch here
    leaves different exchange bypasses open and changes which reactions are found essential.

    Modifies a copy; the original is untouched.
    """
    out = model.copy()
    close_model_in_place(out)
    return out


def close_model_in_place(model: cobra.Model) -> list[str]:
    """Close the unit exchange reactions of ``model`` in place; return their ids.

    In-place variant of :func:`close_model` for callers that already hold a working copy.
    """
    closed: list[str] = []
    for rxn in model.reactions:
        if _is_unit_exchange(rxn):
            rxn.bounds = (0.0, 0.0)
            closed.append(rxn.id)
    return closed
