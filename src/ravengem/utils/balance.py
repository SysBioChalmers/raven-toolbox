"""Check the elemental balance of reactions, distinguishing *unbalanced* from
*unknown* (missing formula).

cobra's ``reaction.check_mass_balance()`` silently treats a missing formula as
empty, so a reaction can look "unbalanced" — or even balanced — when the truth is
that the data is incomplete. This module checks for missing formulas first and
returns a graded status
per reaction (``balanced`` / ``unbalanced`` / ``unknown``) plus the element
imbalance — over a batch, as structured data.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import cobra


@dataclass(frozen=True)
class ElementalBalance:
    """Balance result for one reaction.

    Attributes
    ----------
    reaction_id
        ID of the reaction.
    status
        ``"balanced"`` — elements balance;
        ``"unbalanced"`` — they do not (see ``imbalance``);
        ``"unknown"`` — at least one metabolite has no formula, so it cannot be
        determined (cobra would silently miscount these).
    imbalance
        Element → net coefficient (products − reactants), only for
        ``"unbalanced"``; empty otherwise. Charge is not included.
    """

    reaction_id: str
    status: str
    imbalance: dict[str, float] = field(default_factory=dict)


def get_elemental_balance(
    model: cobra.Model, reactions=None
) -> list[ElementalBalance]:
    """Check whether reactions are elementally balanced.
    Parameters
    ----------
    reactions
        Reaction IDs/objects to check; default all reactions. (Boundary
        reactions exchange mass with the environment and will read as
        ``unbalanced`` — filter them out if that is not wanted.)

    Returns
    -------
    list of ElementalBalance
        One entry per checked reaction, in model order.
    """
    if reactions is None:
        rxns = list(model.reactions)
    else:
        if isinstance(reactions, (str, cobra.Reaction)):
            reactions = [reactions]
        rxns = [
            r if isinstance(r, cobra.Reaction) else model.reactions.get_by_id(r)
            for r in reactions
        ]

    results: list[ElementalBalance] = []
    for rxn in rxns:
        if any(not met.formula for met in rxn.metabolites):
            results.append(ElementalBalance(rxn.id, "unknown"))
            continue
        imbalance = {
            element: amount
            for element, amount in rxn.check_mass_balance().items()
            if element != "charge"
        }
        if imbalance:
            results.append(ElementalBalance(rxn.id, "unbalanced", imbalance))
        else:
            results.append(ElementalBalance(rxn.id, "balanced"))
    return results
