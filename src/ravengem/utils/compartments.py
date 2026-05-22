"""Select reactions or metabolites by compartment.

Port of RAVEN ``getMetsInComp.m`` / ``getRxnsInComp.m``.

These are thin over cobra — cobra gives ``metabolite.compartment`` and the
``reaction.compartments`` set directly — but cobra has no first-class "objects
in compartment" accessor, and the ``include_partial`` distinction (a reaction
*fully contained* in a compartment vs merely *touching* it, i.e. a transport
reaction) is the one fiddly bit worth encapsulating. Unlike RAVEN, which returns
a boolean mask + names, these return the objects (the Pythonic, directly-usable
form).
"""
from __future__ import annotations

import cobra
from cobra import Metabolite, Reaction


def _check_compartment(model: "cobra.Model", compartment: str) -> None:
    if compartment not in model.compartments:
        raise ValueError(
            f"Compartment {compartment!r} is not in the model "
            f"(have: {sorted(model.compartments)})."
        )


def get_metabolites_in_compartment(
    model: "cobra.Model", compartment: str
) -> list[Metabolite]:
    """Return the metabolites assigned to ``compartment``.

    Port of RAVEN ``getMetsInComp.m``.
    """
    _check_compartment(model, compartment)
    return [met for met in model.metabolites if met.compartment == compartment]


def get_reactions_in_compartment(
    model: "cobra.Model", compartment: str, *, include_partial: bool = False
) -> list[Reaction]:
    """Return the reactions in ``compartment``.

    Port of RAVEN ``getRxnsInComp.m``.

    Parameters
    ----------
    include_partial
        If False (default), only reactions whose metabolites are *entirely*
        within ``compartment``. If True, also include reactions that merely
        touch it (e.g. transport reactions spanning several compartments).
    """
    _check_compartment(model, compartment)
    if include_partial:
        return [rxn for rxn in model.reactions if compartment in rxn.compartments]
    return [rxn for rxn in model.reactions if rxn.compartments == {compartment}]
