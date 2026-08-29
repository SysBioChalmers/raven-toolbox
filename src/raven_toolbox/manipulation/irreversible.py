"""Convert reversible reactions to an irreversible (forward + reverse) form.

cobrapy's own ``convert_to_irreversible`` was removed, so this is a genuine
implementation rather than a wrapper.
"""
from __future__ import annotations

import copy

import cobra


def convert_to_irreversible(model: cobra.Model) -> list[str]:
    """Split non-exchange reversible reactions into a forward + reverse pair.
    For each non-exchange reaction with ``lb < 0``:

    - The original reaction is kept as the forward direction. Its
      lower bound is clamped to 0.
    - A new reaction with the same ID plus a ``_REV`` suffix is added,
      representing the reverse direction. Its stoichiometry is the
      negation of the original, its bounds are ``(0, -original_lb)``,
      and it inherits the name (with " (reversible)" appended), the
      gene-protein rule, the subsystem, the annotations and the notes
      of the original. MATLAB's ``convertToIrrev`` copies the same
      per-reaction fields across (``eccodes``, ``rxnMiriams``,
      ``subSystems``, ``rxnNotes``, ...); dropping them would, among
      other things, leave the reverse reaction without an EC code and
      so without a kcat downstream.

    Exchange reactions (boundary reactions) are never split, regardless
    of their bounds, matching MATLAB behavior where exchange reactions
    are explicitly excluded from ``convertToIrrev``.

    Parameters
    ----------
    model
        A cobra.Model, mutated in place.

    Returns
    -------
    list of str
        Sorted IDs of newly added reverse reactions (the ones ending in
        ``_REV``). The forward reactions retain their original IDs.
    """
    reverse_rxns_to_add: list[cobra.Reaction] = []
    forward_updates: list[cobra.Reaction] = []

    for rxn in model.reactions:
        if rxn.boundary:
            continue
        if rxn.lower_bound >= 0:
            continue

        original_lb = rxn.lower_bound

        rev_rxn = cobra.Reaction(
            id=f"{rxn.id}_REV",
            name=(f"{rxn.name} (reversible)" if rxn.name else f"{rxn.id}_REV"),
        )
        rev_rxn.lower_bound = 0.0
        rev_rxn.upper_bound = -original_lb
        rev_rxn.add_metabolites({m: -c for m, c in rxn.metabolites.items()})
        rev_rxn.gene_reaction_rule = rxn.gene_reaction_rule
        rev_rxn.subsystem = rxn.subsystem
        rev_rxn.annotation = copy.deepcopy(rxn.annotation)
        rev_rxn.notes = copy.deepcopy(rxn.notes)

        reverse_rxns_to_add.append(rev_rxn)
        forward_updates.append(rxn)

    for rxn in forward_updates:
        rxn.lower_bound = 0.0

    if reverse_rxns_to_add:
        model.add_reactions(reverse_rxns_to_add)

    return sorted(r.id for r in reverse_rxns_to_add)
