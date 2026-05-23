"""Sort a model's identifiers alphabetically.

Port of RAVEN ``sortIdentifiers.m``. cobra's ``DictList.sort`` reorders one list
(and rebuilds its lookup index), but there is no single "sort the whole model"
call; this provides it. Useful for deterministic, diff-friendly output.
"""
from __future__ import annotations

import cobra


def sort_identifiers(model: cobra.Model) -> cobra.Model:
    """Sort reactions, metabolites and genes alphabetically by ID, in place.

    Returns the same (mutated) model for convenience. Compartments are a plain
    dict and are emitted sorted by writers as needed.
    """
    model.reactions.sort(key=lambda r: r.id)
    model.metabolites.sort(key=lambda m: m.id)
    model.genes.sort(key=lambda g: g.id)
    return model
