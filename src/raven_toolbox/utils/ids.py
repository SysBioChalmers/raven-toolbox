"""Generate fresh sequential identifiers in a model's existing numbering scheme.

Ported from RAVEN MATLAB: manipulation/generateNewIds.m.
"""
from __future__ import annotations

import logging
from typing import Literal

import cobra

logger = logging.getLogger(__name__)

_ENTITY_ATTR = {"reactions": "reactions", "metabolites": "metabolites"}


def generate_new_ids(
    model: cobra.Model,
    entity_type: Literal["reactions", "metabolites"],
    prefix: str,
    *,
    quantity: int = 1,
    num_length: int = 4,
) -> list[str]:
    """Generate ``quantity`` new ids, sequentially numbered after the model's own.

    Existing ids are filtered to those starting with ``prefix``, the prefix is
    stripped, and the numeric width and starting point are both taken from
    whatever the model already has for that prefix -- matching
    generateNewIds.m exactly, quirks included: the highest existing number is
    found by sorting the stripped ids as strings, not as integers, and
    ``num_length`` is overridden by the length of that string whenever at
    least one matching id already exists. For a model whose own ids for this
    prefix are all the same width (the normal case, e.g. every id zero-padded
    to 4 digits), a string sort of zero-padded numbers agrees with a numeric
    sort; it can disagree if a prefix mixes widths (``rxn_1`` alongside
    ``rxn_0002``), the same as it would in RAVEN.

    Parameters
    ----------
    model
        A cobra model.
    entity_type
        Which id namespace to look in and generate for: ``"reactions"`` or
        ``"metabolites"``.
    prefix
        Prefix for all generated ids, e.g. ``"s_"`` or ``"r_"``.
    quantity
        Number of new ids to generate.
    num_length
        Width of the zero-padded numeric part, e.g. 4 gives ids like
        ``r_0001``. Ignored in favour of the model's own existing width
        whenever at least one id with this prefix is already present.

    Returns
    -------
    list[str]
        The generated ids, e.g. ``["r_0001", "r_0002"]``.
    """
    if entity_type not in _ENTITY_ATTR:
        raise ValueError(f"entity_type must be 'reactions' or 'metabolites', got {entity_type!r}")
    existing_ids = [e.id for e in getattr(model, _ENTITY_ATTR[entity_type])]

    stripped = sorted(eid[len(prefix):] for eid in existing_ids if eid.startswith(prefix))

    if stripped:
        last = stripped[-1]
        num_length = len(last)
        try:
            last_id = int(last)
        except ValueError:
            last_id = 0
    else:
        last_id = 0
        logger.info(
            "No %s ids with prefix %r currently exist in the model. "
            "The first new id will be %r.",
            entity_type, prefix, f"{prefix}{1:0{num_length}d}",
        )

    return [f"{prefix}{k + last_id:0{num_length}d}" for k in range(1, quantity + 1)]
