"""Apply named "condition" presets — bound diffs, biomass edits, etc.

A *condition* is a YAML file (or pre-parsed dict) describing a set of
deterministic modifications to apply to a cobra model: a prelude that
resets exchange bounds, edits to a cofactor pseudoreaction
(metabolite removals + automatic charge rebalancing), a per-reaction
bounds diff, and an optional biomass-stoichiometry delta.

Yeast-GEM was the first consumer; the schema is documented in
:func:`apply_condition` and is meant to be reusable for any GEM that
wants to keep its condition presets as data rather than code.
"""
from raven_toolbox.conditions.apply import (
    DEFAULT_RESET_EXCHANGES_UPPER_BOUND,
    apply_condition,
    load_condition,
    set_reaction_bounds,
)

__all__ = [
    "DEFAULT_RESET_EXCHANGES_UPPER_BOUND",
    "apply_condition",
    "load_condition",
    "set_reaction_bounds",
]
