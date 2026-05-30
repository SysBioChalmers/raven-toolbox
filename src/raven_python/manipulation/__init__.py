"""Generic cobra.Model structural transforms that cobrapy does not cover cleanly:
reaction building from equations, batch GPR / bound changes, irreversibility splitting,
isozyme expansion, compartment merge / copy, and model merging by name."""
from .add import add_reactions_from_equations
from .change import change_gene_reaction_rules, change_reaction_equations
from .expand import expand_model
from .irreversible import convert_to_irreversible
from .merge import merge_models
from .parameters import set_variance_bounds
from .remove import remove_genes, remove_metabolites
from .simplify import (
    constrain_reversible_reactions,
    find_duplicate_reactions,
    group_linear_reactions,
    remove_dead_end_reactions,
    remove_duplicate_reactions,
)
from .transfer import add_reactions_from_model
from .transport import add_transport_reactions

__all__ = [
    "add_reactions_from_equations",
    "add_reactions_from_model",
    "add_transport_reactions",
    "change_gene_reaction_rules",
    "change_reaction_equations",
    "constrain_reversible_reactions",
    "convert_to_irreversible",
    "expand_model",
    "find_duplicate_reactions",
    "group_linear_reactions",
    "merge_models",
    "remove_dead_end_reactions",
    "remove_duplicate_reactions",
    "remove_genes",
    "remove_metabolites",
    "set_variance_bounds",
]
