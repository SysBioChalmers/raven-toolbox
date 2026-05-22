"""Generic cobra.Model structural transforms that cobrapy does not cover cleanly.

Hosts RAVEN structural operations such as ``convertToIrrev``
(:func:`convert_to_irreversible`) and ``expandModel`` (:func:`expand_model`),
adopted from geckopy. See PLAN.md sections 2.1b and 7.
"""
from .add import add_reactions_from_equations
from .change import change_gene_reaction_rules, change_reaction_equations
from .expand import expand_model
from .irreversible import convert_to_irreversible
from .parameters import set_parameters
from .remove import remove_genes, remove_metabolites

__all__ = [
    "add_reactions_from_equations",
    "change_gene_reaction_rules",
    "change_reaction_equations",
    "convert_to_irreversible",
    "expand_model",
    "remove_genes",
    "remove_metabolites",
    "set_parameters",
]
