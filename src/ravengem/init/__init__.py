"""Context-specific model extraction (tINIT / ftINIT).

Phase 4c (tINIT):
* :func:`run_init` — the INIT MILP (``runINIT``).
* :func:`score_reactions_from_genes` / :func:`gene_scores_from_expression` — gene →
  reaction scoring (``scoreComplexModel`` core; RNA-seq is the common upstream).
* :func:`get_init_model` — the tINIT pipeline (``getINITModel`` core).

ftINIT (Phase 4d) follows.
"""
from ravengem.init.build import InitModelResult, get_init_model
from ravengem.init.init import InitResult, run_init
from ravengem.init.score import gene_scores_from_expression, score_reactions_from_genes

__all__ = [
    "InitModelResult",
    "InitResult",
    "gene_scores_from_expression",
    "get_init_model",
    "run_init",
    "score_reactions_from_genes",
]
