"""Context-specific model extraction (tINIT / ftINIT).

tINIT:
* :func:`run_init` — the classic INIT MILP.
* :func:`score_reactions_from_genes` / :func:`gene_scores_from_expression` —
  gene → reaction scoring (RNA-seq is the common upstream).
* :func:`get_init_model` — the tINIT pipeline (dead-end removal + ``run_init``).

ftINIT (faster, staged):
* :func:`run_ftinit` — the single-step ftINIT MILP (continuous indicators for
  positive-score reactions; binaries only on negatives — the speedup over ``run_init``).
* :func:`ftinit` — the full pipeline (``prep_init_model`` → staged ``run_ftinit`` →
  ``fill_tasks`` → ``remove_low_score_genes``).
"""
from ravengem.init.build import InitModelResult, get_init_model
from ravengem.init.ftinit import FtInitResult, ftinit, run_ftinit
from ravengem.init.genes import remove_low_score_genes
from ravengem.init.init import InitResult, run_init
from ravengem.init.merge import group_rxn_scores, merge_linear
from ravengem.init.prep import PrepData, ReactionMasks, classify_reactions, prep_init_model
from ravengem.init.score import gene_scores_from_expression, score_reactions_from_genes
from ravengem.init.steps import InitStep, get_init_steps
from ravengem.init.taskfill import TaskFillResult, fill_tasks

__all__ = [
    "FtInitResult",
    "InitModelResult",
    "InitResult",
    "InitStep",
    "PrepData",
    "ReactionMasks",
    "TaskFillResult",
    "classify_reactions",
    "fill_tasks",
    "ftinit",
    "gene_scores_from_expression",
    "get_init_model",
    "get_init_steps",
    "group_rxn_scores",
    "merge_linear",
    "prep_init_model",
    "remove_low_score_genes",
    "run_ftinit",
    "run_init",
    "score_reactions_from_genes",
]
