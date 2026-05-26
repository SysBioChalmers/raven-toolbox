"""Context-specific model extraction (tINIT / ftINIT).

Phase 4c (tINIT):
* :func:`run_init` — the INIT MILP (``runINIT``).
* :func:`score_reactions_from_genes` / :func:`gene_scores_from_expression` — gene →
  reaction scoring (``scoreComplexModel`` core; RNA-seq is the common upstream).
* :func:`get_init_model` — the tINIT pipeline (``getINITModel`` core).

ftINIT (Phase 4d):
* :func:`run_ftinit` — the single-step ftINIT MILP (``ftINITInternalAlg``), with
  continuous indicators for positive-score reactions (the speedup over ``run_init``).
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
