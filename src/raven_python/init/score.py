"""Score reactions from gene scores via the GPR.

Maps per-gene scores (e.g. expression-derived: present → positive, absent → negative)
to per-reaction scores by walking each reaction's GPR: genes joined by **OR**
(isozymes) are combined with ``isozyme_scoring`` (default ``max``); genes joined by
**AND** (complexes) with ``complex_scoring`` (default ``min``). Genes missing from
``gene_scores`` are *omitted*; a reaction with no genes — or whose genes are all
missing — gets ``no_gene_score`` (default −2). These reaction scores feed
:func:`raven_python.init.run_init` and :func:`raven_python.init.ftinit`.

Upstream — the omics-data → gene-score step (thresholding, expression levels) — lives
in :mod:`raven_python.omics`; this function takes gene scores as given.
"""
from __future__ import annotations

import ast
import math
from collections.abc import Mapping

import cobra

from raven_python.utils.gpr import resolve_aggregators


def gene_scores_from_expression(
    expression: Mapping[str, float],
    reference: Mapping[str, float] | float,
    *,
    factor: float = 5.0,
    max_score: float = 10.0,
    min_score: float = -5.0,
) -> dict[str, float]:
    """Gene scores from RNA-seq/array expression, RAVEN's ``5·ln(level/reference)``.

    This is tINIT's usual entry point (RNA-seq is the common case; single-cell and
    HPA are alternative upstream sources). ``reference`` is either a per-gene
    reference level (e.g. the cross-sample mean) or a single threshold for all genes:
    a gene expressed above its reference scores positive, below it negative. The
    score is clamped to ``[min_score, max_score]``; non-positive level/reference (and
    missing reference) → ``min_score`` (RAVEN maps these NaNs to -5).
    """
    scores: dict[str, float] = {}
    for gene, level in expression.items():
        ref = reference if isinstance(reference, (int, float)) else reference.get(gene)
        if not level or not ref or level <= 0 or ref <= 0:
            scores[gene] = min_score
        else:
            scores[gene] = max(min(factor * math.log(level / ref), max_score), min_score)
    return scores


def _score_node(node, gene_scores: Mapping[str, float], iso, cplx) -> float | None:
    if isinstance(node, ast.Name):
        return gene_scores.get(node.id)  # None if the gene has no score
    if isinstance(node, ast.BoolOp):
        agg = iso if isinstance(node.op, ast.Or) else cplx
        vals = [s for v in node.values if (s := _score_node(v, gene_scores, iso, cplx)) is not None]
        return agg(vals) if vals else None
    return None


def score_reactions_from_genes(
    model: cobra.Model,
    gene_scores: Mapping[str, float],
    *,
    isozyme_scoring: str = "max",
    complex_scoring: str = "min",
    no_gene_score: float = -2.0,
) -> dict[str, float]:
    """Return ``{reaction_id: score}`` from per-gene scores via each reaction's GPR."""
    iso, cplx = resolve_aggregators(isozyme_scoring, complex_scoring)

    scores: dict[str, float] = {}
    for rxn in model.reactions:
        body = rxn.gpr.body
        if body is None or not rxn.genes:
            scores[rxn.id] = no_gene_score
        else:
            value = _score_node(body, gene_scores, iso, cplx)
            scores[rxn.id] = no_gene_score if value is None else float(value)
    return scores
