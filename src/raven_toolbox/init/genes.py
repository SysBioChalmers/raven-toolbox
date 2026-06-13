"""Prune low-scoring genes from a model — the last ftINIT step.

Drop negative-scoring genes from each reaction's GPR, while
respecting enzyme structure — genes joined by **OR** (isozymes) are candidates for
removal, but at least one must remain (the least-negative if all are negative);
genes joined by **AND** (complex subunits) are *not* removed individually, though a
whole complex can be dropped as one isozyme alternative if its (aggregated) score is
negative. Operates on cobra's GPR AST recursively, so nested rules like
``G1 and (G2 or G3) and G4`` prune the inner isozyme group correctly.
"""
from __future__ import annotations

import ast
from collections.abc import Mapping

import cobra
from cobra.manipulation import remove_genes

from raven_toolbox.utils.gpr import resolve_aggregators


def _prune(node, scores, iso, cplx) -> tuple[str | None, float | None]:
    """Return (pruned GPR string, aggregate score) for an AST node, or (None, None)."""
    if isinstance(node, ast.Name):
        return node.id, scores.get(node.id)  # None = unscored (NaN: never removed)
    if not isinstance(node, ast.BoolOp):
        return None, None

    pruned = [_prune(v, scores, iso, cplx) for v in node.values]
    children: list[tuple[str, float | None]] = [(s, sc) for s, sc in pruned if s is not None]

    if isinstance(node.op, ast.And):  # complex: keep every subunit, prune nested ORs
        kept = children
    else:  # OR / isozymes: drop negative-scoring alternatives, keep at least one
        kept = [(s, sc) for s, sc in children if sc is None or sc >= 0]
        if not kept:  # all negative → keep the least-negative (every sc is non-None here)
            kept = [max(children, key=lambda c: c[1] if c[1] is not None else float("-inf"))]

    parts = [s for s, _ in kept]
    score_vals = [sc for _, sc in kept if sc is not None]
    agg = (cplx if isinstance(node.op, ast.And) else iso)
    score = agg(score_vals) if score_vals else None
    op = " and " if isinstance(node.op, ast.And) else " or "
    text = parts[0] if len(parts) == 1 else "(" + op.join(parts) + ")"
    return text, score


def remove_low_score_genes(
    model: cobra.Model,
    gene_scores: Mapping[str, float],
    *,
    isozyme_scoring: str = "max",
    complex_scoring: str = "min",
) -> tuple[cobra.Model, list[str]]:
    """Remove negative-scoring genes from GPRs (RAVEN ``removeLowScoreGenes``).

    ``gene_scores`` maps gene id → score; genes absent from it are treated as unscored
    (never removed). Returns ``(new_model, removed_gene_ids)`` — genes dropped from
    *every* rule they were in (and thus from the model). ``isozyme_scoring`` /
    ``complex_scoring`` aggregate alternative/subunit scores (``max``/``min`` default).

    When all isozyme alternatives are negative the least-negative one is kept
    **deterministically** (first on a tie), unlike RAVEN's random tie-break — same
    quality, reproducible.
    """
    iso, cplx = resolve_aggregators(isozyme_scoring, complex_scoring)

    out = model.copy()
    for rxn in out.reactions:
        body = rxn.gpr.body
        if body is None or not rxn.genes:
            continue
        pruned, _ = _prune(body, gene_scores, iso, cplx)
        if pruned is not None:
            rxn.gene_reaction_rule = pruned

    used = {g.id for rxn in out.reactions for g in rxn.genes}
    removed = sorted(g.id for g in out.genes if g.id not in used)
    if removed:
        remove_genes(out, removed, remove_reactions=False)
    return out, removed
