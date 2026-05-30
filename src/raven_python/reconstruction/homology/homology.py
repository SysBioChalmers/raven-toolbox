"""Build a draft model from template models + homology hits.

Key behaviour:

* clear ``bidirectional`` / ``best_hits_only`` parameters control the hit-filtering
  strictness (cleaner than a single overloaded "strictness" knob);
* GPR rewriting works on cobra's AST, not regex;
* explicit ``complex_policy`` decides what happens to AND-subunits that lack an
  ortholog (drop, keep, drop-the-reaction);
* best-hit selection is bitscore-based;
* the ortholog map is a DataFrame; provenance is structured.
"""
from __future__ import annotations

import ast
import warnings
from dataclasses import dataclass, field

import cobra
import pandas as pd

from raven_python.manipulation.merge import merge_models
from raven_python.reconstruction.homology.hits import validate_hits


@dataclass
class HomologyResult:
    """Result of :func:`get_model_from_homology`.

    Attributes
    ----------
    model
        The draft ``cobra.Model``.
    gene_map
        ``{model_id: {template_gene: [new_gene, ...]}}`` ortholog mapping used.
    """

    model: cobra.Model
    gene_map: dict = field(default_factory=dict)


class _Unmapped:
    """A GPR leaf gene with no ortholog in the new organism."""

    __slots__ = ("gene",)

    def __init__(self, gene: str):
        self.gene = gene


def _rewrite_node(node, ortho: dict, policy: str, model_id: str):
    """Rewrite a GPR AST node, substituting template genes by their orthologs.

    Returns a GPR sub-expression string, ``None`` (nothing survives), or an
    ``_Unmapped`` for a bare unmapped leaf (the parent decides what to do).
    """
    if isinstance(node, ast.Name):
        new_genes = ortho.get(node.id)
        if new_genes:
            return new_genes[0] if len(new_genes) == 1 else "(" + " or ".join(new_genes) + ")"
        return _Unmapped(node.id)

    if isinstance(node, ast.BoolOp):
        children = [_rewrite_node(c, ortho, policy, model_id) for c in node.values]
        if isinstance(node.op, ast.Or):
            # An isozyme branch with no ortholog is simply absent.
            parts = [c for c in children if isinstance(c, str)]
            if not parts:
                return None
            return parts[0] if len(parts) == 1 else "(" + " or ".join(parts) + ")"
        # And: apply the complex policy to unmapped subunits.
        parts = []
        for child in children:
            if isinstance(child, str):
                parts.append(child)
            elif isinstance(child, _Unmapped):
                if policy == "flag":
                    parts.append(f"OLD_{model_id}_{child.gene}")
                elif policy == "drop":
                    return None  # incomplete complex -> reaction unsupported
                # policy == "keep": drop the unmapped subunit
            else:  # None (a dead sub-branch)
                if policy == "drop":
                    return None
        if not parts:
            return None
        return parts[0] if len(parts) == 1 else "(" + " and ".join(parts) + ")"

    return None


def _rewrite_gpr(rxn, ortho: dict, policy: str, model_id: str):
    """Return the rewritten GPR string, or None if the reaction is unsupported."""
    if not rxn.gene_reaction_rule:
        return None
    # A reaction is only transferred if at least one of its genes has an ortholog.
    if not any(g.id in ortho for g in rxn.genes):
        return None
    result = _rewrite_node(rxn.gpr.body, ortho, policy, model_id)
    if isinstance(result, str):
        return result
    return None


def _strictness_to_params(strictness, bidirectional, best_hits_only, complex_policy, map_direction):
    """Map RAVEN's strictness 1/2/3 onto the clearer parameters (compat)."""
    if strictness is None:
        return bidirectional, best_hits_only, complex_policy, map_direction
    if strictness == 1:
        return True, False, complex_policy, map_direction
    if strictness == 2:
        return False, False, complex_policy, map_direction
    if strictness == 3:
        return True, True, complex_policy, map_direction
    raise ValueError(f"strictness must be 1, 2 or 3, got {strictness}")


def _ortholog_map(
    hits, model_for, model_ids, *, bidirectional, best_hits_only, score, map_direction,
    model_genes, max_evalue, min_align_len, min_identity,
):
    """Build {model_id: {template_gene: [new_gene, ...]}} from the hits table."""
    h = hits[
        (hits.evalue <= max_evalue)
        & (hits.align_len >= min_align_len)
        & (hits.identity >= min_identity)
    ]

    if best_hits_only:
        ascending = score == "evalue"
        h = h.sort_values(score, ascending=ascending)
        h = h.groupby(["from_id", "to_id", "from_gene"], sort=False).head(1)

    # Directional views, normalised to (model_id, new_gene, template_gene).
    fwd = (
        h[h.from_id == model_for][["to_id", "from_gene", "to_gene"]]
        .rename(columns={"to_id": "model_id", "from_gene": "new_gene", "to_gene": "template_gene"})
    )
    rev = (
        h[h.to_id == model_for][["from_id", "from_gene", "to_gene"]]
        .rename(columns={"from_id": "model_id", "from_gene": "template_gene", "to_gene": "new_gene"})
    )
    fwd = fwd[fwd.model_id.isin(model_ids)]
    rev = rev[rev.model_id.isin(model_ids)]

    if bidirectional:
        pairs = fwd.merge(rev, on=["model_id", "new_gene", "template_gene"], how="inner")
    elif map_direction == "new_to_old":
        pairs = fwd
    else:
        pairs = rev
    pairs = pairs[["model_id", "new_gene", "template_gene"]].drop_duplicates()
    if pairs.empty:
        return {}

    # Keep only template genes that actually exist in their model.
    pairs = pairs[pairs.apply(lambda r: r.template_gene in model_genes.get(r.model_id, ()), axis=1)]

    ortho: dict = {}
    for model_id, template_gene, new_gene in zip(pairs.model_id, pairs.template_gene, pairs.new_gene, strict=True):
        ortho.setdefault(model_id, {}).setdefault(template_gene, [])
        if new_gene not in ortho[model_id][template_gene]:
            ortho[model_id][template_gene].append(new_gene)
    for per_model in ortho.values():
        for genes in per_model.values():
            genes.sort()
    return ortho


def _apply_preferred_order(ortho: dict, order: list[str]) -> dict:
    """Each new gene's reactions come from the first model (in order) that maps it."""
    winner: dict = {}  # new_gene -> winning model_id
    for model_id in order:
        for new_genes in ortho.get(model_id, {}).values():
            for ng in new_genes:
                winner.setdefault(ng, model_id)
    pruned: dict = {mid: {} for mid in ortho}
    for model_id, per_model in ortho.items():
        for template_gene, new_genes in per_model.items():
            kept = [ng for ng in new_genes if winner.get(ng) == model_id]
            if kept:
                pruned[model_id][template_gene] = kept
    return pruned


def get_model_from_homology(
    models,
    hits: pd.DataFrame,
    model_for: str,
    *,
    preferred_order=None,
    bidirectional: bool = True,
    best_hits_only: bool = False,
    map_direction: str = "new_to_old",
    score: str = "bitscore",
    complex_policy: str = "flag",
    only_genes_in_models: bool = False,
    max_evalue: float = 1e-30,
    min_align_len: int = 200,
    min_identity: float = 40,
    strictness: int | None = None,
) -> HomologyResult:
    """Build a draft model for ``model_for`` by transferring reactions from templates.

    ``strictness`` (1/2/3) is a legacy alias for ``bidirectional`` / ``best_hits_only``.
    """
    if isinstance(models, cobra.Model):
        models = [models]
    if complex_policy not in ("flag", "keep", "drop"):
        raise ValueError(f"complex_policy must be flag/keep/drop, got {complex_policy!r}")
    if map_direction not in ("new_to_old", "old_to_new"):
        raise ValueError(f"map_direction must be new_to_old/old_to_new, got {map_direction!r}")
    bidirectional, best_hits_only, complex_policy, map_direction = _strictness_to_params(
        strictness, bidirectional, best_hits_only, complex_policy, map_direction
    )
    validate_hits(hits)

    model_by_id = {m.id: m for m in models}
    model_ids = list(model_by_id)
    model_genes = {mid: {g.id for g in m.genes} for mid, m in model_by_id.items()}
    all_model_genes = set().union(*model_genes.values()) if model_genes else set()

    # Sanity: each template should overlap the hits by >=5% of its genes.
    for mid, genes in model_genes.items():
        in_hits = genes & (set(hits.from_gene) | set(hits.to_gene))
        if genes and len(in_hits) < 0.05 * len(genes):
            warnings.warn(
                f"<5% of genes in template '{mid}' appear in the hits table; "
                "check that the FASTA and model use the same gene identifiers.",
                stacklevel=2,
            )

    if only_genes_in_models:
        hits = hits[hits.from_gene.isin(all_model_genes) | hits.to_gene.isin(all_model_genes)]

    ortho = _ortholog_map(
        hits, model_for, model_ids, bidirectional=bidirectional, best_hits_only=best_hits_only,
        score=score, map_direction=map_direction,
        model_genes=model_genes, max_evalue=max_evalue, min_align_len=min_align_len,
        min_identity=min_identity,
    )

    order = [str(x) for x in preferred_order] if preferred_order else model_ids
    if preferred_order and len(models) > 1:
        ortho = _apply_preferred_order(ortho, order)

    # Build a per-template model holding only the transferred reactions with rewritten GPRs.
    transferred = []
    for mid in order:
        model = model_by_id.get(mid)
        if model is None:
            continue
        per_model = ortho.get(mid, {})
        m = model.copy()
        keep: dict[str, str] = {}
        for rxn in m.reactions:
            new_gpr = _rewrite_gpr(rxn, per_model, complex_policy, mid)
            if new_gpr is not None:
                keep[rxn.id] = new_gpr
        m.remove_reactions([r for r in m.reactions if r.id not in keep], remove_orphans=True)
        for rid, gpr in keep.items():
            r = m.reactions.get_by_id(rid)
            r.gene_reaction_rule = gpr
            r.notes = {"note": "Included by get_model_from_homology", "confidence_score": 2,
                       "homology_source": mid}
        if m.reactions:
            transferred.append(m)

    if transferred:
        draft = merge_models(transferred, match_by="name")
    else:
        draft = cobra.Model()
    draft.id = model_for
    draft.name = "Generated by get_model_from_homology using " + ", ".join(model_ids)

    # Drop OLD_ placeholder genes that ended up orphaned (none survive in OR branches by construction).
    orphan_genes = [g for g in draft.genes if not g.reactions]
    for g in orphan_genes:
        draft.genes.remove(g)

    return HomologyResult(model=draft, gene_map=ortho)
