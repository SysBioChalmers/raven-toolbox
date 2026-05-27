"""tINIT model building (port of the core of RAVEN ``getINITModel``, Phase 4c).

The high-level tINIT pipeline: turn expression-derived scores into reaction scores
(via the GPR), drop reactions that cannot carry flux, then run the INIT MILP to
extract a context-specific model.

Scope: this ports the algorithmic core. The HPA/array/single-cell *data ingestion*
of ``getINITModel`` (parsing expression files, thresholding, the single-cell
beta-fit) is upstream and belongs in Phase 5 ``omics`` — here you pass gene scores
(see :func:`gene_scores_from_expression`, the common RNA-seq path) or reaction
scores directly. Automatic task-essential-reaction discovery and gap-filling to make
tasks feasible are shared with ftINIT (4d) and deferred: pass ``essential_rxns``
explicitly if you have them.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import cobra
from cobra.flux_analysis import find_blocked_reactions

from ravengem.init.init import run_init
from ravengem.init.score import score_reactions_from_genes


@dataclass
class InitModelResult:
    """Result of :func:`get_init_model`."""

    model: cobra.Model
    reaction_scores: dict[str, float]
    deleted_dead_end_reactions: list[str]
    deleted_in_init: list[str]
    met_production: dict[str, bool]
    objective: float


def get_init_model(
    ref_model: cobra.Model,
    *,
    rxn_scores: Mapping[str, float] | None = None,
    gene_scores: Mapping[str, float] | None = None,
    isozyme_scoring: str = "max",
    complex_scoring: str = "min",
    no_gene_score: float = -2.0,
    essential_rxns: Iterable[str] | None = None,
    present_mets: Iterable[str] | None = None,
    prod_weight: float = 0.5,
    allow_excretion: bool = True,
    no_rev_loops: bool = False,
    remove_dead_ends: bool = True,
    eps: float = 1.0,
    big_m: float | None = None,
    mip_gap: float | None = None,
    time_limit: float | None = None,
) -> InitModelResult:
    """Extract a context-specific model with tINIT.

    Provide either ``rxn_scores`` (reaction id → score) or ``gene_scores`` (gene id →
    score, converted via the GPR with :func:`score_reactions_from_genes`). Reactions
    that cannot carry flux (with exchanges open) are removed first unless
    ``remove_dead_ends=False``; ``essential_rxns`` are kept regardless. The remaining
    model is passed to :func:`run_init`.
    """
    if (rxn_scores is None) == (gene_scores is None):
        raise ValueError("Provide exactly one of rxn_scores or gene_scores.")

    model = ref_model.copy()
    essential = set(essential_rxns or [])
    if gene_scores is not None:
        scores = score_reactions_from_genes(
            model, gene_scores, isozyme_scoring=isozyme_scoring,
            complex_scoring=complex_scoring, no_gene_score=no_gene_score,
        )
    else:
        scores = dict(rxn_scores)

    deleted_dead_end: list[str] = []
    if remove_dead_ends:
        # Determine dead ends under the same regime run_init will use: with
        # allow_excretion, any metabolite may be excreted (add a demand for each),
        # so only reactions that cannot carry flux even then are truly dead.
        probe = model.copy()
        original_ids = {r.id for r in model.reactions}
        if allow_excretion:
            has_boundary = {m.id for r in probe.boundary for m in r.metabolites}
            for met in list(probe.metabolites):
                if met.id not in has_boundary:
                    probe.add_boundary(met, type="demand")
        blocked = set(find_blocked_reactions(probe, open_exchanges=True))
        deleted_dead_end = sorted((blocked & original_ids) - essential)
        model.remove_reactions(deleted_dead_end, remove_orphans=True)

    result = run_init(
        model, scores,
        present_mets=present_mets,
        essential_rxns=essential & {r.id for r in model.reactions},
        prod_weight=prod_weight,
        allow_excretion=allow_excretion,
        no_rev_loops=no_rev_loops,
        eps=eps,
        big_m=big_m,
        mip_gap=mip_gap,
        time_limit=time_limit,
    )
    return InitModelResult(
        model=result.model,
        reaction_scores=scores,
        deleted_dead_end_reactions=deleted_dead_end,
        deleted_in_init=result.deleted_reactions,
        met_production=result.met_production,
        objective=result.objective,
    )
