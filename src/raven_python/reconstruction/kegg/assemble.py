"""Shared assembly of a draft model from a KO→genes mapping.

Both KEGG runtime paths end the same way: having decided which genes belong to
which KO — from organism annotations (3b.4) or from HMM hits (3b.5) — they map
KO→reaction against the gene-free reference model, OR-join the genes into each
reaction's GPR, keep gene-backed reactions (plus spontaneous ones when allowed),
and apply the ``keep*`` quality filters. That common tail lives here.
"""
from __future__ import annotations

import cobra
import pandas as pd

_DOMAINS = {"eukaryotes", "prokaryotes"}


def flag_set(rxn_flags: pd.DataFrame | None, column: str) -> set[str]:
    """Reaction ids whose ``column`` flag is truthy (handles bool or TSV strings)."""
    if rxn_flags is None or column not in rxn_flags:
        return set()
    mask = rxn_flags[column].map(lambda v: str(v).strip().lower() in ("true", "1"))
    return set(rxn_flags.loc[mask, "reaction"])


def assemble_model_from_ko_genes(
    reference_model: cobra.Model,
    ko_reaction: pd.DataFrame,
    ko_to_genes: dict[str, list[str]],
    *,
    rxn_flags: pd.DataFrame | None = None,
    keep_spontaneous: bool = True,
    keep_undefined_stoich: bool = True,
    keep_incomplete: bool = True,
    keep_general: bool = False,
    model_id: str | None = None,
    model_name: str | None = None,
    note: str | None = None,
) -> tuple[cobra.Model, dict[str, list[str]]]:
    """Build a draft model from a ``{ko: [gene, ...]}`` assignment.

    Returns ``(model, gpr_map)`` where ``gpr_map`` is the kept reactions' gene
    lists, so callers can add gene annotations afterwards.
    """
    rxn_to_kos: dict[str, set[str]] = {}
    for ko, rid in zip(ko_reaction["ko"], ko_reaction["reaction"], strict=True):
        rxn_to_kos.setdefault(rid, set()).add(ko)

    spontaneous = flag_set(rxn_flags, "spontaneous")
    drop_if = {
        "undefined_stoich": (keep_undefined_stoich, flag_set(rxn_flags, "undefined_stoich")),
        "incomplete": (keep_incomplete, flag_set(rxn_flags, "incomplete")),
        "general": (keep_general, flag_set(rxn_flags, "general")),
    }

    gpr_map: dict[str, list[str]] = {}
    spontaneous_kept: set[str] = set()
    for rxn in reference_model.reactions:
        rid = rxn.id
        # Quality filters first: dropped even if it would have genes.
        if any(not keep_flag and rid in flagged for keep_flag, flagged in drop_if.values()):
            continue
        genes = sorted({g for ko in rxn_to_kos.get(rid, ()) for g in ko_to_genes.get(ko, ())})
        if genes:
            gpr_map[rid] = genes
        elif rid in spontaneous and keep_spontaneous:
            spontaneous_kept.add(rid)

    keep = set(gpr_map) | spontaneous_kept
    model = reference_model.copy()
    if model_id is not None:
        model.id = model_id
    if model_name is not None:
        model.name = model_name
    model.remove_reactions(
        [r for r in model.reactions if r.id not in keep], remove_orphans=True
    )
    for rid, genes in gpr_map.items():
        model.reactions.get_by_id(rid).gene_reaction_rule = " or ".join(genes)
    if note is not None:
        for rid in keep:
            model.reactions.get_by_id(rid).notes["note"] = note
    return model, gpr_map
