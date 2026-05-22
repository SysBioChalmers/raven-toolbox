"""Remove reactions, metabolites, or genes from a model.

Port of RAVEN ``removeReactions.m`` / ``removeMets.m`` / ``removeGenes.m``.

cobra already covers the core of all three, so this module delegates to it and
adds only the genuinely cobra-absent behaviour:

* ``remove_reactions`` — cobra's ``Model.remove_reactions(remove_orphans=True)``
  drops orphaned metabolites **and** genes together; RAVEN exposes *separable*
  flags. We keep them separable (``remove_orphan_metabolites`` /
  ``remove_orphan_genes`` independent), plus the cobra-trivial reaction removal.
* ``remove_metabolites`` — cobra matches metabolites only by ID; RAVEN's
  ``isNames`` lets you delete a metabolite in **every compartment at once** by
  name. That name resolution is the value here.
* ``remove_genes`` — cobra's ``cobra.manipulation.remove_genes`` already rewrites
  GPRs through the boolean AST (correctly: removing one gene of ``A and B``
  empties the rule, of ``A or B`` keeps the other), which is what RAVEN does via
  ``eval``. The gap is RAVEN's default of **constraining** flux-blocked reactions
  to zero rather than deleting them — gene-knockout semantics. We expose a
  ``blocked_reactions`` policy: ``"remove"``, ``"constrain"``, or ``"keep"``.
"""
from __future__ import annotations

from typing import Iterable, Union

import cobra
from cobra import Gene, Metabolite, Reaction
from cobra.manipulation import remove_genes as _cobra_remove_genes


def _as_list(obj) -> list:
    if isinstance(obj, (str, Reaction, Metabolite, Gene)):
        return [obj]
    return list(obj)


def remove_reactions(
    model: "cobra.Model",
    reactions: Union[str, Reaction, Iterable],
    *,
    remove_orphan_metabolites: bool = False,
    remove_orphan_genes: bool = False,
) -> None:
    """Remove reactions, with *separable* orphan cleanup.

    Port of RAVEN ``removeReactions.m``. Unlike cobra's coupled
    ``remove_orphans``, the metabolite and gene cleanups are independent flags.
    """
    model.remove_reactions(_as_list(reactions), remove_orphans=False)

    if remove_orphan_metabolites:
        orphan_mets = [m for m in model.metabolites if not m.reactions]
        if orphan_mets:
            model.remove_metabolites(orphan_mets)

    if remove_orphan_genes:
        orphan_genes = [g for g in model.genes if not g.reactions]
        for gene in orphan_genes:
            model.genes.remove(gene)


def remove_metabolites(
    model: "cobra.Model",
    metabolites: Union[str, Metabolite, Iterable],
    *,
    by_name: bool = False,
    destructive: bool = False,
) -> None:
    """Remove metabolites, optionally matching by name across all compartments.

    Port of RAVEN ``removeMets.m``.

    Parameters
    ----------
    by_name
        If True, ``metabolites`` are metabolite *names*; every metabolite with a
        matching name is removed, regardless of compartment (RAVEN ``isNames``).
        If False, they are IDs/objects, resolved via cobra.
    destructive
        Passed to cobra: if True, also remove every reaction the metabolite
        participates in (RAVEN's ``removeUnusedRxns`` is similar but only drops
        reactions left empty — use cobra's ``prune_unused_reactions`` for that).
    """
    if by_name:
        wanted = set(_as_list(metabolites))
        targets = [m for m in model.metabolites if m.name in wanted]
    else:
        targets = model.metabolites.get_by_any(_as_list(metabolites))
    if targets:
        model.remove_metabolites(targets, destructive=destructive)


def remove_genes(
    model: "cobra.Model",
    genes: Union[str, Gene, Iterable],
    *,
    blocked_reactions: str = "remove",
    remove_orphan_metabolites: bool = False,
) -> list[str]:
    """Remove genes and handle reactions left unable to carry flux.

    Port of RAVEN ``removeGenes.m``. GPR rewriting (with correct AND/OR
    semantics) and gene deletion are done by cobra; this adds RAVEN's choice of
    what to do with reactions whose GPR becomes empty (no enzyme left):

    * ``"remove"`` — delete them (cobra's default; RAVEN ``removeBlockedRxns=true``).
    * ``"constrain"`` — keep them but set bounds to ``(0, 0)`` (RAVEN default).
    * ``"keep"`` — leave them with an empty GPR and unchanged bounds.

    Returns
    -------
    list of str
        IDs of the reactions that became flux-blocked (had a GPR, now empty).
    """
    if blocked_reactions not in ("remove", "constrain", "keep"):
        raise ValueError(
            f"blocked_reactions must be 'remove', 'constrain', or 'keep', "
            f"got {blocked_reactions!r}"
        )

    # Resolve to gene IDs that are actually in the model (RAVEN filters likewise).
    requested = [g.id if isinstance(g, Gene) else g for g in _as_list(genes)]
    present = [gid for gid in requested if gid in model.genes]
    if not present:
        return []

    # Reactions touched by these genes that currently have a GPR.
    affected = set()
    for gid in present:
        affected.update(r.id for r in model.genes.get_by_id(gid).reactions)
    had_gpr = {rid for rid in affected if model.reactions.get_by_id(rid).gene_reaction_rule}

    # cobra rewrites GPRs (AST) and removes the gene objects; we manage reactions.
    _cobra_remove_genes(model, present, remove_reactions=False)

    blocked = [
        rid for rid in had_gpr if not model.reactions.get_by_id(rid).gene_reaction_rule
    ]

    if blocked_reactions == "remove":
        model.remove_reactions(blocked, remove_orphans=remove_orphan_metabolites)
    elif blocked_reactions == "constrain":
        for rid in blocked:
            model.reactions.get_by_id(rid).bounds = (0, 0)

    return sorted(blocked)
