"""Remove metabolites or genes from a model.

Partial port of RAVEN ``removeMets.m`` / ``removeGenes.m``.

``removeReactions`` is **not** ported: once orphan-metabolite and orphan-gene
cleanup are kept coupled (as decided — they are not separated), it is exactly
``cobra.Model.remove_reactions(reactions, remove_orphans=...)`` with nothing to
add. Use cobra's method directly.

The two functions kept here delegate the core to cobra and add only the
cobra-absent behaviour:

* ``remove_metabolites`` — cobra matches metabolites by ID; RAVEN's ``isNames``
  deletes a metabolite in **every compartment at once** by name. That name
  resolution is the *sole* reason this wrapper exists (see the note on it).
* ``remove_genes`` — cobra's ``cobra.manipulation.remove_genes`` already rewrites
  GPRs through the boolean AST (removing one gene of ``A and B`` empties the
  rule, of ``A or B`` keeps the other) — exactly RAVEN's intent, without its
  ``eval``. The gap is RAVEN's default of **constraining** flux-blocked reactions
  to zero instead of deleting them; exposed as ``blocked_reactions``.
"""
from __future__ import annotations

from collections.abc import Iterable

import cobra
from cobra import Gene, Metabolite
from cobra.manipulation import remove_genes as _cobra_remove_genes


def _as_list(obj) -> list:
    if isinstance(obj, (str, Metabolite, Gene)):
        return [obj]
    return list(obj)


def remove_metabolites(
    model: cobra.Model,
    metabolites: str | Metabolite | Iterable,
    *,
    by_name: bool = False,
    destructive: bool = False,
) -> None:
    """Remove metabolites, optionally matching by name across all compartments.

    Partial port of RAVEN ``removeMets.m``.

    Parameters
    ----------
    by_name
        If True, ``metabolites`` are metabolite *names*; every metabolite with a
        matching name is removed, regardless of compartment (RAVEN ``isNames``).
        If False, they are IDs/objects, resolved via cobra.
    destructive
        Passed to cobra: if True, also remove every reaction the metabolite
        participates in.

    Note
    ----
    With ``by_name=False`` this is just ``model.remove_metabolites`` — so the
    ``by_name`` cross-compartment deletion is the only thing this adds over
    cobra. It is likely a rare need (one usually knows the compartment-specific
    IDs); if it proves unused, drop this wrapper and call cobra directly.
    """
    if by_name:
        wanted = set(_as_list(metabolites))
        targets = [m for m in model.metabolites if m.name in wanted]
    else:
        targets = model.metabolites.get_by_any(_as_list(metabolites))
    if targets:
        model.remove_metabolites(targets, destructive=destructive)


def remove_genes(
    model: cobra.Model,
    genes: str | Gene | Iterable,
    *,
    blocked_reactions: str = "remove",
    remove_orphans: bool = False,
) -> list[str]:
    """Remove genes and handle reactions left unable to carry flux.

    Port of RAVEN ``removeGenes.m``. GPR rewriting (with correct AND/OR
    semantics) and gene deletion are done by cobra; this adds RAVEN's choice of
    what to do with reactions whose GPR becomes empty (no enzyme left):

    * ``"remove"`` — delete them (cobra's default; RAVEN ``removeBlockedRxns=true``).
    * ``"constrain"`` — keep them but set bounds to ``(0, 0)`` (RAVEN default).
    * ``"keep"`` — leave them with an empty GPR and unchanged bounds.

    ``remove_orphans`` (only meaningful with ``blocked_reactions="remove"``)
    passes through to cobra: drop metabolites *and* genes orphaned by the removal.

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
        model.remove_reactions(blocked, remove_orphans=remove_orphans)
    elif blocked_reactions == "constrain":
        for rid in blocked:
            model.reactions.get_by_id(rid).bounds = (0, 0)

    return sorted(blocked)
