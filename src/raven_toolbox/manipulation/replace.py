"""Merge a metabolite into another already in the model — port of RAVEN's
``replaceMets``.

Useful when a model has ended up with two records for the same real-world
compound under different names or ids (RAVEN's own example: ``"oxygen"`` and
``"o2"``). Reassigns stoichiometry from the old metabolite to the
replacement, copies the replacement's identity/annotation onto it, and
cleans up any reactions that became duplicates as a result.
"""
from __future__ import annotations

import cobra

from raven_toolbox.manipulation.simplify import remove_duplicate_reactions

__all__ = ["replace_metabolite"]

_NOTES_FIELDS = ("inchis", "deltaG", "metFrom")


def _copy_identity(dst: cobra.Metabolite, src: cobra.Metabolite) -> None:
    """Copy name/formula/charge/annotation/notes from src onto dst.

    Matches replaceMets.m's unconditional copy of metNames/metFormulas/
    metMiriams/metCharges/metDeltaG/inchis/metSmiles from the replacement
    onto the to-be-replaced metabolite (annotation covers MIRIAM xrefs and
    smiles; notes covers RAVEN's own inchis/deltaG/metFrom fields — see
    raven_toolbox.io.yaml's field mapping).
    """
    dst.name = src.name
    dst.formula = src.formula
    dst.charge = src.charge
    dst.annotation = dict(src.annotation)
    for key in _NOTES_FIELDS:
        if key in src.notes:
            dst.notes[key] = src.notes[key]
        else:
            dst.notes.pop(key, None)


def _print_touching(reactions, verbose: bool) -> None:
    if not verbose:
        return
    print("\n\nThe following reactions contain the to-be-replaced metabolite as reactant:")
    print("\n".join(sorted(r.id for r in reactions)))


def replace_metabolite(
    model: cobra.Model,
    metabolite: str,
    replacement: str,
    *,
    identifiers: bool = False,
    verbose: bool = False,
) -> list[str]:
    """Replace ``metabolite`` with ``replacement``, an existing model metabolite.

    Parameters
    ----------
    model:
        Model to modify in place.
    metabolite:
        Name (or id, with ``identifiers=True``) of the metabolite to replace.
    replacement:
        Name (or id) of the metabolite to replace it with. Its identity
        (name, formula, charge, annotation, and RAVEN's inchi/deltaG/metFrom
        notes) is copied onto every replaced metabolite.
    identifiers:
        If True, ``metabolite``/``replacement`` are metabolite ids rather
        than names (default False). Ids are unique, so this is always a
        direct 1:1 merge, and the survivor (``replacement``) keeps its own
        existing identity unchanged — ``replaceMets.m`` copies the
        replacement's identity onto the *old* metabolite's row before
        deleting that row, which has no observable effect (it never
        survives to be seen), so this doesn't replicate that dead step.
        Names are not unique: every metabolite named ``metabolite`` is
        renamed to ``replacement``'s name, and then *any* metabolite in the
        model sharing a (name, compartment) with a renamed one is merged
        into it too — not just the ones being explicitly replaced, matching
        ``replaceMets.m`` exactly (its own docstring example: a model with
        both "oxygen" and "o2" as names). Here the identity copy *is*
        observable (a renamed metabolite that turns out not to collide with
        anything survives with the replacement's copied identity), so it is
        applied.
    verbose:
        Print the ids of reactions touching the to-be-replaced metabolite(s).

    Returns
    -------
    list[str]
        Ids of duplicate reactions removed as a result of the merge (only
        reactions touching the replacement are checked — pre-existing
        duplicates elsewhere in the model are left alone).

    Raises
    ------
    ValueError
        If ``metabolite`` or ``replacement`` cannot be found.
    """
    if identifiers:
        if replacement not in model.metabolites:
            raise ValueError(f"The replacement metabolite {replacement!r} cannot be found in the model.")
        if metabolite not in model.metabolites:
            raise ValueError(f"The to-be-replaced metabolite {metabolite!r} cannot be found in the model.")
        rep = model.metabolites.get_by_id(replacement)
        old = model.metabolites.get_by_id(metabolite)

        _print_touching(old.reactions, verbose)

        # Sums into any pre-existing replacement coefficient in the same
        # reaction, rather than overwriting it -- replaceMets.m does the
        # latter (model.S(repIdx,rxnsWithMet) = originalStoch, a plain
        # assignment). Both agree unless a reaction already references both
        # metabolites, an edge case unlikely to occur in a real model (the
        # same compound recorded twice in one reaction is itself a modelling
        # smell); summing is the more correct choice if it ever does.
        for rxn in list(old.reactions):
            coef = rxn.get_coefficient(old)
            rxn.add_metabolites({old: -coef, rep: coef})
        model.remove_metabolites([old])

        return remove_duplicate_reactions(model, restrict_to_metabolites=[rep])

    reps = [m for m in model.metabolites if m.name == replacement]
    if not reps:
        raise ValueError(f"The replacement metabolite {replacement!r} cannot be found in the model.")
    olds = [m for m in model.metabolites if m.name == metabolite]
    if not olds:
        raise ValueError(f"The to-be-replaced metabolite {metabolite!r} cannot be found in the model.")

    _print_touching({r for m in olds for r in m.reactions}, verbose)

    rep = reps[0]
    for m in olds:
        _copy_identity(m, rep)

    # Re-scan the whole model: any metabolites now sharing a (name, compartment)
    # -- not just the ones just renamed -- are duplicates and get merged.
    groups: dict[tuple[str, str], list[cobra.Metabolite]] = {}
    for m in model.metabolites:
        groups.setdefault((m.name, m.compartment), []).append(m)

    to_delete: list[cobra.Metabolite] = []
    survivors: list[cobra.Metabolite] = []
    seen_reps: set[tuple[str, str]] = set()
    for r in reps:
        key = (r.name, r.compartment)
        if key in seen_reps:
            continue
        seen_reps.add(key)
        group = groups.get(key, [r])
        survivor, *duplicates = group
        # Track the metabolite that actually ends up holding this group's
        # stoichiometry (not necessarily r itself -- the first-encountered
        # member of the group survives, same as remove_duplicate_reactions'
        # own survivor rule, and r can be the one that gets merged away).
        # Reactions touching it are what may have become duplicates.
        survivors.append(survivor)
        for dup in duplicates:
            for rxn in list(dup.reactions):
                coef = rxn.get_coefficient(dup)
                rxn.add_metabolites({dup: -coef, survivor: coef})
        to_delete += duplicates

    if to_delete:
        model.remove_metabolites(to_delete)

    return remove_duplicate_reactions(model, restrict_to_metabolites=survivors)
