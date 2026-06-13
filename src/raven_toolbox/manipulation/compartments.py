"""Compartment manipulation — merge all compartments into one, or copy reactions to a
new compartment (ports of RAVEN's ``mergeCompartments`` and ``copyToComps``).

Both functions are useful **independently of** :func:`raven_toolbox.localization.predict_localization`:
``merge_compartments`` flattens a multi-compartment model for a simplified analysis
(e.g. checking whether the network can in principle make a metabolite, with no
compartment topology in the way); ``copy_to_compartment`` is a building block for
constructing dual-localised pathways. cobra has no equivalents.
"""
from __future__ import annotations

from collections.abc import Iterable

import cobra

# Compartments produced by merge_compartments (RAVEN uses 's' for "system").
_MERGED_COMPARTMENT = "s"


def merge_compartments(
    model: cobra.Model,
    *,
    merged_id: str = _MERGED_COMPARTMENT,
    merged_name: str = "system",
    drop_single_metabolite_reactions: bool = True,
    deduplicate_reactions: bool = True,
) -> tuple[cobra.Model, list[str], list[str]]:
    """Merge every metabolite of ``model`` into one ``merged_id`` compartment.

    Returns ``(model_copy, deleted_single_met_reactions, deduplicated_reactions)``. The
    returned model is a deep copy of the input. Use cases:

    * Check whether the network can produce/consume a metabolite at all (compartment
      topology is often what makes a model look blocked).
    * Simplify a model for visualisation or an analysis that doesn't care about
      compartments.
    * As a pre-step for localisation when the user does want RAVEN's
      "start from scratch" workflow (call :func:`merge_compartments` then
      :func:`raven_toolbox.localization.predict_localization` with the full reaction list).

    Metabolites that already share a base id (e.g. ``glc__D_c`` and ``glc__D_e`` both
    map to ``glc__D``) collapse into one entity in the merged compartment; their
    stoichiometric contributions are summed per reaction. Reactions that end up with
    only one metabolite (e.g. ``A[c] → A[m]`` becomes ``A → A`` = nothing) are deleted
    by default (RAVEN's ``deleteRxnsWithOneMet``). Reactions that become identical
    after merging are deduplicated (one survives).
    """
    out = model.copy()

    # 1. For each metabolite, derive a base id (strip the trailing _<compartment>).
    #    Two mets in different compartments sharing the base id collapse to one.
    new_to_old: dict[str, list[cobra.Metabolite]] = {}
    for m in list(out.metabolites):
        base = _base_id(m)
        new_to_old.setdefault(base, []).append(m)

    # 2. Build the merged metabolites and rewrite reactions.
    canonical: dict[str, cobra.Metabolite] = {}
    for base, mets in new_to_old.items():
        proto = mets[0]
        new_met = cobra.Metabolite(base, name=proto.name, compartment=merged_id,
                                    formula=proto.formula, charge=proto.charge)
        new_met.notes = dict(proto.notes or {})
        canonical[base] = new_met

    # Rewrite all reactions: replace each metabolite with its canonical, summing
    # coefficients where multiple original mets collapse to one.
    rewritten: dict[str, dict[str, float]] = {}
    for r in list(out.reactions):
        new_stoich: dict[cobra.Metabolite, float] = {}
        for m, coeff in list(r.metabolites.items()):
            canon = canonical[_base_id(m)]
            new_stoich[canon] = new_stoich.get(canon, 0.0) + coeff
        # Drop zero net coefficients (substrate + product of the same base met cancel).
        new_stoich = {m: c for m, c in new_stoich.items() if c != 0.0}
        rewritten[r.id] = {m.id: c for m, c in new_stoich.items()}

    # Now build a fresh model with the canonical mets + rewritten reactions; the
    # cobra in-place rewrite would require careful constraint surgery, so a clean
    # rebuild is simpler and less error-prone.
    merged = cobra.Model(out.id or "merged")
    merged.compartments = {merged_id: merged_name}
    merged.add_metabolites(list(canonical.values()))
    deleted_single: list[str] = []
    deduplicated: list[str] = []
    seen_signatures: dict[tuple, str] = {}
    keep_reactions: list[cobra.Reaction] = []
    for r in out.reactions:
        stoich = rewritten[r.id]
        if drop_single_metabolite_reactions and len(stoich) <= 1:
            deleted_single.append(r.id)
            continue
        if not stoich:  # everything cancelled
            deleted_single.append(r.id)
            continue
        sig = (frozenset(stoich.items()), bool(r.lower_bound < 0), bool(r.upper_bound > 0))
        if deduplicate_reactions and sig in seen_signatures:
            deduplicated.append(r.id)
            continue
        seen_signatures[sig] = r.id
        new_r = cobra.Reaction(r.id, name=r.name, lower_bound=r.lower_bound,
                                upper_bound=r.upper_bound)
        new_r.add_metabolites({merged.metabolites.get_by_id(mid): c for mid, c in stoich.items()})
        new_r.gene_reaction_rule = r.gene_reaction_rule
        if r.subsystem:
            new_r.subsystem = r.subsystem
        new_r.notes = dict(r.notes or {})
        keep_reactions.append(new_r)
    merged.add_reactions(keep_reactions)
    return merged, deleted_single, deduplicated


def copy_to_compartment(
    model: cobra.Model,
    reactions: Iterable[str],
    target_compartment: str,
    *,
    target_compartment_name: str | None = None,
    delete_original: bool = False,
    id_suffix: str | None = None,
) -> tuple[cobra.Model, list[str], list[str]]:
    """Copy a set of reactions into ``target_compartment``. RAVEN's ``copyToComps``.

    Returns ``(model_copy, new_reaction_ids, new_metabolite_ids)``. Use cases:

    * Build a dual-localised pathway (e.g. duplicate glycolysis into a peroxisome).
    * Mirror a curated subsystem into an additional compartment as a draft to refine.
    * Set up the input for a flux comparison between alternate compartmentalisations.

    Each copied reaction is given the id ``"<orig_id>_<id_suffix>"`` (default
    ``id_suffix=target_compartment``); each metabolite it touches is mapped to (or
    created in) ``target_compartment`` with the same suffix convention. ``delete_original=True``
    moves the reactions instead of copying.
    """
    out = model.copy()
    suffix = id_suffix if id_suffix is not None else target_compartment
    if target_compartment not in out.compartments:
        out.compartments = {**out.compartments,
                             target_compartment: target_compartment_name or target_compartment}

    preexisting_met_ids = {x.id for x in out.metabolites}
    new_rxn_ids: list[str] = []
    for rid in list(reactions):
        if rid not in out.reactions:
            raise ValueError(f"reaction {rid!r} not in model")
        src = out.reactions.get_by_id(rid)
        new_id = f"{rid}_{suffix}"
        if new_id in out.reactions:
            continue  # already copied; idempotent
        new_stoich: dict[cobra.Metabolite, float] = {}
        for m, coeff in src.metabolites.items():
            target_met = _met_in_compartment(out, m, target_compartment, suffix=suffix)
            new_stoich[target_met] = coeff
        new_r = cobra.Reaction(new_id, name=src.name,
                                lower_bound=src.lower_bound, upper_bound=src.upper_bound)
        new_r.add_metabolites(new_stoich)
        new_r.gene_reaction_rule = src.gene_reaction_rule
        if src.subsystem:
            new_r.subsystem = src.subsystem
        new_r.notes = dict(src.notes or {})
        out.add_reactions([new_r])
        new_rxn_ids.append(new_id)
        if delete_original:
            out.remove_reactions([src.id], remove_orphans=False)

    new_met_ids = [m.id for m in out.metabolites if m.id not in preexisting_met_ids]
    return out, new_rxn_ids, new_met_ids


# ----------------------------------------------------------------- helpers

def _base_id(m: cobra.Metabolite) -> str:
    """Strip the trailing ``_<compartment>`` suffix from a metabolite id (if present)."""
    if m.compartment and m.id.endswith(f"_{m.compartment}"):
        return m.id[: -(len(m.compartment) + 1)]
    return m.id


def _met_in_compartment(model: cobra.Model, source: cobra.Metabolite,
                        compartment: str, *, suffix: str | None = None) -> cobra.Metabolite:
    """Return (creating if needed) the copy of ``source`` in ``compartment``.

    The new metabolite id is ``"<base>_<suffix>"`` (default ``suffix=compartment``).
    Already-existing copies are reused.
    """
    if source.compartment == compartment:
        return source
    base = _base_id(source)
    new_id = f"{base}_{suffix if suffix is not None else compartment}"
    if new_id in model.metabolites:
        return model.metabolites.get_by_id(new_id)
    new_met = cobra.Metabolite(new_id, name=source.name, compartment=compartment,
                                formula=source.formula, charge=source.charge)
    new_met.notes = dict(source.notes or {})
    model.add_metabolites([new_met])
    return new_met
