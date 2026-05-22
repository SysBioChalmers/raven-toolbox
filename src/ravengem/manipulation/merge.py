"""Merge several models into one.

Port of RAVEN ``mergeModels.m``. cobra's ``Model.merge`` is pairwise and matches
everything strictly by ID; this merges **N** models and unifies metabolites by
**name[compartment]** (so the same compound under different IDs in two models
becomes one), while — like RAVEN — adding **all** reactions without de-duplication
(a reaction whose ID already exists is renamed ``id_<sourceid>``). Genes are
unified by ID. Provenance (which source model each object came from) is recorded
in ``notes['origin']``.

The bulk of RAVEN's function is struct field-padding and manual S-matrix
assembly, none of which is needed on ``cobra.Model``.
"""
from __future__ import annotations

import copy
from typing import Iterable

import cobra
from cobra import Metabolite, Model, Reaction


def _unique_id(existing, base: str, suffix: str) -> str:
    """Return base, or base_suffix (then base_suffix_2, ...) if it collides."""
    if base not in existing:
        return base
    candidate = f"{base}_{suffix}"
    n = 2
    while candidate in existing:
        candidate = f"{base}_{suffix}_{n}"
        n += 1
    return candidate


def merge_models(
    models: Iterable["cobra.Model"],
    *,
    match_by: str = "name",
    track_origin: bool = True,
) -> "cobra.Model":
    """Merge models into a single new model.

    Port of RAVEN ``mergeModels.m``.

    Parameters
    ----------
    models
        The models to merge (two or more). A single model is returned as a copy.
    match_by
        How metabolites are unified across models: ``"name"`` (default) treats
        metabolites with the same *name and compartment* as identical (IDs
        ignored); ``"id"`` matches by metabolite ID.
    track_origin
        If True (default), record the source model's ``id`` in each reaction's,
        metabolite's, and gene's ``notes['origin']``.

    Returns
    -------
    cobra.Model
        A new merged model (``id="MERGED"``). Reactions are **not** de-duplicated
        — matching RAVEN, every reaction from every model is kept, with ID
        collisions renamed ``id_<sourceid>``.
    """
    models = list(models)
    if not models:
        raise ValueError("merge_models requires at least one model.")
    if match_by not in ("name", "id"):
        raise ValueError(f"match_by must be 'name' or 'id', got {match_by!r}")
    if len(models) == 1:
        return models[0].copy()

    merged = Model("MERGED")
    comp_names: dict[str, str] = {}
    met_lookup: dict = {}  # name/comp or id key -> merged Metabolite

    def met_key(met: Metabolite):
        return (met.name, met.compartment) if match_by == "name" else met.id

    def ensure_metabolite(src: Metabolite, origin: str) -> Metabolite:
        key = met_key(src)
        if key in met_lookup:
            return met_lookup[key]
        new_id = _unique_id(merged.metabolites, src.id, origin)
        new_met = Metabolite(
            new_id, name=src.name, compartment=src.compartment,
            formula=src.formula, charge=src.charge,
        )
        new_met.annotation = copy.deepcopy(src.annotation)
        new_met.notes = copy.deepcopy(src.notes)
        if track_origin:
            new_met.notes.setdefault("origin", origin)
        merged.add_metabolites([new_met])
        met_lookup[key] = new_met
        return new_met

    for model in models:
        origin = model.id or "model"
        comp_names.update(model.compartments)
        genes_before = {g.id for g in merged.genes}

        for rxn in model.reactions:
            new_id = _unique_id(merged.reactions, rxn.id, origin)
            new_rxn = Reaction(new_id, name=rxn.name)
            new_rxn.bounds = rxn.bounds
            new_rxn.subsystem = rxn.subsystem
            merged.add_reactions([new_rxn])
            new_rxn.add_metabolites(
                {ensure_metabolite(m, origin): coef for m, coef in rxn.metabolites.items()}
            )
            if rxn.gene_reaction_rule:
                new_rxn.gene_reaction_rule = rxn.gene_reaction_rule
            new_rxn.annotation = copy.deepcopy(rxn.annotation)
            new_rxn.notes = copy.deepcopy(rxn.notes)
            if track_origin:
                new_rxn.notes.setdefault("origin", origin)

        if track_origin:
            for gene in merged.genes:
                if gene.id not in genes_before:
                    gene.notes.setdefault("origin", origin)

    merged._compartments.update(comp_names)
    return merged
