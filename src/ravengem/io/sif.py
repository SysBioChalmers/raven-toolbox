"""Export a model to Cytoscape SIF (Simple Interaction Format).

Port of RAVEN ``exportModelToSIF.m``. cobra has no SIF / network-graph export, so
this is genuinely cobra-absent. Three graph types are supported:

* ``"rc"`` reaction–compound: each reaction linked to its metabolites;
* ``"rr"`` reaction–reaction: reactions linked when they share a metabolite;
* ``"cc"`` compound–compound: each substrate linked to the products of the
  reactions it feeds (computed on an irreversible copy, as RAVEN does, to avoid
  spurious double links from reversible reactions).

A SIF line is ``source <tab> graph_type <tab> target1 <tab> target2 ...``.
"""
from __future__ import annotations

from pathlib import Path
from typing import Mapping, Union

import cobra

from ravengem.manipulation.irreversible import convert_to_irreversible

_GRAPH_TYPES = ("rc", "rr", "cc")


def _edges(model, graph_type):
    """Yield (source_object, [target_objects]) per the graph type."""
    if graph_type == "rc":
        for rxn in model.reactions:
            yield rxn, list(rxn.metabolites)
    elif graph_type == "rr":
        for rxn in model.reactions:
            neighbours = {r for met in rxn.metabolites for r in met.reactions}
            neighbours.discard(rxn)
            yield rxn, list(neighbours)
    else:  # cc — on an irreversible copy
        irrev = model.copy()
        convert_to_irreversible(irrev)
        for met in irrev.metabolites:
            products: set = set()
            for rxn in met.reactions:
                if rxn.get_coefficient(met) < 0:  # met is a substrate here
                    products.update(m for m, c in rxn.metabolites.items() if c > 0)
            yield met, list(products)


def export_model_to_sif(
    model: "cobra.Model",
    path: Union[str, Path],
    graph_type: str = "rc",
    *,
    reaction_labels: Mapping[str, str] | None = None,
    metabolite_labels: Mapping[str, str] | None = None,
) -> None:
    """Write ``model`` to a Cytoscape SIF file.

    Port of RAVEN ``exportModelToSIF.m``.

    Parameters
    ----------
    graph_type
        ``"rc"`` (reaction–compound, default), ``"rr"`` (reaction–reaction), or
        ``"cc"`` (compound–compound).
    reaction_labels, metabolite_labels
        Optional ``{id: label}`` maps overriding the node labels (default: IDs).
    """
    if graph_type not in _GRAPH_TYPES:
        raise ValueError(f"graph_type must be one of {_GRAPH_TYPES}, got {graph_type!r}")

    rlabels = reaction_labels or {}
    mlabels = metabolite_labels or {}

    def label(obj) -> str:
        if isinstance(obj, cobra.Reaction):
            return rlabels.get(obj.id, obj.id)
        return mlabels.get(obj.id, obj.id)

    with open(path, "w", encoding="utf-8") as handle:
        for source, targets in _edges(model, graph_type):
            src = label(source)
            names = sorted({label(t) for t in targets} - {src})
            if names:
                handle.write(f"{src}\t{graph_type}\t" + "\t".join(names) + "\n")
