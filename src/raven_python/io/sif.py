"""Export a model to Cytoscape SIF (Simple Interaction Format).

Three graph types are supported:

* ``"rc"`` reaction–compound: each reaction linked to its metabolites;
* ``"rr"`` reaction–reaction: reactions linked when they share a metabolite;
* ``"cc"`` compound–compound: each substrate linked to the products of the
  reactions it feeds (computed on an irreversible copy, as RAVEN does, to avoid
  spurious double links from reversible reactions).

A SIF line is ``source <tab> graph_type <tab> target1 <tab> target2 ...``.
"""
from __future__ import annotations

import warnings
from collections import Counter
from collections.abc import Mapping
from pathlib import Path

import cobra

from raven_python.manipulation.irreversible import convert_to_irreversible

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
    model: cobra.Model,
    path: str | Path,
    graph_type: str = "rc",
    *,
    reaction_labels: Mapping[str, str] | None = None,
    metabolite_labels: Mapping[str, str] | None = None,
) -> None:
    """Write ``model`` to a Cytoscape SIF file.

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

    # Warn when the label maps collapse multiple distinct ids onto the same
    # label: target-side dedup runs on labels, so the collision silently merges
    # two nodes into one edge. Only check the ids actually mapped (cobra default
    # labels are ids, which can't collide).
    for kind, lmap in (("reaction", rlabels), ("metabolite", mlabels)):
        duplicates = [lab for lab, n in Counter(lmap.values()).items() if n > 1]
        if duplicates:
            warnings.warn(
                f"{kind}_labels maps multiple ids to the same label(s) "
                f"({duplicates[:5]}{'…' if len(duplicates) > 5 else ''}); "
                "SIF nodes are keyed by label, so those nodes will collapse.",
                stacklevel=2,
            )

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
