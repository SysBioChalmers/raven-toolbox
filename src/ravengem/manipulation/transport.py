"""Add transport reactions between compartments.

Port of RAVEN ``addTransport.m``. cobra has **no** transport-reaction primitive,
so this is genuinely cobra-absent: for each metabolite it matches the species by
*name* across compartments (the source in ``from_compartment`` and its same-named
twin in each target compartment), optionally creating the target metabolite, and
builds a ``-1 from / +1 to`` reaction with a sequential ``tr_0001`` ID.
"""
from __future__ import annotations

import re
from typing import Iterable, Union

import cobra
from cobra import Metabolite, Reaction

from ravengem.manipulation.add import _new_met_id


def _transport_id_factory(model: "cobra.Model", prefix: str):
    pattern = re.compile(rf"^{re.escape(prefix)}(\d+)$")
    used = [int(m.group(1)) for r in model.reactions if (m := pattern.match(r.id))]
    counter = max(used) + 1 if used else 1

    def next_id() -> str:
        nonlocal counter
        while f"{prefix}{counter:04d}" in model.reactions:
            counter += 1
        rid = f"{prefix}{counter:04d}"
        counter += 1
        return rid

    return next_id


def add_transport_reactions(
    model: "cobra.Model",
    from_compartment: str,
    to_compartments: Union[str, Iterable[str]],
    metabolite_names: Union[str, Iterable[str], None] = None,
    *,
    reversible: bool = True,
    only_to_existing: bool = True,
    id_prefix: str = "tr_",
) -> list[Reaction]:
    """Add transport reactions from one compartment to one or more others.

    Port of RAVEN ``addTransport.m``.

    Parameters
    ----------
    from_compartment
        Source compartment id.
    to_compartments
        Target compartment id(s).
    metabolite_names
        Names of metabolites to transport. Default: every metabolite in
        ``from_compartment``.
    reversible
        If True (default), bounds span the cobra configuration default
        (reversible); otherwise lower bound 0.
    only_to_existing
        If True (default), only transport a metabolite into a target
        compartment where a same-named metabolite already exists. If False,
        create the missing target metabolite (copying name/formula/charge/
        annotation from the source) before adding the transport.
    id_prefix
        Prefix for the sequential reaction IDs (``tr_0001``, ...).

    Returns
    -------
    list of cobra.Reaction
        The transport reactions added, in creation order.
    """
    # cobra's `model.compartments` only lists compartments that have metabolites;
    # include registered-but-empty ones so transport can target an empty compartment.
    known = set(model.compartments) | set(model._compartments)
    if from_compartment not in known:
        raise ValueError(f"Compartment {from_compartment!r} is not in the model.")
    if isinstance(to_compartments, str):
        to_compartments = [to_compartments]
    else:
        to_compartments = list(to_compartments)
    for comp in to_compartments:
        if comp not in known:
            raise ValueError(f"Compartment {comp!r} is not in the model.")

    source = {m.name: m for m in model.metabolites if m.compartment == from_compartment}
    if metabolite_names is None:
        names = list(source)
    else:
        names = [metabolite_names] if isinstance(metabolite_names, str) else list(metabolite_names)
        missing = [n for n in names if n not in source]
        if missing:
            raise ValueError(
                f"Metabolites not found in compartment {from_compartment!r}: {missing}"
            )

    cfg = cobra.Configuration()
    bounds = (cfg.lower_bound, cfg.upper_bound) if reversible else (0.0, cfg.upper_bound)
    from_name = model.compartments.get(from_compartment) or from_compartment
    next_id = _transport_id_factory(model, id_prefix)

    added: list[Reaction] = []
    for to_comp in to_compartments:
        to_name = model.compartments.get(to_comp) or to_comp
        targets = {m.name: m for m in model.metabolites if m.compartment == to_comp}
        for name in names:
            src = source[name]
            dst = targets.get(name)
            if dst is None:
                if only_to_existing:
                    continue
                dst = Metabolite(
                    _new_met_id(model, "m"),
                    name=name,
                    compartment=to_comp,
                    formula=src.formula,
                    charge=src.charge,
                )
                dst.annotation = dict(src.annotation)
                model.add_metabolites([dst])
                targets[name] = dst

            rxn = Reaction(next_id())
            rxn.name = f"{name} transport, {from_name}-{to_name}"
            rxn.bounds = bounds
            model.add_reactions([rxn])
            rxn.add_metabolites({src: -1, dst: 1})
            added.append(rxn)

    return added
