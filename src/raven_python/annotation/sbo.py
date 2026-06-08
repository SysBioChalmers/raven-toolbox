"""SBO term assignment.

Port of the generic core of yeast-GEM's ``code/missingFields/addSBOterms.m``.
Defaults reproduce the yeast-GEM behaviour; overrides parameterise the
pseudo-metabolite / pseudo-reaction name conventions, the transport
detector, and a yeast-specific bug-compat flag.
"""
from __future__ import annotations

from collections.abc import Callable, Iterable

import cobra

DEFAULT_BIOMASS_MET_NAMES: frozenset[str] = frozenset(
    {"biomass", "DNA", "RNA", "protein", "carbohydrate",
     "lipid", "cofactor", "ion"}
)
DEFAULT_BIOMASS_RXN_NAME: str = "biomass pseudoreaction"
DEFAULT_NGAM_RXN_NAME: str = "non-growth associated maintenance reaction"


def add_sbo_terms(
    model: cobra.Model,
    *,
    biomass_met_names: Iterable[str] = DEFAULT_BIOMASS_MET_NAMES,
    biomass_met_suffixes: Iterable[str] = (" backbone", " chain"),
    biomass_rxn_name: str = DEFAULT_BIOMASS_RXN_NAME,
    ngam_rxn_name: str = DEFAULT_NGAM_RXN_NAME,
    pseudoreaction_name_substrings: Iterable[str] = ("pseudoreaction", "SLIME rxn"),
    transport_detector: Callable[[cobra.Model], set[str]] | None = None,
    only_last_reaction_for_pseudo: bool = False,
) -> cobra.Model:
    """Assign SBO terms to metabolites and reactions in-place.

    Metabolite SBO assignment
        - SBO:0000649 (Biomass) for metabolites whose ``name`` is in
          ``biomass_met_names`` or ends with any of
          ``biomass_met_suffixes`` (default: ``" backbone"`` / ``" chain"``).
        - SBO:0000247 (Simple chemical) otherwise.

    Reaction SBO assignment
        - SBO:0000176 (Metabolic reaction) default.
        - Single-reactant reactions (exchange / sink / demand):
            * extracellular → SBO:0000627 (exchange)
            * coefficient < 0 → SBO:0000632 (sink)
            * else → SBO:0000628 (demand)
        - Transport reactions → SBO:0000655 (default detector: same
          metabolite name appearing in ≥ 2 compartments).
        - The reaction whose ``name`` matches ``biomass_rxn_name`` →
          SBO:0000629.
        - The reaction whose ``name`` matches ``ngam_rxn_name`` →
          SBO:0000630.
        - Other reactions whose ``name`` contains any
          ``pseudoreaction_name_substrings`` → SBO:0000395.

    "fill" semantic: SBO is written to ``annotation['sbo']`` only when
    that key is missing or empty — mirrors RAVEN's
    ``editMiriam(..., 'fill')`` mode.

    Parameters
    ----------
    transport_detector
        Optional callable that takes the model and returns the set of
        reaction ids classified as transport. When ``None``, the default
        same-met-name-in-two-compartments heuristic is used. Pass a
        custom detector to use a different convention (e.g. RAVEN's
        ``getTransportRxns`` via the MATLAB engine).
    only_last_reaction_for_pseudo
        **Bug-compatibility flag**, off by default. The legacy yeast-GEM
        MATLAB ``addSBOterms.m`` contains a typo (``for i = numel(...)``
        rather than ``for i = 1:numel(...)``) that causes the pseudo-
        reaction SBO assignments to run only on the very last reaction
        in the model. When ``True``, this flag replicates that bug so
        callers can migrate without altering the committed model. Leave
        ``False`` for new uses.
    """
    biomass_met_names = frozenset(biomass_met_names)
    biomass_met_suffixes = tuple(biomass_met_suffixes)
    pseudoreaction_name_substrings = tuple(pseudoreaction_name_substrings)
    if transport_detector is None:
        transport_detector = _default_transport_detector

    _assign_metabolite_sbo(model, biomass_met_names, biomass_met_suffixes)
    _assign_reaction_sbo(
        model,
        biomass_rxn_name=biomass_rxn_name,
        ngam_rxn_name=ngam_rxn_name,
        pseudo_substrings=pseudoreaction_name_substrings,
        transport_detector=transport_detector,
        only_last_reaction_for_pseudo=only_last_reaction_for_pseudo,
    )
    return model


# --- internals ---------------------------------------------------------

def _assign_metabolite_sbo(
    model: cobra.Model,
    biomass_met_names: frozenset[str],
    biomass_met_suffixes: tuple[str, ...],
) -> None:
    for met in model.metabolites:
        if met.name in biomass_met_names or met.name.endswith(biomass_met_suffixes):
            sbo = "SBO:0000649"
        else:
            sbo = "SBO:0000247"
        _fill_sbo(met, sbo)


def _assign_reaction_sbo(
    model: cobra.Model,
    *,
    biomass_rxn_name: str,
    ngam_rxn_name: str,
    pseudo_substrings: tuple[str, ...],
    transport_detector: Callable[[cobra.Model], set[str]],
    only_last_reaction_for_pseudo: bool,
) -> None:
    """Compute the SBO for each reaction via successive overrides, then
    fill into the model in a single pass — mirrors the MATLAB
    ``rxnSBO``-array + ``editMiriam(..., 'fill')`` two-step.
    """
    transport_set = transport_detector(model)
    rxns = list(model.reactions)

    sbo_for_rxn: dict[str, str] = {}
    for rxn in rxns:
        sbo = "SBO:0000176"
        if len(rxn.metabolites) == 1:
            (met,) = rxn.metabolites
            coef = rxn.metabolites[met]
            if met.compartment == "e" or model.compartments.get(
                met.compartment
            ) == "extracellular":
                sbo = "SBO:0000627"
            elif coef < 0:
                sbo = "SBO:0000632"
            else:
                sbo = "SBO:0000628"
        sbo_for_rxn[rxn.id] = sbo

    # Transport override
    for rxn_id in transport_set:
        sbo_for_rxn[rxn_id] = "SBO:0000655"

    # Pseudoreaction override. The legacy bug-compat flag scopes this
    # loop to the last reaction only; the default walks the whole list.
    pseudo_targets = [rxns[-1]] if only_last_reaction_for_pseudo else rxns
    for rxn in pseudo_targets:
        if rxn.name == biomass_rxn_name:
            sbo_for_rxn[rxn.id] = "SBO:0000629"
        elif rxn.name == ngam_rxn_name:
            sbo_for_rxn[rxn.id] = "SBO:0000630"
        elif any(s in rxn.name for s in pseudo_substrings):
            sbo_for_rxn[rxn.id] = "SBO:0000395"

    for rxn in rxns:
        _fill_sbo(rxn, sbo_for_rxn[rxn.id])


def _fill_sbo(entity, sbo: str) -> None:
    """Set ``annotation['sbo']`` only if it is missing or empty."""
    if not entity.annotation.get("sbo"):
        entity.annotation["sbo"] = sbo


def _default_transport_detector(model: cobra.Model) -> set[str]:
    """Same-met-name-in-two-compartments heuristic for transport detection.

    Cheap analogue of RAVEN's ``getTransportRxns``. Pass a custom
    callable to :func:`add_sbo_terms` to override.
    """
    out: set[str] = set()
    for rxn in model.reactions:
        by_name: dict[str, set[str | None]] = {}
        for met in rxn.metabolites:
            by_name.setdefault(met.name, set()).add(met.compartment)
        if any(len(comps) >= 2 for comps in by_name.values()):
            out.add(rxn.id)
    return out
