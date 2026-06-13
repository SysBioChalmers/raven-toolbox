"""Configuration types for the biomass module."""
from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Literal

#: Strategies for computing a component's mass contribution from the
#: stoichiometry of its pseudoreaction. See :func:`sum_biomass`.
#:
#: ``"mw"``                 multiply each substrate coefficient by the
#:                          molecular weight derived from its
#:                          ``formula`` (g/mol), divide by 1000 to get
#:                          g/gDW. Suits carbohydrate / ion / cofactor.
#: ``"mw_minus_2h"``        as ``"mw"`` but with each MW reduced by
#:                          2.016 g/mol (two protons released per
#:                          charged tRNA on protein-pseudoreaction
#:                          substrates).
#: ``"mw_minus_water"``     as ``"mw"`` but with each MW reduced by
#:                          18.015 g/mol (one water released per
#:                          polymerisation step; RNA / DNA).
#: ``"grams"``              substrate coefficients are already in
#:                          g/gDW (e.g. the lipid-backbone
#:                          pseudoreaction). No MW lookup needed.
MassStrategy = Literal["mw", "mw_minus_2h", "mw_minus_water", "grams"]


@dataclass(frozen=True)
class BiomassComponent:
    """One component of the biomass equation (protein, carbohydrate, …).

    Attributes
    ----------
    name
        Canonical short name — also the ``model.metabolites`` name of
        the metabolite that the pseudoreaction *produces*. Used by
        :func:`rescale_pseudoreaction` to identify the product side.
    pseudoreaction_name
        ``model.reactions[*].name`` of the pseudoreaction.
        cobrapy reactions are looked up by ``name`` (not ``id``)
        because that is how the yeast-GEM convention names them.
    mass_strategy
        How to convert the pseudoreaction's substrates into a mass
        fraction. See :data:`MassStrategy` for the choices and their
        offsets.
    """

    name: str
    pseudoreaction_name: str
    mass_strategy: MassStrategy = "mw"


@dataclass(frozen=True)
class BiomassConfig:
    """Container for the per-organism biomass layout.

    Attributes
    ----------
    biomass_rxn
        ``model.reactions[*].id`` of the top-level biomass
        pseudoreaction (the one whose flux is the growth rate).
    proton_met
        ``model.metabolites[*].id`` of the cytosolic H+ metabolite —
        used by :func:`rescale_pseudoreaction` to keep each
        pseudoreaction charge-balanced after rescaling.
    components
        Ordered tuple of :class:`BiomassComponent` describing every
        non-top-level pseudoreaction that contributes to mass.
    """

    biomass_rxn: str
    proton_met: str
    components: tuple[BiomassComponent, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        # Freeze the components into a tuple even if a list was passed.
        if not isinstance(self.components, tuple):
            object.__setattr__(self, "components", tuple(self.components))

    @classmethod
    def from_components(
        cls,
        biomass_rxn: str,
        proton_met: str,
        components: Iterable[BiomassComponent],
    ) -> BiomassConfig:
        return cls(biomass_rxn=biomass_rxn, proton_met=proton_met,
                   components=tuple(components))

    def get(self, component_name: str) -> BiomassComponent:
        """Look up a component by name (raises ``KeyError`` if missing)."""
        for c in self.components:
            if c.name == component_name:
                return c
        raise KeyError(f"BiomassConfig has no component named {component_name!r}")
