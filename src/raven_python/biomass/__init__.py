"""Biomass equation manipulation — growth-associated maintenance, amino-
acid ratios, component scaling, and biomass-fraction reporting.

The yeast-GEM port (see yeast-GEM/code/python/PORTING_PLAN.md) was the
first consumer; the API is parameterised by a :class:`BiomassConfig`
so other GEMs can describe their own component layout.

A typical caller assembles ``BiomassConfig`` once (often from a
project-level YAML) and passes it to every operation:

.. code-block:: python

    cfg = BiomassConfig(
        biomass_rxn="r_4041",
        proton_met="s_0794",
        components=[
            BiomassComponent("protein", "protein pseudoreaction",
                             mass_strategy="mw_minus_2h"),
            BiomassComponent("carbohydrate", "carbohydrate pseudoreaction"),
            BiomassComponent("lipid_backbone", "lipid backbone pseudoreaction",
                             mass_strategy="grams"),
            BiomassComponent("RNA", "RNA pseudoreaction",
                             mass_strategy="mw_minus_water"),
            BiomassComponent("DNA", "DNA pseudoreaction",
                             mass_strategy="mw_minus_water"),
            BiomassComponent("ion", "ion pseudoreaction"),
            BiomassComponent("cofactor", "cofactor pseudoreaction"),
        ],
    )
    fractions = sum_biomass(model, cfg)            # {'protein': 0.46, ...}
    scale_biomass(model, cfg, "protein", 0.50, balance_out="carbohydrate")
    set_gam(model, 80, biomass_rxn=cfg.biomass_rxn,
            cofactor_met_names=("ATP", "ADP", "H2O", "H+", "phosphate"))
"""
from raven_python.biomass.config import BiomassComponent, BiomassConfig
from raven_python.biomass.gam import set_gam
from raven_python.biomass.scale import (
    rescale_pseudoreaction,
    scale_biomass,
    sum_biomass,
)

__all__ = [
    "BiomassComponent",
    "BiomassConfig",
    "rescale_pseudoreaction",
    "scale_biomass",
    "set_gam",
    "sum_biomass",
]
