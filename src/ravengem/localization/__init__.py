"""Subcellular localisation (Phase 7) — predictor-agnostic, partial-update-friendly.

See [docs/localization_design.md](../../docs/localization_design.md) for the design
rationale (critical review of RAVEN's ``predictLocalization`` + the redesign).
"""
from ravengem.localization.predict import (
    LocalizationProposal,
    LocalizationResult,
    apply_localization,
    predict_localization,
)
from ravengem.localization.scores import (
    LocalizationScores,
    load_deeploc,
    load_wolfpsort,
)

__all__ = [
    "LocalizationProposal",
    "LocalizationResult",
    "LocalizationScores",
    "apply_localization",
    "load_deeploc",
    "load_wolfpsort",
    "predict_localization",
]
