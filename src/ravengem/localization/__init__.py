"""Sub-cellular localisation — predictor-agnostic, partial-update friendly.

:func:`predict_localization` is the MILP entry point;
:func:`load_wolfpsort` / :func:`load_deeploc` parse predictor outputs into the
``gene × compartment`` :class:`LocalizationScores` DataFrame the algorithm consumes.
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
