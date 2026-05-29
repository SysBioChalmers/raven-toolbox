"""Model comparison utilities.

Two flavours:

* :func:`compare_models` — N-model presence-matrix overview (RAVEN's
  ``compareMultipleModels`` analogue). "How do these models relate?"
* :func:`diff_models` — strict two-model semantic-equality diff for CI
  gates. "Are these two models the same?"
"""
from raven_python.comparison.compare import ModelComparison, compare_models
from raven_python.comparison.diff import (
    DEFAULT_ANNOTATION_KEYS,
    DiffReport,
    diff_models,
)

__all__ = [
    "DEFAULT_ANNOTATION_KEYS",
    "DiffReport",
    "ModelComparison",
    "compare_models",
    "diff_models",
]
