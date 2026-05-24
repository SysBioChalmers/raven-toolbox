"""RAVEN template-based MILP gap-filling (port of ``fillGaps``, Phase 4b).

* :func:`fill_gaps` — connectivity gap-fill: add template reactions so blocked draft
  reactions can carry flux (``fillGaps(..., useModelConstraints=false)``). No cobra
  equivalent.
* :func:`gapfill_to_objective` — targeted gap-fill to an objective lower bound
  (``fillGaps(..., useModelConstraints=true)``); name-matching analogue of
  ``cobra.flux_analysis.gapfill``.
"""
from ravengem.gapfilling.fill import GapFillResult, fill_gaps, gapfill_to_objective

__all__ = ["GapFillResult", "fill_gaps", "gapfill_to_objective"]
