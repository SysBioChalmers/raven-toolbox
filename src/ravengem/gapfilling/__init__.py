"""RAVEN connectivity gap-filling (port of ``fillGaps``, Phase 4b).

:func:`connect_blocked_reactions` adds the fewest (lowest-penalty) template reactions
so reactions that are *blocked* in a draft can carry flux — RAVEN's
``fillGaps(..., useModelConstraints=false)``, which has no cobra equivalent. The other
RAVEN mode (fill to make the objective feasible) is ``cobra.flux_analysis.gapfill``
(see the PLAN.md cheatsheet), so it is not re-wrapped here.
"""
from ravengem.gapfilling.fill import GapFillResult, connect_blocked_reactions

__all__ = ["GapFillResult", "connect_blocked_reactions"]
