"""Connectivity gap-filling against template models.

:func:`connect_blocked_reactions` adds the fewest (lowest-penalty) template reactions so
reactions blocked in a draft can carry flux. For the other gap-fill flavour (fill until
the objective is feasible) use ``cobra.flux_analysis.gapfill``.
"""
from ravengem.gapfilling.fill import GapFillResult, connect_blocked_reactions

__all__ = ["GapFillResult", "connect_blocked_reactions"]
