"""Gap-filling against template models.

Three complementary strategies are provided:

- :func:`connect_blocked_reactions` (connectivity MILP) — minimum-penalty template
  reactions so blocked draft reactions can carry flux.
- :func:`fill_gaps_fast_lp` (LP, fast) — L1-norm LP per blocked reaction; no MILP solver
  required; scales better to large models (fastGapFill / swiftGapFill).
- :func:`analyse_topology` (no solver) — BFS metabolite-producibility pre-screening;
  identifies unreachable metabolites and prunes the candidate reaction pool.
- :func:`fill_gaps_kumar_milp` (objective MILP) — growth-floor MILP with directionality
  reversal repair (Kumar et al. 2007).

For the objective-feasibility flavour without directionality repair, also see
``cobra.flux_analysis.gapfill``.
"""
from raven_toolbox.gapfilling.fast_lp import FastLPResult, fill_gaps_fast_lp
from raven_toolbox.gapfilling.fill import GapFillResult, connect_blocked_reactions
from raven_toolbox.gapfilling.kumar_milp import KumarGapFillResult, fill_gaps_kumar_milp
from raven_toolbox.gapfilling.topological import TopologicalAnalysisResult, analyse_topology

__all__ = [
    "GapFillResult",
    "connect_blocked_reactions",
    "FastLPResult",
    "fill_gaps_fast_lp",
    "TopologicalAnalysisResult",
    "analyse_topology",
    "KumarGapFillResult",
    "fill_gaps_kumar_milp",
]
