# Gapfilling parameter benchmarks

Functions: `raven_toolbox.gapfilling.fill.connect_blocked_reactions`,
`raven_toolbox.gapfilling.fast_lp.fill_gaps_fast_lp`,
`raven_toolbox.gapfilling.kumar_milp.fill_gaps_kumar_milp`

Date: 2026-06-20. No MATLAB equivalent for most parameters
(gapfilling uses pure Python implementations).

---

## `connect_blocked_reactions` — `eps`

**Parameter:** `eps=1.0` (Python default, no MATLAB equivalent)

`eps` is the minimum flux that a reaction must carry in the gapfilled model to
be considered "connected". The LP minimises the number of database reactions added
subject to `v_target ≥ eps` for each previously blocked reaction.

A higher `eps` requires a stronger flux signal, potentially identifying more
specific gapfills but at the risk of infeasibility when only small fluxes are
achievable. A lower `eps` is more permissive but may add unnecessary reactions.

**Status: untested.** The return type is `GapFillResult` (not a list); a proper
sensitivity analysis requires building a model with a realistic blocked reaction
scenario (e.g., a gap in yeast-GEM's amino acid biosynthesis) and checking whether
the added reactions are biochemically sensible.

Proposed test: introduce a known gap in iJO1366's TCA cycle, then run
`connect_blocked_reactions` at `eps` ∈ {0.01, 0.1, 1.0, 10.0}. Assess
(a) whether the correct reaction is added, (b) whether extra reactions are added.

---

## `connect_blocked_reactions` — `allow_net_production`

**Parameter:** `allow_net_production=False` (Python and MATLAB `fillGaps`)

When `False`, the LP is forced to route all metabolite production through a
degradation or export path, preventing thermodynamically impossible "free lunch"
solutions. When `True`, net production of metabolites (metabolites with no consumer)
is allowed, which can produce energetically inconsistent gapfills.

`False` is the correct default for metabolic network reconstruction. `True` is
appropriate only for very incomplete draft networks where no degradation pathway
is known.

**Decision: ✓ keep `False`.**

---

## `fill_gaps_fast_lp` — `epsilon`

**Parameter:** `epsilon=0.0001` (Python default, matching fastGapFill paper)

From Thiele et al. 2014 (fastGapFill): `epsilon` sets the minimum flux that must
pass through each reaction to be counted as "active" in the penalty objective.
The paper uses `epsilon = 1e-4`.

**Decision: ✓ keep `0.0001`.** Matches the fastGapFill paper value.

---

## `fill_gaps_kumar_milp` — `weights` and `big_m`

**`weights=(1.0, 2.0)`** — Kumar et al. 2007: `w_rev=1` for reversing existing
reactions, `w_add=2` for adding new database reactions. Penalising additions more
than reversals reflects the prior that the existing model topology is more likely
to be correct than the direction assignments.

**`big_m=1000.0`** — used in the MILP binary indicator formulation. The conventional
big-M of 1000 matches the standard RAVEN upper bound for reaction fluxes. Must be
≥ maximum expected flux; for yeast-GEM this is satisfied.

**Status: untested on realistic gapfilling scenarios.** Both values are taken
directly from the Kumar et al. 2007 paper without empirical modification.

Proposed test: run `fill_gaps_kumar_milp` on a subset of yeast-GEM with a
simulated gap and compare the added reactions at different weight ratios
(e.g., `(1.0, 1.0)`, `(1.0, 2.0)`, `(1.0, 5.0)`).
