# Gapfilling parameter benchmarks

Functions: `raven_toolbox.gapfilling.fill.connect_blocked_reactions`,
`raven_toolbox.gapfilling.fast_lp.fill_gaps_fast_lp`,
`raven_toolbox.gapfilling.kumar_milp.fill_gaps_kumar_milp`

Date: 2026-06-20. No MATLAB equivalent for most parameters
(gapfilling uses pure Python implementations).

---

## `connect_blocked_reactions` — `eps`

**Parameter:** `eps=1.0` (Python default, no MATLAB equivalent)

`eps` serves two roles: (1) it filters which blocked reactions are gap-fillable —
only reactions whose maximum feasible flux (after merging the template) exceeds `eps`
are targeted; (2) it sets the minimum flux the reaction must carry after gap-filling
(`lower_bound = eps`). The comparison is strict (`fva.max > eps`).

**Benchmark (2026-06-21): synthetic model — supply-limited gap**

Model: A import (max 0.5 units), A→B (irreversible, blocked: B dead-end).
Template: EX_B (B export). Max feasible flux through R_AB with template = 0.5.

| `eps` | Gap filled? | Added reactions | FVA_max > eps? |
|---|---|---|---|
| `0.01` | ✓ yes | `['EX_B']` | 0.5 > 0.01 ✓ |
| `0.1` | ✓ yes | `['EX_B']` | 0.5 > 0.1 ✓ |
| `0.5` | ✗ no | `[]` | 0.5 > 0.5 ✗ (strict) |
| `1.0` | ✗ no | `[]` | 0.5 > 1.0 ✗ |

The condition is strictly greater than (`>`), so a reaction at exactly `eps` capacity
is NOT targeted.

**Decision: ✓ keep `eps=1.0`.** For RAVEN-convention models where exchange reactions
are unconstrained (bounds ±1000) and all real metabolic fluxes easily exceed 1.0,
eps=1.0 is a reliable noise filter. Edge case: if nutrient supply exchanges are
tightly constrained to < 1 mmol/gDW/h (experimental flux constraints, enzyme-
constrained models), eps=1.0 will incorrectly skip gap-fillable reactions. In that
case lower eps to 0.01 or 0.001.

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

**Benchmark (2026-06-21): synthetic model — two equi-capacity repair options**

Model: A import (max 10 units), R1: A→B, R2: C→B (irreversible), biomass: C→{}.
Biomass blocked because C has no source. Two repairs:
- **Reversal of R2** (B→C): A→B→C→biomass via R2 reversed. Cost = `w_rev`.
- **Add EX_C** (C import, max 10): C→biomass directly. Cost = `w_db`.

Both repairs give max biomass = 10. min_growth auto-set to 1.0 (10% of FBA max).

### `weights=(w_rev, w_db)` — reversal vs DB-addition cost ratio

| `weights` | Repair chosen | Biomass |
|---|---|---|
| `(1.0, 2.0)` (default) | R2 reversal | 10.0 |
| `(2.0, 1.0)` | EX_C addition | 10.0 |
| `(1.0, 1.0)` | R2 reversal (tie-break) | 10.0 |

**Decision: ✓ keep `weights=(1.0, 2.0)`.** Correctly implements Kumar 2007 preference
ordering: reversals (weaker biochemical assumption) are preferred over adding entirely
new reactions. The 2:1 ratio ensures a clear preference in ambiguous cases.

### `big_m=1000.0` — reversal capacity cap

The big-M controls how much flux a reversed reaction can carry: `rxn.lower_bound`
is set to `-big_m` when the reversal binary is 1. Same synthetic model, EX_A capped
at 10, so correct reversal capacity requires `big_m ≥ 10`.

| `big_m` | Repair chosen | Biomass |
|---|---|---|
| `5.0` (too small) | R2 reversal | 5.0 (capped: `big_m < 10`) |
| `10.0` (exact) | R2 reversal | 10.0 |
| `1000.0` (default) | R2 reversal | 10.0 |

**Decision: ✓ keep `big_m=1000.0`.** Matches RAVEN-convention model bounds (±1000),
ensuring reversed reactions are not artificially capped below their real capacity.
Note: for enzyme-constrained models where reaction bounds exceed ±1000, `big_m`
must be increased to match the model's maximum bound magnitude.
