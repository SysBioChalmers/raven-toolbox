# MATLAB RAVEN ↔ raven-toolbox — differences

Record of functionality in MATLAB RAVEN that raven-toolbox deliberately does not have, plus
raven-toolbox functionality still to be back-ported into MATLAB RAVEN.

---

## Pending back-port: functionality-constrained compartment-assignment MILP

raven-toolbox's `localization/` now hosts **two** compartment-assignment algorithms:
`predict_localization` (score-driven MILP) and **`assign_compartments`** (`localization/assign.py`) —
the *functionality-constrained* MILP consolidated from the retired `edkerk/assignCompartments` repo.
Over the score-driven version, `assign_compartments` adds a **biomass/growth floor**, **big-M flux
gating** (a placement carries flux or scores nothing; sound via a tightened integrality tolerance),
optional **gap-fill coupling** (universal-DB candidates added only when biomass feasibility needs them),
and **sound reaction-level multi-localisation** (ε-flux activity coupling forbids dead placements —
design in [multi_localization_design.md](multi_localization_design.md); benchmark in
[assign_compartments on yeast-GEM](../studies/assign_compartments_yeast_gem.md)).

MATLAB RAVEN has **no equivalent**: `core/predictLocalization.m` is a *simulated-annealing* heuristic
(one gene → one compartment, no biomass constraint, no flux gating). The score adapters this needs are
already in RAVEN (`parseScores`, `getUniProtScores`, `defaultCompartmentMap`).

**Port plan.** Add a new `assignCompartments.m` (MILP via RAVEN's `optimizeProb` / `getMILPParams`,
Gurobi/GLPK) that **coexists** with `predictLocalization.m` (mirroring the Python coexistence), reusing
`parseScores` for the `gene × compartment` scores. Verify the `.m` filename does not clash with a COBRA
Toolbox function before committing. Tests under `testing/` mirroring `tests/test_localization_assign*.py`.
Reference implementation: `src/raven_toolbox/localization/assign.py`.

---

## Functionality in MATLAB RAVEN not in raven-toolbox

Principled omissions — present in MATLAB RAVEN, deliberately **not** ported to
raven-toolbox.

* **MetaCyc-based reconstruction** (`external/metacyc/*`, `getMetaCycModelForOrganism`).
  Not ported, and proposed for **removal from MATLAB RAVEN** as well: MetaCyc's one
  representative sequence per enzyme gives intrinsically low gene-calling precision
  (~64 % of reaction assignments wrong at the default cutoff; no cutoff rescues it).
  Full evidence and the MATLAB removal list in
  [IMPROVEMENTS.md § R-MetaCyc](improvements.md).
* **Dynamic FBA.** Not ported — maintained Python packages already cover it
  ([`dfba`](https://pypi.org/project/dfba/), [`reframed`](https://pypi.org/project/reframed/),
  [`mewpy`](https://pypi.org/project/mewpy/)).
* **Metabolomics-based scoring in ftINIT** (the 4d.6 production-bonus block).
  `ftinit(metabolomics=…)` raises `NotImplementedError`. The linear merge eliminates
  degree-2 detected metabolites, so it would need RAVEN's producer-group-mapping +
  `mon`/`vnrbm`/`vnrvm`/`vnim` negative-producer force-flux block — the most intricate
  MILP in ftINIT, for its least-used input.
