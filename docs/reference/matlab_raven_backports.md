# MATLAB RAVEN ↔ raven-toolbox — differences

Record of functionality in MATLAB RAVEN that raven-toolbox deliberately does not have, plus
raven-toolbox functionality still to be back-ported into MATLAB RAVEN.

---

## Pending back-port: functionality-constrained compartment-assignment MILP

raven-toolbox's `localization/` now hosts **two** compartment-assignment algorithms:
`predict_localization` (score-driven MILP) and **`assign_compartments`** (`localization/certify.py`) —
the *functionality-constrained* method, seeded by the port of the retired `edkerk/assignCompartments`
repo. Over the score-driven version it adds a **biomass/growth floor** enforced by certification:
placement is decided by a flux-free score MILP and the result is confirmed by a real FBA on the
materialised model. It also adds optional **gap-fill coupling** (universal-DB candidates added only when
biomass feasibility needs them) and **sound reaction-level multi-localisation** (a second compartment is
kept only if a loopless FVA on the materialised model shows it carries real flux — design in
[multi_localization_design.md](multi_localization_design.md)).

MATLAB RAVEN has **no equivalent**: `core/predictLocalization.m` is a *simulated-annealing* heuristic
(one gene → one compartment, no biomass constraint, no flux gating). The score adapters this needs are
already in RAVEN (`parseScores`, `getUniProtScores`, `defaultCompartmentMap`).

**Port plan.** Add a new `assignCompartments.m` (MILP via RAVEN's `optimizeProb` / `getMILPParams`,
Gurobi/GLPK) that **coexists** with `predictLocalization.m` (mirroring the Python coexistence), reusing
`parseScores` for the `gene × compartment` scores. Verify the `.m` filename does not clash with a COBRA
Toolbox function before committing. Tests under `testing/` mirroring `tests/test_localization_assign*.py`.
Reference implementation: `src/raven_toolbox/localization/assign.py`.

---

## Pending back-port INTO raven-toolbox: order-insensitive grRule comparison in model diff

MATLAB RAVEN's `comparison/diffModels.m` (the port of `comparison/diff.py::diff_models`, added in
[RAVEN #686](https://github.com/SysBioChalmers/RAVEN/pull/686)) compares grRules as **logic**, not
text: each rule is expanded to disjunctive normal form and the genes within each isozyme, and the
isozymes themselves, are sorted before comparison, so `a and b` equals `b and a` and `(a or b)`
equals `(b or a)`.

raven-toolbox's `diff_models` is weaker here. `_normalise_gpr` (`comparison/diff.py`) only
lowercases and collapses whitespace — its own docstring notes "a more robust comparator would parse
to a GPR AST and compare structures; this is the cheap heuristic". So it reports two logically
identical rules that differ only in operand order as a difference.

**Back-port.** Replace the string `_normalise_gpr` with an AST-based comparison. cobra already
exposes the parsed GPR (`Reaction.gpr` / `GPR`), and the repo has `gpr_to_dnf`
(`manipulation/expand.py`) — canonicalise by DNF-expanding, sorting the genes within each clause and
sorting the clauses, then compare. This is a small change and makes the two implementations agree.
Reference: `comparison/diffModels.m::canonicalGpr` in MATLAB RAVEN.

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
