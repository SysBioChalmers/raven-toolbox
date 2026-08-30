# MATLAB RAVEN ↔ raven-toolbox — differences

Record of functionality in MATLAB RAVEN that raven-toolbox deliberately does not have, plus
raven-toolbox functionality still to be back-ported into MATLAB RAVEN.

---

## Back-ported but diverged: functionality-constrained compartment-assignment MILP

raven-toolbox's `localization/` hosts **two** compartment-assignment algorithms:
`predict_localization` (score-driven MILP) and **`assign_compartments`** (`localization/certify.py`) —
the *functionality-constrained* method, seeded by the port of the retired `edkerk/assignCompartments`
repo. Over the score-driven version it adds a **biomass/growth floor** enforced by certification:
placement is decided by a flux-free score MILP and the result is confirmed by a real FBA on the
materialised model. It also adds optional **gap-fill coupling** (universal-DB candidates added only when
biomass feasibility needs them) and **sound reaction-level multi-localisation** (a second compartment is
kept only if a loopless FVA on the materialised model shows it carries real flux — design in
[multi_localization_design.md](multi_localization_design.md)).

MATLAB RAVEN's default branch has no equivalent (`core/predictLocalization.m` is a
*simulated-annealing* heuristic: one gene → one compartment, no biomass constraint, no flux gating),
**but `develop3` already has `localization/assignCompartments.m`** — a port of raven-toolbox's
*earlier* design, from around the original `assign_compartments` (PR #58): a single MILP with the
biomass floor and flux gating fused directly into the placement problem (`bigM`-gated flux variables,
`minGrowth`), mono-localisation only. It predates raven-toolbox's rework (PR #62 onward) to a
flux-free placement MILP + separate real-FBA certification, and has none of the certification/feedback
loop, gap-fill coupling, or sound multi-localisation described above. The two `assign_compartments` are
no longer equivalent — the MATLAB side needs a re-sync, not a fresh port.

**Re-sync plan.** Rework `assignCompartments.m` to match the current design: drop the fused
biomass/flux-gating constraints from the placement MILP, add a separate real-FBA certification pass
(`optimizeProb` on the materialised model) plus the confinement-repair and feedback-loop steps, and the
multi-localisation flux-activity coupling ([multi_localization_design.md](multi_localization_design.md)
path 2). Verify the `.m` filename does not clash with a COBRA Toolbox function before committing (it
already exists on `develop3`, so this is a rework in place, not a new file). Tests under `testing/`
mirroring `tests/test_localization_certify.py`. Reference implementation:
`src/raven_toolbox/localization/certify.py` (the algorithm) and
`src/raven_toolbox/localization/assign.py` (`AssignmentProposal` / `apply_assignment` materialisation).

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
