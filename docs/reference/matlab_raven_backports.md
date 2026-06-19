# MATLAB RAVEN ↔ raven-toolbox — differences

Record of functionality in MATLAB RAVEN that raven-toolbox deliberately does not have.
All planned backports from raven-toolbox into MATLAB RAVEN have been completed.

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
