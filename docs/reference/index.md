# Reference

Conceptual and API reference for raven-toolbox.

- **[Roadmap](roadmap.md)** — the phased development plan: what ships in which
  order, what each phase depends on, and the decisions still open. The item-level
  backlog behind it is [Open work](todo.md).

- **[RAVEN ↔ raven-toolbox migration map](migration.md)** — the function-by-function map
  from MATLAB RAVEN to raven-toolbox (and cobrapy where appropriate). Start here if you're
  porting RAVEN code.
- **[MATLAB RAVEN back-port proposals](matlab_raven_backports.md)** — improvements
  raven-toolbox makes that are candidates to back-port into the MATLAB toolbox.
- **[Improvements over RAVEN](improvements.md)** — the full catalogue of correctness /
  ergonomics improvements (the `IMPROVEMENTS.md` master list).
- **[API reference](api/index.md)** — every public function and class, generated from the
  docstrings.
- **[COBRA vs RAVEN comparison](cobra_raven_comparison.md)** — feature-by-feature comparison
  of the COBRA Toolbox and RAVEN Toolbox, identifying gaps and overlap.
- **[Evidence-aware transport scoring — design & plan](transport_evidence_scoring.md)** — a cross-repo
  (RAVEN + raven-toolbox) plan to replace the blanket transport penalty in localisation with
  transporter-evidence-weighted costs; carrier-general, organism-agnostic, local-binary-based.
- **[Sound reaction-level multi-localisation — design](multi_localization_design.md)** — why naive
  multi-localisation admits *dead* placements, and the ε-flux activity-coupling formulation
  `assign_compartments` uses to forbid them (solver-independent).

:::{admonition} Moved to raven-docs
:class: note

The YAML model format spec, tuned parameter defaults, and the CHRR/ACHR flux
sampling algorithm reference now live on
[raven-docs](https://github.com/edkerk/raven-docs) — they document
both RAVEN MATLAB and raven-toolbox side by side, not just this package. The
gap-filling literature review was dropped (MATLAB-scoped, and never
discussed raven-toolbox's own gap-filling module).
:::

```{toctree}
:hidden:

roadmap
migration
matlab_raven_backports
improvements
api/index
cobra_raven_comparison
transport_evidence_scoring
multi_localization_design
```
