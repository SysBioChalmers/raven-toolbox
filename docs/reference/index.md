# Reference

Conceptual and API reference for raven-toolbox.

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

migration
matlab_raven_backports
improvements
api/index
cobra_raven_comparison
```
