# Reference

Conceptual and API reference for raven-toolbox.

- **[RAVEN ↔ raven-toolbox migration map](migration.md)** — the function-by-function map
  from MATLAB RAVEN to raven-toolbox (and cobrapy where appropriate). Start here if you're
  porting RAVEN code.
- **[YAML model format](yaml_format.md)** — the shared YAML schema produced and consumed
  by cobrapy, raven-toolbox, and RAVEN MATLAB, with a fully-annotated example and the
  field-order / quoting rules.
- **[MATLAB RAVEN back-port proposals](matlab_raven_backports.md)** — improvements
  raven-toolbox makes that are candidates to back-port into the MATLAB toolbox.
- **[Improvements over RAVEN](improvements.md)** — the full catalogue of correctness /
  ergonomics improvements (the `IMPROVEMENTS.md` master list).
- **[API reference](api/index.md)** — every public function and class, generated from the
  docstrings.
- **[Gap-filling algorithms — literature review](gap_filling_algorithms.md)** — survey of
  published gap-filling methods with formulation, trade-offs, and implementation recommendations.
- **[COBRA vs RAVEN comparison](cobra_raven_comparison.md)** — feature-by-feature comparison
  of the COBRA Toolbox and RAVEN Toolbox, identifying gaps and overlap.

```{toctree}
:hidden:

migration
yaml_format
matlab_raven_backports
improvements
api/index
gap_filling_algorithms
cobra_raven_comparison
```
