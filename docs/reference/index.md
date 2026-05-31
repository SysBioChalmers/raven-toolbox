# Reference

Conceptual and API reference for raven-python.

- **[RAVEN ↔ raven-python migration map](migration.md)** — the function-by-function map
  from MATLAB RAVEN to raven-python (and cobrapy where appropriate). Start here if you're
  porting RAVEN code.
- **[YAML model format](yaml_format.md)** — the shared YAML schema produced and consumed
  by cobrapy, raven-python, and RAVEN MATLAB, with a fully-annotated example and the
  field-order / quoting rules.
- **[MATLAB RAVEN back-port proposals](matlab_raven_backports.md)** — improvements
  raven-python makes that are candidates to back-port into the MATLAB toolbox.
- **[Improvements over RAVEN](improvements.md)** — the full catalogue of correctness /
  ergonomics improvements (the `IMPROVEMENTS.md` master list).
- **[API reference](api/index.md)** — every public function and class, generated from the
  docstrings.

```{toctree}
:hidden:

migration
yaml_format
matlab_raven_backports
improvements
api/index
```
