# Task checking parameter benchmarks

Functions: `raven_toolbox.tasks.check.check_tasks`,
`raven_toolbox.tasks.check.find_task_essential_reactions`

Date: 2026-06-20.

---

## `close_boundaries` — whether to close all exchange reactions before evaluating tasks

**Parameter:** `close_boundaries=True` (Python default, MATLAB implied)

When `True`, all exchange reaction bounds are set to zero before evaluating each
metabolic task. This ensures tasks are evaluated against the model's internal
biochemistry only, not against artefactual uptake of media components.

Setting `False` would allow arbitrary exchange fluxes, making tasks trivially
satisfiable by importing the required metabolites directly from the boundary.

**Decision: ✓ keep `True`.** This is the only semantically correct default for
task-based model validation; `False` would make tasks uninformative.

---

## `find_task_essential_reactions` — `tol`

**Parameter:** `tol=1e-8` (Python default, no direct MATLAB equivalent)

`tol` is the minimum task objective flux that must be achieved after deleting a
candidate essential reaction. A reaction is essential if deleting it causes the
task flux to fall below `tol`.

`1e-8` matches Gurobi's default primal feasibility tolerance. Setting `tol` too
high (e.g., `0.001`) would flag reactions as non-essential when the task objective
can still be satisfied at a tiny but numerically non-zero flux.

**Decision: ✓ keep `1e-8`.** Consistent with solver tolerances. No empirical test
needed — this is a threshold parameter whose correct value is determined by the
solver tolerance, not by biological considerations.
