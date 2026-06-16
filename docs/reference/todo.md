# Open work

What's still on the books. See [raven_migration.md](migration.md) for the function-by-
function port status (most of it is done); see [IMPROVEMENTS.md](improvements.md) for the
catalogue of raven-toolbox improvements that should also be back-ported into MATLAB RAVEN.

## Smaller items

* [known_issues.md](known_issues.md) — backlog of low-priority edge cases / robustness gaps /
  dead code from the full-codebase review. None affects correctness on well-formed inputs.
* [IMPROVEMENTS.md](improvements.md) — items marked 💡 *proposed* are candidates to
  implement (and back-port).

## Not planned

* **Metabolomics-based scoring in tINIT / ftINIT.** Intentionally not implemented;
  [`init.ftinit`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/ftinit.py)
  raises `NotImplementedError` for a non-empty `metabolomics` argument (the most intricate
  MILP piece, for the least-used input).

## Upstream blockers (not raven-toolbox work, but worth tracking)

* `optlang.hybrid_interface.Configuration.clone()` bug — blocks HiGHS at any scale (CI catches
  it in `tests/test_init_solvers.py`).
* GLPK's MIP solve ignores `configuration.timeout` at genome scale — blocks GLPK on large MILPs.
* Both documented in [init_solver_benchmark.md](../studies/init_solver_benchmark.md) with concrete fix
  suggestions.
