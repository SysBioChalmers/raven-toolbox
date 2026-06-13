# Open work

What's still on the books. See [raven_migration.md](migration.md) for the function-by-
function port status (most of it is done); see [IMPROVEMENTS.md](improvements.md) for the
catalogue of raven-toolbox improvements that should also be back-ported into MATLAB RAVEN.

## Major

### Metabolomics-based scoring in ftINIT

The metabolomics-detected metabolite production-reward block in [`init.ftinit`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/ftinit.py)
currently raises `NotImplementedError` if a non-empty `metabolomics` argument is passed. The
linear merge eliminates degree-2 detected metabolites, so it needs producer-group-mapping +
negative-producer force-flux constraints — the most intricate MILP piece, for the least-used
input. Worth doing only when a real user request lands.

## Infrastructure

* **Binary ZIP releases** for BLAST/DIAMOND (Phase 3a). The runtime resolver in
  [`binaries.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/binaries.py) is ready; the registry is empty until ZIPs are
  published as GitHub release assets.
* **KEGG data artefact releases.** See [maintaining_kegg_data.md](../maintenance/maintaining_kegg_data.md).

## Smaller items

* [known_issues.md](known_issues.md) — backlog of low-priority edge cases / robustness gaps /
  dead code from the full-codebase review. None affects correctness on well-formed inputs.
* [IMPROVEMENTS.md](improvements.md) — items marked 💡 *proposed* are candidates to
  implement (and back-port).

## Upstream blockers (not raven-toolbox work, but worth tracking)

* `optlang.hybrid_interface.Configuration.clone()` bug — blocks HiGHS at any scale (CI catches
  it in `tests/test_init_solvers.py`).
* GLPK's MIP solve ignores `configuration.timeout` at genome scale — blocks GLPK on large MILPs.
* Both documented in [init_solver_benchmark.md](../studies/init_solver_benchmark.md) with concrete fix
  suggestions.
