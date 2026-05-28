# Open work

What's still on the books. See [raven_migration.md](raven_migration.md) for the function-by-
function port status (most of it is done); see [IMPROVEMENTS.md](../IMPROVEMENTS.md) for the
catalogue of ravengem improvements that should also be back-ported into MATLAB RAVEN.

## Major

### Visualization (`plotting/`)

Not started. RAVEN has limited plotting (`drawMap` etc., MATLAB-bound). For ravengem the most
useful targets are:

* Pathway maps / Escher integration for context-specific models.
* Omics overlay (gene-score / expression heatmaps on the reaction set).
* Flux distribution overlays.

cobrapy + Escher already covers a lot here — likely a thin integration layer rather than a port.

### Metabolomics-based scoring in ftINIT

The metabolomics-detected metabolite production-reward block in [`init.ftinit`](../src/ravengem/init/ftinit.py)
currently raises `NotImplementedError` if a non-empty `metabolomics` argument is passed. The
linear merge eliminates degree-2 detected metabolites, so it needs producer-group-mapping +
negative-producer force-flux constraints — the most intricate MILP piece, for the least-used
input. Worth doing only when a real user request lands.

## Infrastructure

* **CI for Phase 1 / foundation.** Pytest + lint workflow on GitHub Actions. Currently runs
  locally only.
* **Binary ZIP releases** for BLAST/DIAMOND (Phase 3a). The runtime resolver in
  [`binaries.py`](../src/ravengem/binaries.py) is ready; the registry is empty until ZIPs are
  published as GitHub release assets.
* **KEGG data artefact releases.** See [maintaining_kegg_data.md](maintaining_kegg_data.md).

## Smaller items

* [known_issues.md](known_issues.md) — backlog of low-priority edge cases / robustness gaps /
  dead code from the full-codebase review. None affects correctness on well-formed inputs.
* [IMPROVEMENTS.md](../IMPROVEMENTS.md) — items marked 💡 *proposed* are candidates to
  implement (and back-port).

## Upstream blockers (not ravengem work, but worth tracking)

* `optlang.hybrid_interface.Configuration.clone()` bug — blocks HiGHS at any scale (CI catches
  it in `tests/test_init_solvers.py`).
* GLPK's MIP solve ignores `configuration.timeout` at genome scale — blocks GLPK on large MILPs.
* Both documented in [init_solver_benchmark.md](init_solver_benchmark.md) with concrete fix
  suggestions.
