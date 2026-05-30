# Installation

raven-python requires **Python ≥ 3.11**.

## From a source checkout (development)

```bash
git clone https://github.com/SysBioChalmers/raven-python
cd raven-python
pip install -e ".[dev]"
```

## Optional extras

The core install pulls in cobra, numpy, pandas, scipy, ruamel.yaml, requests and tqdm.
Some features need extra packages, exposed as optional extras:

| Extra | Pulls in | Needed for |
| --- | --- | --- |
| `excel` | `openpyxl` | {func}`raven_python.io.excel.export_to_excel` |
| `plotting` | `matplotlib` | visualisation (Phase 6, in progress) |
| `dev` | `pytest`, `pytest-cov`, `ruff` | running the test-suite / linting |

```bash
pip install -e ".[excel,plotting,dev]"
```

## Solvers

Linear and most MILP work runs on the open-source **GLPK** that ships with cobra.
Genome-scale **(f)tINIT** MILPs currently require **Gurobi** for tractable solve times —
see the [solver benchmark](studies/init_solver_benchmark.md) for the Gurobi vs HiGHS vs GLPK
comparison and portability notes.

## External binaries

Homology and KEGG-HMM reconstruction shell out to **BLAST+ / DIAMOND / HMMER**. These are
not Python packages; raven-python resolves them from a version-pinned, SHA256-verified
release registry (see {mod}`raven_python.binaries`). The
[binary maintenance guide](maintenance/maintaining_binaries.md) covers building and
publishing those release ZIPs.
