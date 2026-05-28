# ravengem

[![CI](https://github.com/SysBioChalmers/ravengem/actions/workflows/ci.yml/badge.svg)](https://github.com/SysBioChalmers/ravengem/actions/workflows/ci.yml)

**Reconstruction, Analysis and Visualisation of Metabolic Networks — in Python.**

`ravengem` is the Python counterpart of the
[RAVEN Toolbox 2](https://github.com/SysBioChalmers/RAVEN) (MATLAB). It builds on
[**cobrapy**](https://github.com/opencobra/cobrapy) for everything cobrapy already does
well (simulation, standard analyses, SBML I/O, model manipulation) and adds the
functionality that's unique to RAVEN:

* **De novo reconstruction** from KEGG and protein homology (BLAST / DIAMOND).
* **Context-specific models** from omics data via **tINIT / ftINIT**, with task-aware
  gap-filling and the linear-merge MILP reduction.
* **Metabolic-task** validation (`check_tasks`, `fitTasks`).
* **Connectivity gap-filling** against template models.
* **Omics integration** — Human Protein Atlas (proteomics + RNA-seq) ingestion.
* **Sub-cellular localisation** prediction by MILP, with partial-update mode and
  pluggable predictors (WoLF PSORT, DeepLoc, …).
* **N-model comparison**; **reporter metabolites**; **FSEOF**; **flux sampling**.
* **YAML I/O** following the cobra standard, plus geckopy's `ec-*` enzyme-constrained
  fields. **SIF** export. **RAVEN-style Excel** export.

The status of every RAVEN function (ported, cheatsheet-mapped to cobra, or explicitly
not ported) is documented function-by-function in
**[docs/raven_migration.md](docs/raven_migration.md)**.

## Design principle

The canonical in-memory object is always a [`cobra.Model`](https://cobrapy.readthedocs.io).
There is no parallel RAVEN struct, no `ravenCobraWrapper`-style adapter. RAVEN-specific
fields that cobra doesn't model natively (`rxnMiriams`, `metDeltaG`,
`rxnConfidenceScores`, …) live in cobra's `annotation` / `notes` dictionaries. This
avoids duplicating cobra's data model and keeps ravengem interoperable with the wider
COBRA ecosystem.

## Status

ravengem has been validated against MATLAB RAVEN on **Human-GEM** (5 Hart2015 cell-line
models, Jaccard 0.975–0.980 — see [docs/humangem_validation.md](docs/humangem_validation.md)).
The functional scope of the original RAVEN toolbox is covered with two principled
omissions:

* **MetaCyc-based reconstruction** is not implemented and is flagged for removal from
  MATLAB RAVEN as well — see [IMPROVEMENTS.md](IMPROVEMENTS.md) under `R-MetaCyc`.
* **Dynamic FBA** is not implemented — several maintained Python packages already cover
  it ([`dfba`](https://pypi.org/project/dfba/), [`reframed`](https://pypi.org/project/reframed/),
  [`mewpy`](https://pypi.org/project/mewpy/)).

What's still open is catalogued in **[docs/todo.md](docs/todo.md)** (visualisation / Phase
6 is the main item).

## Installation (development)

```bash
git clone https://github.com/SysBioChalmers/ravengem
cd ravengem
pip install -e ".[dev]"
```

ravengem requires Python ≥ 3.11. Genome-scale (f)tINIT MILPs currently require **Gurobi**
([details on solver portability](docs/init_solver_benchmark.md)); toy and unit-test work
runs on the open-source GLPK.

## Documentation

See **[docs/README.md](docs/README.md)** for the documentation index.

## Relationship to MATLAB RAVEN

`ravengem` is a derivative work and is released under the same **GPL-3.0-or-later**
license. If you use it in scientific work, please cite the RAVEN 2 paper:

> Wang H, Marcišauskas S, Sánchez BJ, Domenzain I, Hermansson D, Agren R, Nielsen J,
> Kerkhoven EJ. (2018) RAVEN 2.0: A versatile toolbox for metabolic network
> reconstruction and a case study on *Streptomyces coelicolor*. PLoS Comput Biol 14(10):
> e1006541.

## License

[GPL-3.0-or-later](LICENSE)
