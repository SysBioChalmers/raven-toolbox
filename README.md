# raven-toolbox

[![CI](https://github.com/SysBioChalmers/raven-toolbox/actions/workflows/ci.yml/badge.svg)](https://github.com/SysBioChalmers/raven-toolbox/actions/workflows/ci.yml)

**Reconstruction, Analysis and Visualisation of Metabolic Networks — in Python.**

`raven-toolbox` is the Python counterpart of the
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
**[docs/raven_migration.md](docs/reference/migration.md)**.

## Design principle

The canonical in-memory object is always a [`cobra.Model`](https://cobrapy.readthedocs.io).
There is no parallel RAVEN struct, no `ravenCobraWrapper`-style adapter. RAVEN-specific
fields that cobra doesn't model natively (`rxnMiriams`, `metDeltaG`,
`rxnConfidenceScores`, …) live in cobra's `annotation` / `notes` dictionaries. This
avoids duplicating cobra's data model and keeps raven-toolbox interoperable with the wider
COBRA ecosystem.

## Status

raven-toolbox has been validated against MATLAB RAVEN on **Human-GEM** (5 Hart2015 cell-line
models, Jaccard 0.975–0.980 — see [docs/humangem_validation.md](docs/studies/humangem_validation.md)).
The functional scope of the original RAVEN toolbox is covered with two principled
omissions:

* **MetaCyc-based reconstruction** is not implemented and is flagged for removal from
  MATLAB RAVEN as well — see [IMPROVEMENTS.md](IMPROVEMENTS.md) under `R-MetaCyc`.
* **Dynamic FBA** is not implemented — several maintained Python packages already cover
  it ([`dfba`](https://pypi.org/project/dfba/), [`reframed`](https://pypi.org/project/reframed/),
  [`mewpy`](https://pypi.org/project/mewpy/)).

### Not yet implemented

Planned or partial functionality still on the books (full detail in
**[docs/todo.md](docs/reference/todo.md)**):

- [ ] **Metabolomics-based scoring in tINIT / ftINIT** — passing a non-empty `metabolomics`
  argument currently raises `NotImplementedError`.
- [ ] **Published binary release ZIPs** (BLAST / DIAMOND / HMMER) — the resolver in
  `binaries.py` is ready; the registry is empty until the ZIPs are published as release assets.
- [ ] **Published KEGG data-artefact releases** — the build pipeline exists; the artefact
  bundle is not published yet.

The two principled omissions above (MetaCyc reconstruction, dynamic FBA) are **not** on this
list — they are deliberately out of scope, not pending work.

## Installation (development)

```bash
git clone https://github.com/SysBioChalmers/raven-toolbox
cd raven-toolbox
pip install -e ".[dev]"
```

raven-toolbox requires Python ≥ 3.11. Genome-scale (f)tINIT MILPs currently require **Gurobi**
([details on solver portability](docs/studies/init_solver_benchmark.md)); toy and unit-test work
runs on the open-source GLPK.

## Documentation

The documentation is built with Sphinx (MyST Markdown); the source lives in
[docs/](docs/) — see [docs/README.md](docs/README.md) for the layout and local-build
instructions. (A hosted ReadTheDocs site is not yet published.)

## Relationship to MATLAB RAVEN

`raven-toolbox` is an independent Python reimplementation of the
[RAVEN Toolbox 2](https://github.com/SysBioChalmers/RAVEN), released under the permissive
**MIT** license. If you use it in scientific work, please cite the RAVEN 2 paper:

> Wang H, Marcišauskas S, Sánchez BJ, Domenzain I, Hermansson D, Agren R, Nielsen J,
> Kerkhoven EJ. (2018) RAVEN 2.0: A versatile toolbox for metabolic network
> reconstruction and a case study on *Streptomyces coelicolor*. PLoS Comput Biol 14(10):
> e1006541.

## License

[MIT](LICENSE)
