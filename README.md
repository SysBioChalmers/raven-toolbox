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
  pluggable evidence sources (DeepLoc 2, MULocDeep, COMPARTMENTS, UniProt, …).
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
models, Jaccard 0.975–0.980 — see [the Human-GEM validation study](https://github.com/edkerk/raven-docs/blob/main/docs/parameter-tuning/studies/humangem-validation.md) on raven-docs).
The functional scope of the original RAVEN toolbox is covered with three principled
omissions, all deliberately out of scope rather than pending work:

* **MetaCyc-based reconstruction** is not implemented and is flagged for removal from
  MATLAB RAVEN as well — see [IMPROVEMENTS.md](IMPROVEMENTS.md) under `R-MetaCyc`.
* **Dynamic FBA** is not implemented — several maintained Python packages already cover
  it ([`dfba`](https://pypi.org/project/dfba/), [`reframed`](https://pypi.org/project/reframed/),
  [`mewpy`](https://pypi.org/project/mewpy/)).
* **Metabolomics-based scoring in tINIT / ftINIT** is not implemented — passing a
  non-empty `metabolomics` argument raises `NotImplementedError`.

## Installation (development)

```bash
git clone https://github.com/SysBioChalmers/raven-toolbox
cd raven-toolbox
pip install -e ".[dev]"
```

raven-toolbox requires Python ≥ 3.11. Genome-scale (f)tINIT MILPs currently require **Gurobi**
([details on solver portability](https://github.com/edkerk/raven-docs/blob/main/docs/parameter-tuning/studies/init-solver-benchmark.md)
on raven-docs); toy and unit-test work runs on the open-source GLPK.

### External command-line tools (BLAST, DIAMOND, HMMER, MAFFT, CD-HIT)

Some workflows call external tools. **For most users there is nothing to do** —
raven-toolbox downloads each tool it needs automatically the first time it's used.

Which tools a workflow uses:

| Workflow | Tools |
|---|---|
| Homology-based reconstruction | `blastp` + `makeblastdb`, or `diamond` |
| KEGG HMM query (`getKEGGModelForOrganism`) | `hmmsearch` |
| Building the KEGG HMM libraries (maintainers) | `hmmbuild`, `mafft`, `cd-hit` |

Optional:

- **Fetch them up front** instead of on first use:
  ```bash
  raven-toolbox-binaries --set runtime   # blastp, makeblastdb, diamond, hmmsearch
  raven-toolbox-binaries --set build     # hmmbuild, mafft, cd-hit
  ```
- **Use your own install** — if a tool is on your `PATH` it's used instead of a
  download, e.g. `conda install -c bioconda blast diamond hmmer mafft cd-hit`.
- **Disable automatic downloads** (air-gapped / conda-only setups): set
  `RAVEN_PYTHON_AUTOFETCH=0`.

**Windows:** homology reconstruction and the KEGG species model work as-is. To
*build* the KEGG HMM libraries (needs MAFFT/CD-HIT), use **WSL2**.

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
