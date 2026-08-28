# raven-toolbox

**Reconstruction, Analysis and Visualisation of Metabolic Networks — in Python.**

`raven-toolbox` is the Python counterpart of the
[RAVEN Toolbox 2](https://github.com/SysBioChalmers/RAVEN) (MATLAB). It builds on
[**cobrapy**](https://cobrapy.readthedocs.io) for everything cobrapy already does well
(simulation, standard analyses, SBML I/O, model manipulation) and adds the functionality
that is unique to RAVEN:

- **De novo reconstruction** from KEGG and protein homology (BLAST / DIAMOND).
- **Context-specific models** from omics data via **tINIT / ftINIT**, with task-aware
  gap-filling and the linear-merge MILP reduction.
- **Metabolic-task** validation (`check_tasks`, `find_task_essential_reactions`).
- **Connectivity gap-filling** against template models.
- **Omics integration** — Human Protein Atlas (proteomics + RNA-seq) ingestion.
- **Sub-cellular localisation** prediction by MILP, with partial-update mode and pluggable
  predictors (WoLF PSORT, DeepLoc, …).
- **N-model comparison**; **reporter metabolites**; **FSEOF**; **flux sampling**.
- **YAML I/O** following the cobra standard, plus geckopy's `ec-*` enzyme-constrained
  fields; **SIF** export; **RAVEN-style Excel** export.

:::{admonition} Design principle
:class: tip

The canonical in-memory object is always a {class}`cobra.Model`. There is no parallel
RAVEN struct and no `ravenCobraWrapper`-style adapter — RAVEN-specific fields that cobra
doesn't model natively live in cobra's `annotation` / `notes` dictionaries. This keeps
raven-toolbox interoperable with the wider COBRA ecosystem.
:::

## Where to start

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} 🚀 Quickstart
:link: guide/quickstart
:link-type: doc

Install, load a model, and run your first reconstruction / ftINIT extraction.
:::

:::{grid-item-card} 🧭 Coming from MATLAB RAVEN?
:link: reference/migration
:link-type: doc

The function-by-function map from RAVEN to raven-toolbox (and cobrapy).
:::

:::{grid-item-card} 📚 User guide
:link: guide/index
:link-type: doc

Task-oriented how-tos for each capability — reconstruction, tINIT, tasks, omics, …
:::

:::{grid-item-card} 🔍 API reference
:link: reference/api/index
:link-type: doc

Every public function and class, generated from the docstrings.
:::
::::

## Status

raven-toolbox has been validated against MATLAB RAVEN on **Human-GEM** (5 Hart2015 cell-line
models, Jaccard 0.975–0.980 — see
[the Human-GEM validation study](https://github.com/edkerk/raven-docs/blob/main/docs/parameter-tuning/studies/humangem-validation.md)
on raven-docs). The functional scope of the original toolbox is covered, with two principled
omissions: **MetaCyc-based reconstruction** (flagged for removal from MATLAB RAVEN too) and
**dynamic FBA** (well covered by other maintained Python packages). Candidates for back-porting
to MATLAB RAVEN are catalogued in [the improvements list](reference/improvements.md).

```{toctree}
:hidden:
:caption: Getting started

guide/quickstart
installation
```

```{toctree}
:hidden:
:caption: User guide

guide/index
```

```{toctree}
:hidden:
:caption: Reference

reference/index
```

```{toctree}
:hidden:
:caption: Maintenance

maintenance/parameter_defaults
```
