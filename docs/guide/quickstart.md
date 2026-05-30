# Quickstart

This page gets you from a fresh install to loading a model and running the two headline
workflows: a de-novo draft reconstruction and a context-specific extraction.

See [Installation](../installation.md) for the full dependency / solver / binary matrix.

```bash
pip install -e ".[dev]"
```

## Load and save a model

The canonical in-memory object is a {class}`cobra.Model`, so everything cobra can do is
available unchanged. raven-python adds RAVEN's YAML reader/writer (cobra-standard layout
plus RAVEN/GECKO side-fields preserved on `notes`):

```python
from raven_python.io import read_yaml_model, write_yaml_model

model = read_yaml_model("model.yml")     # transparently handles .yml.gz
print(model.summary())

write_yaml_model(model, "out.yml", sort_ids=True)
```

Other I/O and structural edits live in {mod}`raven_python.io` and
{mod}`raven_python.manipulation` — see the
[I/O & manipulation guide](io_and_manipulation.md).

## De-novo reconstruction (homology)

Build a draft for a new organism from a curated template model and a BLAST/DIAMOND
ortholog search:

```python
from raven_python.reconstruction.homology import get_model_from_homology

draft = get_model_from_homology(template_model, ortholog_hits)
```

The KEGG route ({func}`raven_python.reconstruction.kegg.get_kegg_model_for_organism`) and
the BLAST/DIAMOND helpers are covered in the
[reconstruction guide](reconstruction.md).

## Context-specific model (ftINIT)

Extract a tissue/condition-specific model from a reference GEM plus omics-derived gene
scores:

```python
from raven_python.omics import parse_hpa_rna, rna_gene_scores
from raven_python.init import ftinit

rna = parse_hpa_rna("rna_tissue.tsv")
scores = rna_gene_scores(reference_model, rna, tissue="liver")
context_model = ftinit(reference_model, scores)
```

ftINIT, tINIT, the scoring adapters and task-aware gap-filling are detailed in the
[context-specific modeling guide](context_specific.md). Genome-scale (f)tINIT currently
needs **Gurobi** — see the [solver benchmark](../studies/init_solver_benchmark.md).

:::{note}
The snippets above show the entry points; consult each capability guide and the
[API reference](../reference/api/index.md) for the full keyword arguments.
:::
