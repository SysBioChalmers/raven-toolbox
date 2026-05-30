"""Batch curation of metabolites / reactions / genes from data tables.

Port of yeast-GEM's MATLAB ``curateMetsRxnsGenes`` into a generic
DataFrame-driven engine. Other GEM projects (Human-GEM, custom
reconstructions, …) can use the same machinery with their own TSV
layouts; the only required pieces are the data tables and (optionally)
project-specific id prefixes for fresh metabolites and reactions.

Public API:

* :func:`batch_curate` — entrypoint taking pandas DataFrames.
* :func:`batch_curate_from_tsv` — file-path convenience wrapper.
* :class:`CurationResult` — record of what was added / updated.

Schema (mirrors yeast-GEM's ``data/modelCuration/template/`` layout):

- **mets_df**: ``metNames, comps, formula, charge, inchi, metNotes``
  + any number of MIRIAM-annotation columns. Match key is
  ``(name, comp)``.
- **genes_df**: ``genes, geneShortNames`` + MIRIAM columns. Match key
  is ``genes``.
- **rxns_df**: ``rxnNames, grRules, lb, ub, rev, subSystems, eccodes,
  rxnNotes, rxnReferences, rxnConfidenceScores`` + MIRIAM columns.
  Match key is the reaction's *stoichiometry* — same metabolites and
  coefficients ⇒ same reaction.
- **rxns_coeffs_df**: ``rxnNames, metNames, comps, coefficient``. One
  row per ``(reaction, metabolite)`` pair. The ``rxnNames`` column
  links each coefficient back to a row in ``rxns_df``. An optional
  ``index`` first column from the legacy yeast-GEM schema is silently
  ignored.

Everything after the core columns in any of the four tables is
interpreted as a MIRIAM annotation: the column header becomes the
namespace key (``met.annotation[<header>] = <cell>``).
"""
from raven_python.curation.batch import (
    DEFAULT_CORE_GENE_COLUMNS,
    DEFAULT_CORE_MET_COLUMNS,
    DEFAULT_CORE_RXN_COEFFS_COLUMNS,
    DEFAULT_CORE_RXN_COLUMNS,
    CurationResult,
    batch_curate,
    batch_curate_from_tsv,
)

__all__ = [
    "DEFAULT_CORE_GENE_COLUMNS",
    "DEFAULT_CORE_MET_COLUMNS",
    "DEFAULT_CORE_RXN_COEFFS_COLUMNS",
    "DEFAULT_CORE_RXN_COLUMNS",
    "CurationResult",
    "batch_curate",
    "batch_curate_from_tsv",
]
