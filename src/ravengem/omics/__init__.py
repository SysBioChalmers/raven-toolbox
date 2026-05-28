"""Omics integration — HPA proteomics + RNA-seq parsing and gene-scoring adapters.

Entry point for tissue-specific (f)tINIT runs. See :mod:`ravengem.omics.hpa`.
"""
from ravengem.omics.hpa import (
    HPA_LEVEL_SCORES,
    HPAData,
    HPARnaData,
    hpa_gene_scores,
    parse_hpa,
    parse_hpa_rna,
    rna_gene_scores,
)

__all__ = [
    "HPA_LEVEL_SCORES",
    "HPAData",
    "HPARnaData",
    "hpa_gene_scores",
    "parse_hpa",
    "parse_hpa_rna",
    "rna_gene_scores",
]
