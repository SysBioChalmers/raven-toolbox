"""Replace gene identifiers in a model — port of RAVEN's ``renameModelGenes``.

Builds an old-id → new-id mapping from a gene table (e.g. from
:mod:`raven_toolbox.reconstruction.genome_data`, once ported, or supplied
manually) and applies it via :func:`cobra.manipulation.rename_genes`, rather
than reimplementing RAVEN's own regex-based word-boundary GPR rewrite:
cobra's version parses the GPR into an AST and renames the matched gene
nodes directly, so it cannot suffer the partial-match failure mode
(`"gene1"` matching inside `"gene10"`) RAVEN's regex specifically guards
against, and correctly *merges* the underlying gene objects when a rename
target already exists elsewhere in the model (the GPR string still
references it as many times as before — e.g. ``"G1 or G1"`` — renaming
doesn't algebraically simplify a rule, only the gene identity is
deduplicated) — RAVEN's own version has no such check, and can leave a
model with two ``model.genes`` entries carrying the same id.
cobra also already keeps GPR syntax normalised on assignment, so
``removeUnnecessaryParentheses``/``standardizeGrRules`` have nothing to do
here either.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import cobra
import cobra.manipulation
import pandas as pd

__all__ = ["RenameGenesResult", "rename_model_genes"]


@dataclass
class RenameGenesResult:
    """Outcome of a gene-renaming pass.

    Parameters
    ----------
    renamed:
        Number of model genes for which a mapping was found and applied.
    unmapped:
        Model gene ids with no entry in ``from_col`` (or an empty ``to_col``),
        left unchanged. Sorted, unique.
    """

    renamed: int
    unmapped: list[str]


def rename_model_genes(
    model: cobra.Model,
    gene_table: pd.DataFrame | str | Path,
    from_col: str,
    to_col: str,
) -> RenameGenesResult:
    """Rename ``model``'s genes in place using a mapping table.

    Parameters
    ----------
    model:
        Model to rename genes in.
    gene_table:
        Either a DataFrame, or a path to a tab-delimited file to load one from.
    from_col:
        Column whose values match the ids currently in ``model.genes``.
    to_col:
        Column whose values will replace them.

    Returns
    -------
    RenameGenesResult

    Raises
    ------
    ValueError
        If ``from_col`` or ``to_col`` is not a column of ``gene_table``.
    """
    if isinstance(gene_table, (str, Path)):
        gene_table = pd.read_csv(gene_table, sep="\t")

    for col in (from_col, to_col):
        if col not in gene_table.columns:
            raise ValueError(
                f"Column {col!r} not found in gene_table. "
                f"Available columns: {list(gene_table.columns)}"
            )

    # First non-empty (from, to) pair wins for a duplicated from_col value;
    # rows missing either side are skipped.
    mapping: dict[str, str] = {}
    for from_val, to_val in zip(gene_table[from_col], gene_table[to_col], strict=True):
        if pd.isna(from_val) or pd.isna(to_val):
            continue
        old = str(from_val).strip()
        new = str(to_val).strip()
        if old and new and old not in mapping:
            mapping[old] = new

    model_gene_ids = [g.id for g in model.genes]
    rename_dict = {g: mapping[g] for g in model_gene_ids if g in mapping}
    unmapped = sorted(g for g in model_gene_ids if g not in mapping)

    cobra.manipulation.rename_genes(model, rename_dict)

    return RenameGenesResult(renamed=len(rename_dict), unmapped=unmapped)
