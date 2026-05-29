"""Loaders for gene → compartment localisation predictors (WoLF PSORT, DeepLoc, …).

The localisation algorithm in :mod:`raven_python.localization.predict` consumes a
*gene × compartment* score table (:class:`LocalizationScores`) where higher = stronger
evidence. Each predictor produces this differently; loaders here normalise them. The
format is open — a user can build a :class:`LocalizationScores` from any source by
constructing the :class:`pandas.DataFrame` directly.

Each loader normalises each gene's row so the best compartment is 1.0 (RAVEN's
``parseScores`` convention), which lets transport costs be set on a comparable scale.
"""
from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass
class LocalizationScores:
    """Per-gene compartment scores. ``df`` is indexed by ``gene_id`` with one column per
    compartment id; values are floats (higher = stronger evidence for that compartment).

    Genes absent from ``df`` and NaN entries are treated as "no signal" by
    :func:`raven_python.localization.predict_localization` (uniform prior contribution).
    """

    df: pd.DataFrame

    def __post_init__(self) -> None:
        if not isinstance(self.df.index, pd.Index) or self.df.index.name not in (None, "gene_id"):
            # accept but normalise
            self.df = self.df.copy()
            self.df.index.name = "gene_id"

    @property
    def genes(self) -> list[str]:
        return list(self.df.index)

    @property
    def compartments(self) -> list[str]:
        return list(self.df.columns)

    def with_compartments(self, mapping: Mapping[str, str]) -> LocalizationScores:
        """Rename compartment columns via ``{old: new}`` (e.g. predictor labels →
        model compartments). Unmapped columns are kept; multiple sources can be merged
        with ``df.combine_first`` afterwards."""
        return LocalizationScores(self.df.rename(columns=dict(mapping)))


# ----------------------------------------------------------------------- WoLF PSORT

# WoLF PSORT summary lines look like:
#     PROTEIN_ID cyto 13, nucl 7, mito 4
# with comments starting '#' and noisy 'treating ...' lines (which we drop).
_WOLF_COMMA = re.compile(r"[,]\s*")


def load_wolfpsort(path: str | Path) -> LocalizationScores:
    """Parse WoLF PSORT summary output (``runWolfPsortSummary``) into a normalised
    :class:`LocalizationScores`. Rows like ``PROT: treating N X's as ...`` are skipped."""
    rows: dict[str, dict[str, float]] = {}
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "treating " in line:
            continue
        tokens = _WOLF_COMMA.sub(" ", line).split()
        if len(tokens) < 3 or (len(tokens) - 1) % 2 != 0:
            continue  # malformed; skip
        gene = tokens[0]
        comp_scores: dict[str, float] = {}
        for comp, score in zip(tokens[1::2], tokens[2::2], strict=True):
            try:
                comp_scores[comp] = float(score)
            except ValueError:
                continue
        if comp_scores:
            rows[gene] = comp_scores
    df = pd.DataFrame.from_dict(rows, orient="index").fillna(0.0)
    df.index.name = "gene_id"
    return _normalise_rows(LocalizationScores(df))


# ----------------------------------------------------------------------- DeepLoc

def load_deeploc(path: str | Path) -> LocalizationScores:
    """Parse DeepLoc 2 CSV output into a normalised :class:`LocalizationScores`.

    DeepLoc 2's per-protein CSV has columns ``Protein_ID, Localizations, Signals,
    <Compartment1>, <Compartment2>, ...`` where columns 4+ are per-class probabilities.
    The first three metadata columns are dropped; the rest become compartment columns.
    """
    df = pd.read_csv(path)
    if df.shape[1] < 4:
        raise ValueError(f"{path}: expected ≥4 columns from DeepLoc, got {list(df.columns)}")
    gene_col = df.columns[0]            # Protein_ID
    comp_cols = list(df.columns[3:])    # cols 0-2 are Protein_ID/Localizations/Signals metadata
    scores = df.set_index(gene_col)[comp_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    scores.index.name = "gene_id"
    return _normalise_rows(LocalizationScores(scores))


# ----------------------------------------------------------------------- helpers

def _normalise_rows(s: LocalizationScores) -> LocalizationScores:
    """Per-gene row normalisation: best compartment → 1.0 (RAVEN's parseScores convention).

    Rows whose max is ≤0 are left unscaled (no positive evidence to normalise against).
    """
    df = s.df.copy()
    row_max = df.max(axis=1)
    safe = row_max > 0
    df.loc[safe] = df.loc[safe].div(row_max[safe], axis=0)
    return LocalizationScores(df)
