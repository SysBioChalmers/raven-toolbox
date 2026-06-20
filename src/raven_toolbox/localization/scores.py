"""Loaders for gene → compartment localisation predictors (DeepLoc 2, MULocDeep, COMPARTMENTS).

The localisation algorithm in :mod:`raven_toolbox.localization.predict` consumes a
*gene × compartment* score table (:class:`LocalizationScores`) where higher = stronger
evidence. Each predictor/database produces this differently; the loaders here normalise them.
The format is open — a user can build a :class:`LocalizationScores` from any source by
constructing the :class:`pandas.DataFrame` directly.

Predictors label compartments with their own names (``Mitochondrion``, ``Cytoplasm``, …).
Pass ``compartment_map`` (e.g. :data:`DEFAULT_COMPARTMENT_MAP`) to rename them to your model's
compartment ids and collapse synonyms; labels absent from the map are dropped. Without a map the
predictor's own labels are kept (use :meth:`LocalizationScores.with_compartments` to rename later).

Each loader normalises every gene's row so the best compartment is 1.0 (RAVEN's ``parseScores``
convention), which lets transport costs be set on a comparable scale.

Sequence-based predictors are external tools; run them yourself and pass their output files here.
Modern multi-label predictors (DeepLoc 2, MULocDeep) and the COMPARTMENTS evidence database
supersede older single-label callers, so no loader for those is provided.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

__all__ = [
    "LocalizationScores",
    "DEFAULT_COMPARTMENT_MAP",
    "load_deeploc",
    "load_mulocdeep",
    "load_compartments",
]


@dataclass
class LocalizationScores:
    """Per-gene compartment scores. ``df`` is indexed by ``gene_id`` with one column per
    compartment id; values are floats (higher = stronger evidence for that compartment).

    Genes absent from ``df`` and NaN entries are treated as "no signal" by
    :func:`raven_toolbox.localization.predict_localization` (uniform prior contribution).
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


# --------------------------------------------------------------------- compartment mapping

#: Default predictor/database compartment label → model compartment id. Keys are matched
#: case-insensitively. Tuned for yeast/fungal models (e.g. yeast-GEM codes); override per model.
#: Labels not present here (e.g. ``Plastid`` for fungi) are dropped when this map is used.
DEFAULT_COMPARTMENT_MAP: dict[str, str] = {
    "cytoplasm": "c",
    "cytosol": "c",
    "nucleus": "n",
    "nucleoplasm": "n",
    "mitochondrion": "m",
    "mitochondria": "m",
    "mitochondrial": "m",
    "peroxisome": "p",
    "endoplasmic reticulum": "er",
    "golgi apparatus": "g",
    "golgi": "g",
    "vacuole": "v",
    "lysosome/vacuole": "v",
    "lysosome": "v",
    "extracellular": "e",
    "extracellular space": "e",
    "extracellular region": "e",
    "secreted": "e",
    "cell membrane": "ce",
    "plasma membrane": "ce",
    "cell envelope": "ce",
    "lipid particle": "lp",
    "lipid droplet": "lp",
}


def _apply_map(df: pd.DataFrame, compartment_map: Mapping[str, str] | None) -> pd.DataFrame:
    """Rename compartment columns to model ids and collapse synonyms (max). Columns whose label
    is not in ``compartment_map`` are dropped. ``None`` leaves the columns untouched."""
    if not compartment_map:
        return df
    lower = {str(k).strip().lower(): v for k, v in compartment_map.items()}
    renamed = {col: lower[str(col).strip().lower()]
               for col in df.columns if str(col).strip().lower() in lower}
    sub = df[list(renamed)].rename(columns=renamed)
    if sub.shape[1] == 0:
        return sub
    # collapse duplicate target ids (e.g. Lysosome + Vacuole → v) taking the strongest evidence
    return sub.T.groupby(level=0).max().T


# --------------------------------------------------------------------- DeepLoc 2 / MULocDeep

def _load_wide(df: pd.DataFrame, *, id_column, compartment_map, source: str) -> LocalizationScores:
    """Engine for wide *id + per-compartment-probability* tables. Non-numeric metadata columns
    are detected and dropped automatically."""
    id_col = id_column if id_column is not None else df.columns[0]
    wide = df.set_index(id_col).apply(pd.to_numeric, errors="coerce").dropna(axis=1, how="all")
    if wide.shape[1] == 0:
        raise ValueError(f"{source}: no numeric compartment columns found in {list(df.columns)}")
    wide = _apply_map(wide.fillna(0.0), compartment_map)
    wide.index.name = "gene_id"
    return _normalise_rows(LocalizationScores(wide))


def load_deeploc(path: str | Path, *,
                 compartment_map: Mapping[str, str] | None = None) -> LocalizationScores:
    """Parse DeepLoc 2 CSV output into a normalised :class:`LocalizationScores`.

    DeepLoc 2's per-protein CSV has ``Protein_ID, Localizations, Signals`` then one probability
    column per compartment. The first column is the gene id; non-numeric metadata columns
    (``Localizations``, ``Signals``) are dropped automatically.
    """
    return _load_wide(pd.read_csv(path), id_column=None,
                      compartment_map=compartment_map, source=str(path))


def load_mulocdeep(path: str | Path, *,
                   compartment_map: Mapping[str, str] | None = None,
                   id_column: str | None = None,
                   sep: str | None = None) -> LocalizationScores:
    """Parse MULocDeep output into a normalised :class:`LocalizationScores`.

    Expects a wide table: a protein/gene id column plus one numeric probability column per major
    subcellular compartment (MULocDeep predicts 10). ``sep=None`` auto-detects the delimiter;
    pass ``id_column`` if the id is not the first column.
    """
    df = (pd.read_csv(path, sep=None, engine="python") if sep is None
          else pd.read_csv(path, sep=sep))
    return _load_wide(df, id_column=id_column, compartment_map=compartment_map, source=str(path))


# --------------------------------------------------------------------- COMPARTMENTS database

def load_compartments(path: str | Path, *,
                      compartment_map: Mapping[str, str] | None = None,
                      min_confidence: float = 0.0) -> LocalizationScores:
    """Parse a COMPARTMENTS (jensenlab.org) channel file into a normalised
    :class:`LocalizationScores`.

    Targets the tidy ``*_<channel>_full.tsv`` layout (no header): column 1 = protein/gene id,
    column 4 = GO compartment term name, last column = confidence score. Rows are aggregated per
    (gene, compartment) by the maximum confidence; ``min_confidence`` drops weak annotations
    (the integrated/text-mining channels score 0–5).
    """
    raw = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if raw.shape[1] < 5:
        raise ValueError(f"{path}: expected a COMPARTMENTS full TSV (>=5 columns), "
                         f"got {raw.shape[1]}")
    long = pd.DataFrame({
        "gene_id": raw[0].astype(str),
        "compartment": raw[3].astype(str),
        "score": pd.to_numeric(raw[raw.columns[-1]], errors="coerce"),
    }).dropna(subset=["score"])
    long = long[long["score"] >= min_confidence]
    wide = (long.pivot_table(index="gene_id", columns="compartment", values="score",
                             aggfunc="max").fillna(0.0))
    wide = _apply_map(wide, compartment_map)
    wide.index.name = "gene_id"
    return _normalise_rows(LocalizationScores(wide))


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
