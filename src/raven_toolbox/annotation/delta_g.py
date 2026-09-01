"""Persist scalar annotations (Gibbs free energy, etc.) as a side-car CSV.

The committed SBML standard has no slot for thermodynamic ΔG values, so
projects like yeast-GEM keep them in a paired CSV instead. This module
provides the generic load / save round-trip — the file format is a two
column table ``<id>, <value>`` — and exposes the per-call schema so
existing projects can configure their conventions.

The functions operate on cobra ``.notes`` (string round-tripped through
SBML), keeping any other notes untouched.
"""
from __future__ import annotations

import math
from collections.abc import Iterable
from pathlib import Path

import cobra
import pandas as pd


def load_delta_g_csv(
    entities: Iterable,
    path: str | Path,
    *,
    id_column: str = "Var1",
    value_column: str = "Var2",
    note_key: str = "deltaG",
    verbose: bool = False,
) -> int:
    """Record ``note_key`` on each entity from a CSV of ``id → value``.

    Every matched value is recorded exactly as it appears in the CSV,
    including yeast-GEM's own "no measurement" placeholder (``10000000.0``,
    on 777 of its 4102 reaction rows) -- this function does not interpret
    CSV values, matching RAVEN's own ``loadDeltaGCSV``. Callers that need to
    treat a particular value as absent should filter it themselves after
    loading.

    Parameters
    ----------
    entities
        Cobra entities (``model.metabolites`` or ``model.reactions``).
    path
        Path to the CSV.
    id_column, value_column
        Column labels in the CSV. Defaults match the yeast-GEM
        convention (``Var1`` / ``Var2`` from MATLAB's ``array2table``).
    note_key
        Key under which the value is stored on ``entity.notes``.
        Default ``"deltaG"``.
    verbose
        Print a summary of unmatched entity ids.

    Returns
    -------
    The number of entities that were recorded (i.e. matched the CSV and
    carried a value).
    """
    df = pd.read_csv(path)
    if id_column not in df.columns or value_column not in df.columns:
        raise ValueError(
            f"{path} columns {list(df.columns)} do not contain "
            f"both {id_column!r} and {value_column!r}"
        )
    lookup = dict(zip(df[id_column], df[value_column], strict=True))

    stamped = 0
    missing: list[str] = []
    for entity in entities:
        value = lookup.get(entity.id)
        if value is None or (isinstance(value, float) and math.isnan(value)):
            missing.append(entity.id)
            continue
        entity.notes[note_key] = str(value)
        stamped += 1

    if verbose and missing:
        print(
            f"{len(missing)} entity ids were not in {path} "
            f"(e.g. {missing[:3]}); left untouched."
        )
    return stamped


def save_delta_g_csv(
    entities: Iterable,
    path: str | Path,
    *,
    id_column: str = "Var1",
    value_column: str = "Var2",
    note_key: str = "deltaG",
) -> int:
    """Dump ``entity.notes[note_key]`` for each entity to a CSV.

    Entities without ``note_key`` set get ``NaN`` written, preserving
    one-row-per-entity ordering (mirrors MATLAB's ``array2table([ids,
    values])`` behaviour).

    Returns
    -------
    The number of rows written (= number of entities).
    """
    rows: list[tuple[str, float]] = []
    for entity in entities:
        raw = entity.notes.get(note_key)
        if raw is None:
            value: float = math.nan
        else:
            try:
                value = float(raw)
            except ValueError:
                value = math.nan
        rows.append((entity.id, value))
    pd.DataFrame(rows, columns=[id_column, value_column]).to_csv(path, index=False)
    return len(rows)


__all__ = ["load_delta_g_csv", "save_delta_g_csv"]
_ = cobra  # silence "imported but unused" — kept for type-checker/IDE context
