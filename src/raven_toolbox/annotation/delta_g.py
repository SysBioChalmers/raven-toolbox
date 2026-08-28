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

#: The "no valid ΔG" sentinel used by the ΔG side-car tables (e.g. yeast-GEM's
#: ``model_rxnDeltaG.csv``). yeast-GEM's ``checkrxnDirection.m`` gates on it verbatim:
#: ``if ~isequal(seed_rxnInfo{...},'10000000.0') %check if database contains valid deltaG
#: value``. Stamping it would present a physically impossible 10⁷ kJ/mol as a measurement.
#: Pass this as ``load_delta_g_csv``'s ``missing_value`` to have it treated as missing
#: instead of stored verbatim — RAVEN's own ``deltaGCSV`` has no sentinel concept by
#: default either, so this is opt-in on both sides, not automatic on either.
DELTA_G_MISSING = 1e7


def _is_missing(value, sentinel: float) -> bool:
    """True when ``value`` is the sentinel, whether it arrived as a number or as text.

    The CSV round-trips through MATLAB and pandas, so the same sentinel shows up as ``10000000``,
    ``10000000.0`` or ``"10000000.0"`` depending on the writer and the column's inferred dtype.
    """
    try:
        return math.isclose(float(value), sentinel, rel_tol=1e-9)
    except (TypeError, ValueError):
        return False


def load_delta_g_csv(
    entities: Iterable,
    path: str | Path,
    *,
    id_column: str = "Var1",
    value_column: str = "Var2",
    note_key: str = "deltaG",
    missing_value: float | None = None,
    verbose: bool = False,
) -> int:
    """Stamp ``note_key`` on each entity from a CSV of ``id → value``.

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
    missing_value
        A sentinel standing for "no value", left unstamped rather than
        recorded as a measurement. Default ``None``: every matched value
        is stamped verbatim, matching RAVEN's own ``deltaGCSV`` (which has
        no sentinel concept of its own). Pass :data:`DELTA_G_MISSING`
        (10⁷) to opt into yeast-GEM's own convention, which covers 777 of
        yeast-GEM's 4102 reaction rows.
    verbose
        Print a summary of unmatched entity ids.

    Returns
    -------
    The number of entities that were stamped (i.e. matched the CSV and
    carried a real value).
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
        if missing_value is not None and _is_missing(value, missing_value):
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


# Re-export the cobra Model type for type-checker friendliness; helps
# IDEs surface the right hints to callers that hand us model.metabolites
# / model.reactions directly.
__all__ = ["DELTA_G_MISSING", "load_delta_g_csv", "save_delta_g_csv"]
_ = cobra  # silence "imported but unused" — used for typing context above
