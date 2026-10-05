"""Model version, date and descriptive metadata.

RAVEN keeps ``version`` and ``date`` on the model itself (``model.version``,
``model.date``). A cobra model has no such attributes, so raven-toolbox keeps them in
``model.notes``: ``notes["version"]`` and ``notes["metaData"]`` (which holds ``date`` and
the descriptive fields). :func:`write_yaml_model` reads the version from the first and the
date from the second, so setting them by hand means touching both. These two functions are
the one place that does it.
"""
from __future__ import annotations

import datetime

import cobra

from raven_toolbox.io.yaml import _META_ANNOTATION_FIELDS

__all__ = ["get_model_metadata", "set_model_metadata"]


def get_model_metadata(model: cobra.Model) -> dict:
    """The model's metadata as a plain dict: ``version``, ``date`` and any descriptive fields.

    ``version`` is taken from ``notes["version"]``, falling back to ``notes["metaData"]``.
    Fields the model does not carry are absent.
    """
    notes = model.notes or {}
    meta = dict(notes.get("metaData") or {})
    version = notes.get("version", meta.get("version"))
    if version is not None:
        meta["version"] = version
    return meta


def set_model_metadata(
    model: cobra.Model,
    *,
    version: str | None = None,
    date: str | datetime.date | None = None,
    **fields: str,
) -> None:
    """Set the version, date and descriptive fields that the YAML export writes, in place.

    Parameters
    ----------
    version
        Model version, e.g. ``"1.2.0"``.
    date
        Model date, as ``YYYY-MM-DD`` or a ``datetime.date``. The YAML writer stamps today's
        date only when the model has none, so pass one when the date should change.
    **fields
        Descriptive ``metaData`` fields: ``givenName``, ``familyName``, ``authors``,
        ``email``, ``organization``, ``taxonomy``, ``note`` and ``sourceUrl``.

    Arguments left out are not touched.
    """
    unknown = sorted(set(fields) - set(_META_ANNOTATION_FIELDS))
    if unknown:
        raise ValueError(f"unknown metadata field(s) {unknown}; allowed: {list(_META_ANNOTATION_FIELDS)}")
    notes = dict(model.notes or {})
    meta = dict(notes.get("metaData") or {})
    if version is not None:
        notes["version"] = meta["version"] = str(version)
    if date is not None:
        meta["date"] = date.isoformat() if isinstance(date, datetime.date) else str(date)
    meta.update(fields)
    notes["metaData"] = meta
    model.notes = notes
