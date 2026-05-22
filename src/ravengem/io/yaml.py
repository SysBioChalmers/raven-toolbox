"""Read and write RAVEN/Metabolic-Atlas YAML models.

Port of RAVEN ``readYAMLmodel.m`` / ``writeYAMLmodel.m``.

The lens here is important. RAVEN's YAML is deliberately the Human-GEM / Metabolic
Atlas schema, which is **cobra's own YAML format** (``!!omap`` ordered maps,
loaded/written by ``cobra.io.load_yaml_model`` / ``save_yaml_model``). cobra
already reads a real yeast-GEM/Human-GEM file end to end. So this module does
*not* re-implement YAML parsing or the standard reaction/metabolite/gene fields.

It adds only the two things cobra silently drops:

1. **Model identity & provenance.** RAVEN keeps ``id``, ``name``, ``version`` and
   provenance (``defaultLB``/``defaultUB``, authors, ``taxonomy``, ...) in a
   top-level ``metaData:`` block. cobra ignores it, so a cobra-loaded RAVEN model
   has ``id is None``. We map ``metaData.id``/``name`` onto the model and keep the
   whole block in ``model.notes['metaData']`` so it round-trips.
2. **RAVEN-only per-entry fields** that sit as top-level keys on each entry and
   have no slot in cobra's dict: ``smiles``/``inchis``/``deltaG``/``metFrom`` on
   metabolites, ``confidence_score``/``references``/``rxnFrom``/``deltaG`` on
   reactions, ``protein`` on genes. We move them into the object's ``notes`` on
   read and lift them back to top-level on write.

Any other unrecognised top-level sections (e.g. GECKO's ``ec-rxns`` /
``ec-enzymes`` / ``gecko_light``) are preserved verbatim in
``model.notes['_yaml_sections']`` so a read/write round-trip never destroys an
enzyme-constrained model, even though ravengem does not interpret those sections.
"""
from __future__ import annotations

from pathlib import Path
from typing import Union

import cobra
from cobra.io.dict import model_from_dict, model_to_dict
from cobra.io.yaml import yaml as _cobra_yaml  # ruamel round-trip YAML (handles !!omap)

# RAVEN-only fields written as top-level keys on each entry, with no place in
# cobra's metabolite/reaction/gene dict. Moved to/from the object's notes.
_MET_EXTRAS = ("smiles", "inchis", "deltaG", "metFrom")
_RXN_EXTRAS = ("confidence_score", "references", "rxnFrom", "deltaG")
_GENE_EXTRAS = ("protein",)

# cobra's own top-level dict keys (everything else is a foreign section).
_COBRA_TOP_KEYS = frozenset(
    {"metabolites", "reactions", "genes", "compartments", "id", "name", "version"}
)


def _to_plain(obj):
    """Recursively convert ruamel containers/scalars to plain Python types."""
    if isinstance(obj, dict):
        return {str(k): _to_plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_plain(v) for v in obj]
    if isinstance(obj, bool):
        return obj
    if isinstance(obj, int):
        return int(obj)
    if isinstance(obj, float):
        return float(obj)
    return obj if obj is None else str(obj) if not isinstance(obj, str) else obj


def _absorb_entry_extras(entries: list[dict], extras: tuple[str, ...]) -> None:
    """Move RAVEN-only top-level keys into each entry's ``notes`` dict (in place)."""
    for entry in entries:
        notes = dict(entry.get("notes") or {})
        for key in extras:
            if key in entry:
                notes[key] = entry.pop(key)
        if notes:
            entry["notes"] = notes


def _emit_entry_extras(entries: list[dict], extras: tuple[str, ...]) -> None:
    """Lift RAVEN-only keys from each entry's ``notes`` back to top level (in place)."""
    for entry in entries:
        notes = entry.get("notes")
        if not isinstance(notes, dict):
            continue
        for key in extras:
            if key in notes:
                entry[key] = notes.pop(key)
        if not notes:
            entry.pop("notes", None)


def read_yaml_model(path: Union[str, Path]) -> "cobra.Model":
    """Read a RAVEN/Metabolic-Atlas YAML model into a ``cobra.Model``.

    Port of RAVEN ``readYAMLmodel.m``. Standard fields are parsed by cobra;
    RAVEN's ``metaData`` block and RAVEN-only per-entry fields are preserved
    (see module docstring).
    """
    with open(path, encoding="utf-8") as handle:
        raw = _to_plain(_cobra_yaml.load(handle))

    if not isinstance(raw, dict):
        raise ValueError(
            f"{path}: top-level YAML is a {type(raw).__name__}, not a mapping."
        )

    metadata = raw.pop("metaData", None) or {}
    foreign = {k: raw.pop(k) for k in list(raw) if k not in _COBRA_TOP_KEYS}

    _absorb_entry_extras(raw.get("metabolites", []), _MET_EXTRAS)
    _absorb_entry_extras(raw.get("reactions", []), _RXN_EXTRAS)
    _absorb_entry_extras(raw.get("genes", []), _GENE_EXTRAS)

    model = model_from_dict(raw)

    if metadata:
        if metadata.get("id"):
            model.id = metadata["id"]
        if metadata.get("name"):
            model.name = metadata["name"]
        model.notes["metaData"] = metadata
    if foreign:
        model.notes["_yaml_sections"] = foreign

    return model


def write_yaml_model(model: "cobra.Model", path: Union[str, Path]) -> None:
    """Write a ``cobra.Model`` to RAVEN/Metabolic-Atlas YAML.

    Port of RAVEN ``writeYAMLmodel.m``. Emits a ``metaData`` block (model
    identity + any stored provenance), lifts RAVEN-only fields back to per-entry
    top-level keys, and re-emits any preserved foreign sections (e.g. GECKO ec).
    """
    notes = dict(model.notes or {})
    stored_meta = notes.pop("metaData", None) or {}
    foreign = notes.pop("_yaml_sections", None) or {}

    doc = model_to_dict(model)
    # Drop the housekeeping notes we manage ourselves so they don't leak into the
    # serialized model-level notes (cobra's dict has no model-level notes key, but
    # be safe in case of future cobra versions).

    _emit_entry_extras(doc.get("metabolites", []), _MET_EXTRAS)
    _emit_entry_extras(doc.get("reactions", []), _RXN_EXTRAS)
    _emit_entry_extras(doc.get("genes", []), _GENE_EXTRAS)

    metadata = dict(stored_meta)
    if model.id:
        metadata["id"] = model.id
    if model.name:
        metadata["name"] = model.name

    out: dict = {}
    if metadata:
        out["metaData"] = metadata
    out.update(_to_plain(doc))
    for key, value in foreign.items():
        out[key] = value

    with open(path, "w", encoding="utf-8") as handle:
        _cobra_yaml.dump(out, handle)
