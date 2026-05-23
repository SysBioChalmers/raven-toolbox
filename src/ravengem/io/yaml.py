"""Read and write RAVEN/cobrapy YAML models.

Aligned to RAVEN ``writeYAMLmodel.m`` / ``readYAMLmodel.m`` as of the
``feat/geckopy-compat-yaml`` work (commit fa281a1), whose writer emits **cobra's
native ``!!omap`` YAML**. Because the format *is* cobra's, the standard model
content — id, name, compartments, and per-entry id/name/compartment/formula/
charge/bounds/gene_reaction_rule/objective_coefficient/subsystem/metabolites and
the whole ``annotation`` block (which carries ``smiles`` for metabolites,
``ec-code`` for reactions, and all MIRIAM cross-references) — is read and written
by ``cobra.io`` directly.

This module only handles what cobra drops or mishandles:

* **RAVEN-only top-level per-entry keys** that cobra ignores: ``inchis``,
  ``deltaG``, ``metFrom`` and the free-text ``notes`` (metNotes) on metabolites;
  ``confidence_score``, ``references``, ``rxnFrom``, ``deltaG`` and ``notes``
  (rxnNotes) on reactions; ``protein`` on genes. These are stashed in the cobra
  object's ``.notes`` dict on read and lifted back to top-level keys on write.
* **Model-level extras** cobra ignores: ``version``, the ``metaData`` provenance
  block, and the GECKO sections (``gecko_light``/``ec-rxns``/``ec-enzymes``),
  preserved on ``model.notes`` for round-tripping.

The reader also accepts the older RAVEN files (id/name nested in ``metaData``).
"""
from __future__ import annotations

from collections import OrderedDict
from pathlib import Path

import cobra
from cobra.io.dict import model_from_dict, model_to_dict
from cobra.io.yaml import yaml as _cobra_yaml  # ruamel round-trip YAML (handles !!omap)

# RAVEN-only top-level per-entry keys -> the key used inside the cobra object's
# .notes dict. ('notes' is RAVEN's free-text metNotes/rxnNotes; stored under
# 'note' to avoid colliding with the notes container itself.)
_MET_FIELDS = (("inchis", "inchis"), ("deltaG", "deltaG"), ("metFrom", "metFrom"), ("notes", "note"))
_RXN_FIELDS = (
    ("confidence_score", "confidence_score"),
    ("references", "references"),
    ("rxnFrom", "rxnFrom"),
    ("deltaG", "deltaG"),
    ("notes", "note"),
)
_GENE_FIELDS = (("protein", "protein"),)

_COBRA_TOP_KEYS = frozenset({"metabolites", "reactions", "genes", "compartments", "id", "name"})


def _to_plain(obj):
    if isinstance(obj, dict):
        return {str(k): _to_plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_plain(v) for v in obj]
    if isinstance(obj, bool) or obj is None:
        return obj
    if isinstance(obj, int):
        return int(obj)
    if isinstance(obj, float):
        return float(obj)
    return obj if isinstance(obj, str) else str(obj)


def _capture_entry_fields(entries, fields):
    """Pop RAVEN-only top-level keys off each entry into a parallel notes dict.

    Returns a list of ``{notes_key: value}`` dicts aligned with ``entries`` (so
    cobra never sees these keys), to be attached to the built objects afterwards.
    """
    captured = []
    for entry in entries:
        notes = {}
        for yaml_key, notes_key in fields:
            if yaml_key in entry:
                notes[notes_key] = entry.pop(yaml_key)
        captured.append(notes)
    return captured


def read_yaml_model(path: str | Path) -> cobra.Model:
    """Read a RAVEN/cobrapy YAML model into a ``cobra.Model``."""
    with open(path, encoding="utf-8") as handle:
        raw = _to_plain(_cobra_yaml.load(handle))

    if not isinstance(raw, dict):
        raise ValueError(f"{path}: top-level YAML is a {type(raw).__name__}, not a mapping.")

    metadata = raw.pop("metaData", None) or {}
    version = raw.pop("version", None)
    foreign = {k: raw.pop(k) for k in list(raw) if k not in _COBRA_TOP_KEYS}

    met_notes = _capture_entry_fields(raw.get("metabolites", []), _MET_FIELDS)
    rxn_notes = _capture_entry_fields(raw.get("reactions", []), _RXN_FIELDS)
    gene_notes = _capture_entry_fields(raw.get("genes", []), _GENE_FIELDS)

    model = model_from_dict(raw)

    for met, notes in zip(model.metabolites, met_notes, strict=False):
        met.notes = notes
    for rxn, notes in zip(model.reactions, rxn_notes, strict=False):
        rxn.notes = notes
    for gene, notes in zip(model.genes, gene_notes, strict=False):
        gene.notes = notes

    # Legacy files keep id/name inside metaData; restore them if cobra found none.
    if metadata.get("id") and not model.id:
        model.id = metadata["id"]
    if metadata.get("name") and not model.name:
        model.name = metadata["name"]
    if metadata:
        model.notes["metaData"] = metadata
    if version is not None:
        model.notes["version"] = version
    if foreign:
        model.notes["_yaml_sections"] = foreign

    return model


def _emit_entry_fields(entries, fields):
    """Lift RAVEN-only keys out of each entry's ``notes`` dict to top level."""
    for entry in entries:
        notes = entry.pop("notes", None)
        if not isinstance(notes, dict):
            continue
        notes = dict(notes)
        for yaml_key, notes_key in fields:
            if notes_key in notes:
                entry[yaml_key] = notes.pop(notes_key)
        if notes:  # any non-RAVEN notes survive as a notes block (unless 'notes' taken)
            entry.setdefault("notes", notes)


def write_yaml_model(
    model: cobra.Model, path: str | Path, *, sort_ids: bool = False
) -> None:
    """Write a ``cobra.Model`` to RAVEN/cobrapy (``!!omap``) YAML.

    With ``sort_ids=True`` metabolites/reactions/genes/compartments are written
    in alphabetical order (diff-friendly), without modifying ``model``.
    """
    model_notes = dict(model.notes or {})
    stored_meta = model_notes.pop("metaData", None) or {}
    version = model_notes.pop("version", None)
    foreign = model_notes.pop("_yaml_sections", None) or {}

    doc = OrderedDict(_to_plain(model_to_dict(model)))

    if sort_ids:
        for section in ("metabolites", "reactions", "genes"):
            if section in doc:
                doc[section] = sorted(doc[section], key=lambda e: e.get("id", ""))
        if isinstance(doc.get("compartments"), dict):
            doc["compartments"] = dict(sorted(doc["compartments"].items()))

    _emit_entry_fields(doc.get("metabolites", []), _MET_FIELDS)
    _emit_entry_fields(doc.get("reactions", []), _RXN_FIELDS)
    _emit_entry_fields(doc.get("genes", []), _GENE_FIELDS)

    # cobra dict order is metabolites, reactions, genes, id, name, compartments;
    # append version / gecko_light / metaData / ec-* like RAVEN's writer.
    if version is not None:
        doc["version"] = version
    metadata = dict(stored_meta)
    if model.id:
        metadata.setdefault("id", model.id)
    if model.name:
        metadata.setdefault("name", model.name)
    for key in ("gecko_light",):
        if key in foreign:
            doc[key] = foreign.pop(key)
    if metadata:
        doc["metaData"] = metadata
    for key, value in foreign.items():
        doc[key] = value

    with open(path, "w", encoding="utf-8") as handle:
        _cobra_yaml.dump(doc, handle)
