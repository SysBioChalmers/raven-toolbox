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
* **Model-level extras** cobra ignores: ``version`` and the ``metaData``
  provenance block, preserved on ``model.notes`` for round-tripping.
* **The GECKO ec sections** (``ec-rxns``/``ec-enzymes``/``gecko_light``): on
  read, parsed into a typed :class:`~raven_python.io.ec_data.EcData` and
  attached as ``model.ec``; on write, serialised back to top-level sections
  whenever ``model.ec`` is present. This mirrors RAVEN ``readYAMLmodel.m`` /
  ``writeYAMLmodel.m``, which populate the ``model.ec`` struct when the YAML
  defines it. Any other unknown top-level keys are still preserved opaquely
  via ``model.notes['_yaml_sections']``.

Legacy quirks the reader also accepts (silent normalisation):

* older RAVEN files with ``id`` / ``name`` nested in ``metaData``;
* per-metabolite top-level ``smiles`` (lifted into ``annotation['smiles']``);
* per-reaction top-level ``eccodes`` (lifted into ``annotation['ec-code']`` —
  the cobra-standard place where geckopy reads EC numbers);
* very old RAVEN files written as a bare ``-`` sequence of single-key mappings
  rather than one big mapping;
* MATLAB GECKO ecModels whose ``usage_prot_*`` and ``prot_pool_exchange``
  reactions use the older reverse-direction convention (negative lower bound,
  swapped stoichiometry signs); these are flipped on load so consumers always
  see the forward convention.
"""
from __future__ import annotations

import gzip
import warnings
from collections import OrderedDict
from pathlib import Path

import cobra
from cobra.io.dict import model_from_dict, model_to_dict
from cobra.io.yaml import yaml as _cobra_yaml  # ruamel round-trip YAML (handles !!omap)

from raven_python.io.ec_data import (
    EcData,
    ec_data_from_yaml_sections,
    ec_data_to_yaml_sections,
)


def _open_text(path: str | Path, mode: str):
    """Open ``path`` as a text handle, transparently gzipping when it ends ``.gz``."""
    if str(path).endswith(".gz"):
        return gzip.open(path, f"{mode}t", encoding="utf-8")
    return open(path, mode, encoding="utf-8")

# RAVEN-only top-level per-entry keys -> the key used inside the cobra object's
# .notes dict. ('notes' is RAVEN's free-text metNotes/rxnNotes; stored under
# 'note' to avoid colliding with the notes container itself.)
_MET_FIELDS = (("inchis", "inchis"), ("deltaG", "deltaG"), ("metFrom", "metFrom"), ("notes", "note"))
_RXN_FIELDS = (
    ("references", "references"),
    ("rxnFrom", "rxnFrom"),
    ("deltaG", "deltaG"),
    ("confidence_score", "confidence_score"),
    ("notes", "note"),
)
# Legacy YAML keys accepted on READ for reaction notes. Old RAVEN MATLAB writers
# used "rxnNotes"; the canonical key (matching cobrapy and the current MATLAB
# writer) is "notes". When both appear, "notes" wins.
_LEGACY_RXN_KEY_ALIASES = (("rxnNotes", "notes"),)
_GENE_FIELDS = (("protein", "protein"),)

_COBRA_TOP_KEYS = frozenset({"metabolites", "reactions", "genes", "compartments", "id", "name"})

# Top-level keys consumed by the typed EcData layer; not re-emitted as opaque
# `_yaml_sections` on round-trip (the source of truth is ``model.ec``).
_EC_TOP_KEYS = frozenset({"ec-rxns", "ec-enzymes", "gecko_light"})


def _to_plain(obj):
    """Coerce ruamel/numpy scalars to plain Python primitives.

    Dicts are returned as ``OrderedDict`` so that the round-trip dumper
    emits them with the ``!!omap`` tag (ruamel's CommentedMap subclass
    is replaced by a plain ``OrderedDict`` to avoid carrying ruamel's
    own type tags through). Returning a plain ``dict`` instead would
    drop the ``!!omap`` tag and produce files that RAVEN's MATLAB
    reader (a line-based parser keyed on ``!!omap``) cannot load.
    """
    if isinstance(obj, dict):
        return OrderedDict((str(k), _to_plain(v)) for k, v in obj.items())
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
    """Read a RAVEN/cobrapy YAML model into a ``cobra.Model``.

    Convenience wrapper around :func:`model_from_yaml_data` that opens the
    file (transparently un-gzipping ``.gz``) and parses the YAML. Callers
    that need to pre-process the document (e.g. lift further non-standard
    fields that cobra doesn't recognise) can read+normalise themselves and
    call :func:`model_from_yaml_data` with the resulting dict.

    Accepts both the cobra `!!omap` shape and a very old RAVEN shape where
    the document root is a bare ``-`` sequence of single-key mappings; the
    latter is merged into one mapping before parsing.
    """
    with _open_text(path, "r") as handle:
        raw = _to_plain(_cobra_yaml.load(handle))

    if isinstance(raw, list):
        # Very old RAVEN files: a sequence of one-key mappings instead of
        # one big !!omap. Merge into a single dict before parsing.
        merged: dict = {}
        for item in raw:
            if isinstance(item, dict):
                merged.update(item)
        raw = merged

    if not isinstance(raw, dict):
        raise ValueError(
            f"{path}: top-level YAML must be a mapping or a sequence of "
            f"single-key mappings, got {type(raw).__name__}."
        )
    return model_from_yaml_data(raw)


def model_from_yaml_data(raw: dict) -> cobra.Model:
    """Build a ``cobra.Model`` from an already-parsed RAVEN/cobrapy YAML dict.

    Performs three jobs in order:

    1. **cobra-shaped portion:** strips and restores RAVEN-only per-entry
       side-fields onto each entry's ``.notes``; lifts ``id`` / ``name``
       out of legacy ``metaData``; preserves ``version`` and ``metaData``
       on ``model.notes`` for round-trip.
    2. **legacy quirks:** lifts per-metabolite top-level ``smiles`` into
       ``annotation['smiles']`` and per-reaction top-level ``eccodes`` into
       ``annotation['ec-code']`` (older RAVEN/MATLAB GECKO files emitted
       these at the top level); flips the older reverse-direction
       ``usage_prot_*`` / ``prot_pool_exchange`` convention to the
       forward convention.
    3. **GECKO ec sections:** when ``ec-rxns`` / ``ec-enzymes`` are
       present, parses them into a typed :class:`EcData` and attaches
       it as ``model.ec``. Other unknown top-level keys land opaquely on
       ``model.notes['_yaml_sections']`` for round-trip.

    ``raw`` is mutated in place — copy it first if the caller needs the
    original.
    """
    metadata = raw.pop("metaData", None) or {}
    version = raw.pop("version", None)
    foreign = {k: raw.pop(k) for k in list(raw) if k not in _COBRA_TOP_KEYS}

    # Legacy quirk: per-metabolite top-level `smiles` -> annotation.smiles.
    # Done before model_from_dict so cobra sees the annotation in its
    # canonical place. No-op on current files.
    _lift_smiles_to_annotation(raw.get("metabolites"))

    # Legacy quirk: per-reaction top-level `eccodes` -> annotation['ec-code'].
    # EC numbers are standard cobra annotation; older RAVEN/MATLAB files put
    # them at the reaction top level, which hid them from cobra/geckopy (which
    # read annotation['ec-code']). Lift before model_from_dict. No-op on
    # current cobra-shaped files.
    _lift_eccodes_to_annotation(raw.get("reactions"))

    # Normalise legacy reaction-side YAML keys (e.g. RAVEN MATLAB's
    # ``rxnNotes`` -> the canonical ``notes``) before any field capture so
    # the capture step sees a single key per concept.
    for entry in raw.get("reactions", []):
        if not isinstance(entry, dict):
            continue
        for legacy_key, canonical_key in _LEGACY_RXN_KEY_ALIASES:
            if legacy_key in entry and canonical_key not in entry:
                entry[canonical_key] = entry.pop(legacy_key)
            elif legacy_key in entry:
                # Canonical wins; just drop the legacy duplicate.
                entry.pop(legacy_key)

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

    # Legacy quirk: flip reverse-direction usage_prot_* / prot_pool_exchange.
    # No-op when none match.
    _flip_legacy_prot_direction(model)

    # RAVEN convention keeps id/name/version inside metaData; lift them
    # onto the model so cobra-shaped accessors find them too. Older
    # cobra-style files (id/name/version at root) are handled by the
    # cobra reader; the metaData lookups below are pure fallbacks.
    if metadata.get("id") and not model.id:
        model.id = metadata["id"]
    if metadata.get("name") and not model.name:
        model.name = metadata["name"]
    if version is None and metadata.get("version") is not None:
        version = metadata["version"]
    if metadata:
        model.notes["metaData"] = metadata
    if version is not None:
        model.notes["version"] = version

    # Pop the ec sections out of `foreign` and into a typed EcData.
    # The remaining unknown keys round-trip opaquely. Pre-shim RAVEN
    # MATLAB writes wrote `geckoLight: "true"` inside metaData (rather
    # than the current top-level `gecko_light`); honour the legacy
    # placement too — keep the metaData entry untouched (round-trip)
    # and surface it at the top level so EcData picks it up.
    legacy_gecko = metadata.get("geckoLight")
    if legacy_gecko is not None and "gecko_light" not in foreign:
        if isinstance(legacy_gecko, str):
            foreign["gecko_light"] = legacy_gecko.lower() == "true"
        else:
            foreign["gecko_light"] = bool(legacy_gecko)
    ec_sections = {k: foreign.pop(k) for k in list(foreign) if k in _EC_TOP_KEYS}
    ec_data = ec_data_from_yaml_sections(ec_sections)
    if ec_data is not None:
        model.ec = ec_data
    if foreign:
        model.notes["_yaml_sections"] = foreign

    return model


# --------------------------------------------------------------------------- #
# Legacy quirk normalisers
# --------------------------------------------------------------------------- #

def _lift_smiles_to_annotation(metabolites) -> None:
    """Move per-metabolite top-level ``smiles`` into ``annotation['smiles']``.

    Older MATLAB GECKO ecModel writers placed SMILES at the metabolite top
    level; the cobra/raven convention is to nest them inside ``annotation``.
    Normalises in place; no-op when no metabolite carries a top-level
    ``smiles`` key.
    """
    if not isinstance(metabolites, list):
        return
    for met in metabolites:
        if not (isinstance(met, dict) and "smiles" in met):
            continue
        smiles = met.pop("smiles")
        annotation = met.get("annotation")
        if not isinstance(annotation, dict):
            annotation = {}
            met["annotation"] = annotation
        if "smiles" not in annotation and smiles:
            annotation["smiles"] = (
                smiles if isinstance(smiles, list) else [smiles]
            )


def _lift_eccodes_to_annotation(reactions) -> None:
    """Move a reaction's legacy top-level ``eccodes`` into ``annotation['ec-code']``.

    EC numbers are a standard MIRIAM cross-reference, so the cobra/raven
    convention is to carry them inside ``annotation`` under the ``ec-code``
    key — where cobra and geckopy read them — not as a RAVEN-only top-level
    field. Older RAVEN/MATLAB writers emitted a top-level ``eccodes`` (a
    ``;``-joined string or a list of codes); lift it into the canonical
    place. Normalises in place; no-op when no reaction carries a top-level
    ``eccodes`` key. A native ``annotation['ec-code']`` (if already present)
    wins and is left untouched.
    """
    if not isinstance(reactions, list):
        return
    for rxn in reactions:
        if not (isinstance(rxn, dict) and "eccodes" in rxn):
            continue
        codes = _eccodes_to_list(rxn.pop("eccodes"))
        if not codes:
            continue
        annotation = rxn.get("annotation")
        if not isinstance(annotation, dict):
            annotation = {}
            rxn["annotation"] = annotation
        annotation.setdefault("ec-code", codes)


def _eccodes_to_list(value) -> list:
    """Normalise a RAVEN ``eccodes`` value to a list of trimmed code strings.

    Accepts a ``;``-joined string (RAVEN MATLAB's ``getECstring`` form) or an
    already-split list; drops empty tokens.
    """
    if value is None:
        return []
    items = value if isinstance(value, list) else str(value).split(";")
    return [str(s).strip() for s in items if str(s).strip()]


def _flip_legacy_prot_direction(model: cobra.Model) -> None:
    """Flip pre-forward-direction protein reactions in place.

    Older MATLAB GECKO ecModels defined ``usage_prot_*`` and
    ``prot_pool_exchange`` as "reverse" reactions: their flux was
    negative, and the stoichiometry signs were correspondingly swapped.
    The current convention (in both geckopy and recent MATLAB GECKO)
    treats them as ordinary forward reactions with positive flux. When a
    loaded model still uses the older convention we flip the affected
    reactions in place so consumers never have to handle two shapes.

    The signature we look for is any ``usage_prot_*`` or
    ``prot_pool_exchange`` reaction whose lower bound is negative.
    """
    flipped: list[str] = []
    for rxn in model.reactions:
        if not (
            rxn.id.startswith("usage_prot_")
            or rxn.id == "prot_pool_exchange"
        ):
            continue
        if rxn.lower_bound >= -1e-9:
            continue
        rxn.add_metabolites(
            {met: -2.0 * coef for met, coef in rxn.metabolites.items()},
            combine=True,
        )
        lb, ub = rxn.lower_bound, rxn.upper_bound
        rxn.lower_bound = -ub
        rxn.upper_bound = -lb
        flipped.append(rxn.id)
    if flipped:
        warnings.warn(
            f"ecModel uses the older reverse-direction convention for "
            f"{len(flipped)} protein usage/pool reaction(s); flipping to "
            "the current forward convention.",
            stacklevel=3,
        )


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
        # Preserve any remaining (non-RAVEN) notes. The RAVEN free-text note is lifted
        # to the YAML key "notes"; if leftovers also exist, merge them with it under
        # that key (rather than silently dropping the leftovers).
        if notes:
            if "notes" in entry:
                notes["note"] = entry["notes"]
            entry["notes"] = notes


def write_yaml_model(
    model: cobra.Model, path: str | Path, *, sort_ids: bool = False
) -> None:
    """Write a ``cobra.Model`` to RAVEN/cobrapy (``!!omap``) YAML.

    When ``model.ec`` is a populated :class:`EcData`, the ``gecko_light``
    flag and the ``ec-rxns`` / ``ec-enzymes`` top-level sections are
    emitted from it (numpy/ruamel scalars are coerced to plain Python
    primitives en route, so the dumper never sees them).

    With ``sort_ids=True`` metabolites/reactions/genes/compartments are
    written in alphabetical order (diff-friendly), without modifying
    ``model``.
    """
    model_notes = dict(model.notes or {})
    stored_meta = model_notes.pop("metaData", None) or {}
    version = model_notes.pop("version", None)
    foreign = model_notes.pop("_yaml_sections", None) or {}

    doc = OrderedDict(_to_plain(model_to_dict(model)))

    # cobra's model_to_dict serialises model.notes verbatim into doc["notes"],
    # so the three management keys we just lifted out would otherwise also
    # appear nested inside the notes section. Strip them; preserve any other
    # genuine notes the caller stored on the model.
    doc_notes = doc.get("notes")
    if isinstance(doc_notes, dict):
        for key in ("metaData", "version", "_yaml_sections"):
            doc_notes.pop(key, None)
        if not doc_notes:
            doc.pop("notes", None)

    if sort_ids:
        for section in ("metabolites", "reactions", "genes"):
            if section in doc:
                doc[section] = sorted(doc[section], key=lambda e: e.get("id", ""))
        if isinstance(doc.get("compartments"), dict):
            doc["compartments"] = dict(sorted(doc["compartments"].items()))

    _emit_entry_fields(doc.get("metabolites", []), _MET_FIELDS)
    _emit_entry_fields(doc.get("reactions", []), _RXN_FIELDS)
    _emit_entry_fields(doc.get("genes", []), _GENE_FIELDS)

    # cobra's _gene_to_dict always emits `name: ''` because name is a
    # required attribute; RAVEN MATLAB skips empty names. Drop the
    # empty-string entry so the two writers produce identical genes.
    for gene in doc.get("genes", []) or ():
        if isinstance(gene, dict) and gene.get("name") == "":
            gene.pop("name", None)

    # ec sections come from the typed model.ec (when present), not from the
    # opaque foreign-keys stash. Drop any stale ec-* entries in `foreign` so
    # they can't conflict with the EcData-derived ones.
    for ec_key in _EC_TOP_KEYS:
        foreign.pop(ec_key, None)
    ec = getattr(model, "ec", None)
    ec_sections = (
        _to_plain(ec_data_to_yaml_sections(ec))
        if isinstance(ec, EcData) and (ec.n_rxns or ec.n_enzymes)
        else None
    )

    # Final document order (matches RAVEN MATLAB writeYAMLmodel):
    #   metaData, metabolites, reactions, genes, compartments,
    #   gecko_light, ec-rxns, ec-enzymes, <opaque foreign>
    # id/name/version live inside metaData (the RAVEN convention) — they
    # are NOT emitted at root level. Cobra reading such a file recovers
    # the lists and compartments and leaves model.id empty (consistent
    # with how RAVEN MATLAB has always laid these files out).
    metadata = OrderedDict(stored_meta)
    if model.id:
        metadata.setdefault("id", model.id)
    if model.name:
        metadata.setdefault("name", model.name)
    if version is not None:
        metadata.setdefault("version", version)

    # cobra's model_to_dict put id / name at root level; drop them so they
    # don't duplicate the metaData copy.
    doc.pop("id", None)
    doc.pop("name", None)

    ordered = OrderedDict()
    if metadata:
        ordered["metaData"] = metadata
    for key in ("metabolites", "reactions", "genes", "compartments", "notes"):
        if key in doc:
            ordered[key] = doc.pop(key)
    if ec_sections is not None:
        ordered["gecko_light"] = ec_sections["gecko_light"]
        ordered["ec-rxns"] = ec_sections["ec-rxns"]
        ordered["ec-enzymes"] = ec_sections["ec-enzymes"]
    # Carry over any cobra-shaped fields we didn't classify above
    # (defensive: keeps forward compatibility with future cobra additions).
    for key, value in doc.items():
        ordered.setdefault(key, value)
    for key, value in foreign.items():
        ordered.setdefault(key, value)

    with _open_text(path, "w") as handle:
        _cobra_yaml.dump(ordered, handle)
