"""Per-reaction, multi-facet confidence — persisted in the model, ignored by plain cobra.

Attaches a small structured record to a reaction scoring how well-supported each *facet* of it is
(``localization`` in this first increment; ``equation``/``gene_association``/``reversibility`` follow the
same shape). Each facet is a :class:`ConfidenceEntry` — a continuous 0-1 ``score`` plus optional
provenance (a categorical ``level``, the ``basis`` evidence, ``method``/``source``/``note``). A reaction
carries a :class:`ReactionConfidence` (facet → entry) whose ``overall`` is the weakest facet.

**Storage.** The record lives as one JSON blob under ``reaction.notes["raven_confidence"]``. cobra
round-trips ``notes`` losslessly through YAML/JSON and, as an HTML-escaped string, through SBML; the
helpers here write ``json.dumps`` and read ``json.loads(html.unescape(...))``, so the same code works for
either format. Plain cobra ignores the key entirely — a confidence-annotated model still loads and
solves unchanged (a test invariant).

The design and roadmap (equation/gene/reversibility facets, ECO/SBO and Thiele-Palsson mapping) are in
``docs/studies/confidence_tracking.md``. Wire it in by calling :func:`score_localization_confidence` on
an :class:`~raven_toolbox.localization.AssignmentProposal`, and :func:`mark_curated` when a curator
firmly fixes a placement (e.g. after :func:`~raven_toolbox.localization.relocate_reactions`).
"""
from __future__ import annotations

import html
import json
import math
from dataclasses import dataclass, field
from typing import Any

import cobra

__all__ = [
    "ConfidenceEntry",
    "ReactionConfidence",
    "confidence_report",
    "get_confidence",
    "mark_curated",
    "read_confidence",
    "score_localization_confidence",
    "set_confidence",
]

_KEY = "raven_confidence"
_SCHEMA_VERSION = 1
# reserved top-level keys inside the stored blob that are not facets
_RESERVED = frozenset({"schema_version"})


@dataclass
class ConfidenceEntry:
    """One facet's confidence: a continuous ``score`` in [0, 1] plus optional provenance.

    ``level`` is an optional categorical band (``"curated"`` / ``"strong"`` / ``"moderate"`` / ``"weak"``
    / ``"none"``); ``basis`` names the evidence (``"deeploc"``, ``"fba-certified"``, ``"curator"``,
    ``"connectivity"``, ...); ``source`` distinguishes ``"auto"`` from ``"curator:<id>"``.
    """

    score: float
    basis: str = ""
    level: str | None = None
    method: str | None = None
    source: str | None = None
    note: str | None = None
    updated: str | None = None  # ISO date, passed in by the caller (not generated here, to stay pure)

    def __post_init__(self) -> None:
        self.score = float(max(0.0, min(1.0, self.score)))

    def to_dict(self) -> dict[str, Any]:
        d: dict[str, Any] = {"score": round(self.score, 4)}
        for k in ("level", "basis", "method", "source", "note", "updated"):
            v = getattr(self, k)
            if v:
                d[k] = v
        return d

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> ConfidenceEntry:
        return cls(score=float(d.get("score", 0.0)), basis=d.get("basis", ""), level=d.get("level"),
                   method=d.get("method"), source=d.get("source"), note=d.get("note"),
                   updated=d.get("updated"))


@dataclass
class ReactionConfidence:
    """A reaction's confidence across facets (``{facet_name: ConfidenceEntry}``)."""

    facets: dict[str, ConfidenceEntry] = field(default_factory=dict)

    @property
    def overall(self) -> float | None:
        """The weakest facet's score (the honest single number), or ``None`` when there are no facets."""
        return min((e.score for e in self.facets.values()), default=None)

    def to_dict(self) -> dict[str, Any]:
        d: dict[str, Any] = {"schema_version": _SCHEMA_VERSION}
        for name, entry in self.facets.items():
            d[name] = entry.to_dict()
        return d

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> ReactionConfidence:
        facets = {k: ConfidenceEntry.from_dict(v) for k, v in d.items()
                  if k not in _RESERVED and isinstance(v, dict)}
        return cls(facets=facets)


# --------------------------------------------------------------------------- storage (notes JSON blob)

def get_confidence(reaction: cobra.Reaction) -> ReactionConfidence:
    """Read a reaction's :class:`ReactionConfidence` from its notes (empty if none / unreadable).

    Handles both a native dict (YAML/JSON models) and a JSON string (SBML, possibly HTML-escaped).
    """
    raw = reaction.notes.get(_KEY)
    if raw is None:
        return ReactionConfidence()
    if isinstance(raw, str):
        try:
            raw = json.loads(html.unescape(raw))
        except (ValueError, TypeError):
            return ReactionConfidence()
    if not isinstance(raw, dict):
        return ReactionConfidence()
    return ReactionConfidence.from_dict(raw)


def _write(reaction: cobra.Reaction, record: ReactionConfidence) -> None:
    if record.facets:
        reaction.notes[_KEY] = json.dumps(record.to_dict())
    else:
        reaction.notes.pop(_KEY, None)


def set_confidence(reaction: cobra.Reaction, facet: str, entry: ConfidenceEntry) -> None:
    """Merge one facet's :class:`ConfidenceEntry` into the reaction's record and persist it to notes."""
    record = get_confidence(reaction)
    record.facets[facet] = entry
    _write(reaction, record)


def clear_confidence(reaction: cobra.Reaction, facet: str | None = None) -> None:
    """Drop one facet (``facet=...``) or the whole record (``facet=None``) from the reaction's notes."""
    if facet is None:
        reaction.notes.pop(_KEY, None)
        return
    record = get_confidence(reaction)
    record.facets.pop(facet, None)
    _write(reaction, record)


def read_confidence(model: cobra.Model) -> dict[str, ReactionConfidence]:
    """Every annotated reaction id -> its :class:`ReactionConfidence` (reactions with no record omitted)."""
    out: dict[str, ReactionConfidence] = {}
    for r in model.reactions:
        rc = get_confidence(r)
        if rc.facets:
            out[r.id] = rc
    return out


def confidence_report(model: cobra.Model):
    """A :class:`pandas.DataFrame` of confidence: one row per annotated reaction, one column per facet
    score, plus ``overall`` (the weakest facet). Sorted lowest-confidence first — the review queue."""
    import pandas as pd

    records = read_confidence(model)
    facet_names: list[str] = []
    for rc in records.values():
        for name in rc.facets:
            if name not in facet_names:
                facet_names.append(name)
    rows = []
    for rid, rc in records.items():
        row: dict[str, Any] = {"reaction": rid}
        for name in facet_names:
            entry = rc.facets.get(name)
            row[name] = entry.score if entry else None
        row["overall"] = rc.overall
        rows.append(row)
    cols = ["reaction", *sorted(facet_names), "overall"]
    out = pd.DataFrame(rows, columns=cols)
    if not out.empty:
        out = out.sort_values("overall", ascending=True, na_position="last", ignore_index=True)
    return out


# --------------------------------------------------------------------------- localization scorer (P1)

def _level(score: float, *, curated: bool = False) -> str:
    if curated:
        return "curated"
    if score >= 0.8:
        return "strong"
    if score >= 0.5:
        return "moderate"
    if score > 0.0:
        return "weak"
    return "none"


def mark_curated(reaction: cobra.Reaction, *, source: str = "curator", note: str | None = None,
                 updated: str | None = None) -> None:
    """Stamp a placement as **curator-verified**: localization confidence 1.0, ``level="curated"``.

    Call after a curator firmly fixes a reaction's compartment (e.g. via
    :func:`~raven_toolbox.localization.relocate_reactions`), so the decision persists in the model and
    :func:`score_localization_confidence` leaves it untouched on a later automated pass.
    """
    set_confidence(reaction, "localization",
                   ConfidenceEntry(score=1.0, level="curated", basis="curator", source=source,
                                   note=note, updated=updated))


def score_localization_confidence(model, proposal, scores, *,
                                  overwrite_curated: bool = False, updated: str | None = None) -> int:
    """Attach a ``localization`` confidence to every placement in ``proposal`` and return the count.

    The score is the DeepLoc support for the assigned compartment (the strongest scored gene's score
    there): a placement its evidence backs is high-confidence, one placed by connectivity alone
    (``no scored gene``) is low. The ``basis`` records ``deeploc`` and, when the proposal is certified,
    ``fba-certified``. Placements already marked ``curated`` are left alone unless ``overwrite_curated``.

    This is the inverse of the localization signal in
    :func:`~raven_toolbox.localization.curation_priority`: high confidence here == low review priority
    there. ``model`` is the draft the proposal was built on.
    """
    df = scores.df
    certified = bool(getattr(proposal, "certified", False))
    unplaced = set(getattr(proposal, "unplaced_reactions", ()))
    n = 0
    for rid, comps in proposal.placements.items():
        if not comps or rid not in model.reactions:
            continue
        reaction = model.reactions.get_by_id(rid)
        existing = get_confidence(reaction).facets.get("localization")
        if existing is not None and existing.level == "curated" and not overwrite_curated:
            continue
        comp = comps[0]
        genes = [g.id for g in reaction.genes if g.id in df.index]
        if genes and comp in df.columns:
            support = max((_score(df, g, comp) for g in genes), default=0.0)
            basis = "deeploc+fba-certified" if certified else "deeploc"
        else:
            support = 0.0
            basis = "connectivity"  # placed by function/connectivity, no localisation evidence
        note = "placed by connectivity; no scored gene" if rid in unplaced or not genes else None
        set_confidence(reaction, "localization",
                       ConfidenceEntry(score=support, level=_level(support), basis=basis,
                                       method="score_localization_confidence", source="auto",
                                       note=note, updated=updated))
        n += 1
    return n


def _score(df, g: str, c: str) -> float:
    if g not in df.index or c not in df.columns:
        return 0.0
    v = df.at[g, c]
    return 0.0 if v is None or (isinstance(v, float) and math.isnan(v)) else float(v)
