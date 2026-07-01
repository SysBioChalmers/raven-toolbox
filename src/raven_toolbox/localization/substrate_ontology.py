"""ChEBI-ontology substrate matching — the substrate-specific layer of evidence-aware transport scoring.

Two curated raven-data artefacts underlie this (built by ``scripts/build_transporter_data.py``):

* ``tcdb_substrates.tsv`` — TCDB's curated ``TC-ID -> substrate ChEBI(s)`` table;
* ``chebi_relations.tsv.gz`` — the ChEBI ``is_a`` graph plus the protonation/tautomer bridges.

Together they let a model metabolite's ChEBI be matched to the *specific* substrate a transporter is
curated to carry — a graded, distance-decayed refinement of the coarse family→class match. A metabolite
that **is** the curated substrate scores 1.0; a near subtype or protonation/tautomer variant scores a
little less; a far, generic relative (everything is eventually ``is_a`` "organic anion") falls outside
the hop budget and scores 0, so the coarse class layer takes over. This lifts *recall* — it evidences
metabolites the name-keyword classifier (:func:`~.transport_evidence.default_substrate_of`) misses —
without the coarse layer's promiscuity.
"""
from __future__ import annotations

import gzip
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path

from raven_toolbox.data import ensure_data_file

__all__ = ["SubstrateOntology", "load_tc_substrates"]

# ChEBI relations traversed when rolling a metabolite up toward a transporter substrate. ``is_a`` is
# directed (specific -> general); the protonation/tautomer relations are chemical equivalences, so they
# are followed in both directions (a model anion <-> the TCDB acid/base form).
_SYMMETRIC = frozenset({"is_conjugate_base_of", "is_conjugate_acid_of", "is_tautomer_of"})
# Evidence weight by hop distance to the nearest curated substrate; beyond the last entry the relative
# is too generic to trust and the match is dropped (coarse class takes over). Tunable.
_WEIGHT_BY_HOP: tuple[float, ...] = (1.0, 0.9, 0.8, 0.65, 0.5)
_MAX_HOPS = len(_WEIGHT_BY_HOP) - 1


def load_tc_substrates(path: str | Path | None = None) -> dict[str, frozenset[str]]:
    """Load ``tcdb_substrates.tsv`` -> ``{TC-ID: frozenset(ChEBI)}`` (auto-downloaded if ``path`` is None)."""
    path = ensure_data_file("transporters", "tcdb_substrates.tsv") if path is None else Path(path)
    out: dict[str, frozenset[str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        tc, _, chebis = line.partition("\t")
        if chebis:
            out[tc] = frozenset(chebis.split(";"))
    return out


class SubstrateOntology:
    """The ChEBI graph + TCDB substrate table, with a graded metabolite→substrate match.

    Build one with :meth:`load` (downloads/caches both artefacts) and pass it to
    :func:`~.transport_evidence.evidence_aware_transport_cost`. Reachability closures are memoised, so
    scoring a whole model reuses each metabolite ChEBI's roll-up.
    """

    def __init__(self, edges: dict[str, set[str]], tc_substrates: dict[str, frozenset[str]],
                 alt: dict[str, str] | None = None):
        self._edges = edges
        self._tc = tc_substrates
        self._alt = alt or {}  # secondary/deprecated ChEBI id -> connected primary id
        self._reach_cache: dict[str, dict[str, int]] = {}

    @classmethod
    def load(
        cls,
        *,
        relations_path: str | Path | None = None,
        substrates_path: str | Path | None = None,
    ) -> SubstrateOntology:
        relations_path = (relations_path if relations_path is not None
                          else ensure_data_file("transporters", "chebi_relations.tsv.gz"))
        edges: dict[str, set[str]] = defaultdict(set)
        alt: dict[str, str] = {}
        with gzip.open(relations_path, "rt", encoding="utf-8") as fh:
            for line in fh:
                child, rel, parent = line.rstrip("\n").split("\t")
                if rel == "alt_id":
                    alt[child] = parent  # secondary -> primary; applied before any graph walk
                    continue
                edges[child].add(parent)
                if rel in _SYMMETRIC:
                    edges[parent].add(child)
        return cls(dict(edges), load_tc_substrates(substrates_path), alt)

    def substrates_of(self, tc_ids: Iterable[str]) -> frozenset[str]:
        """Union of curated substrate ChEBIs over a set of (full 5-level) TCDB TC-IDs."""
        out: set[str] = set()
        for tc in tc_ids:
            out |= self._tc.get(tc, frozenset())
        return frozenset(out)

    def _reach(self, chebi: str) -> dict[str, int]:
        """``{reachable ChEBI: min hop}`` within the hop budget (self at hop 0), memoised."""
        chebi = self._alt.get(chebi, chebi)  # normalise a secondary id onto its connected primary
        cached = self._reach_cache.get(chebi)
        if cached is not None:
            return cached
        reach = {chebi: 0}
        frontier = {chebi}
        for hop in range(1, _MAX_HOPS + 1):
            nxt = set().union(*(self._edges.get(x, ()) for x in frontier)) - reach.keys()
            if not nxt:
                break
            for node in nxt:
                reach[node] = hop
            frontier = nxt
        self._reach_cache[chebi] = reach
        return reach

    def match(self, metabolite_chebis: Iterable[str], substrate_chebis: Iterable[str]) -> float:
        """Graded [0, 1] match: the strongest (nearest) roll-up from any metabolite ChEBI to any of the
        transporter's curated substrate ChEBIs; 0 when none is within the hop budget."""
        subs = frozenset(self._alt.get(s, s) for s in substrate_chebis)  # normalise secondary ids
        if not subs:
            return 0.0
        best = 0.0
        for cm in metabolite_chebis:
            reach = self._reach(cm)
            hits = [reach[s] for s in subs if s in reach]
            if hits:
                weight = _WEIGHT_BY_HOP[min(hits)]
                if weight > best:
                    best = weight
                    if best >= 1.0:
                        return 1.0
        return best
