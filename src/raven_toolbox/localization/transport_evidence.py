"""Evidence-aware transport scoring — turn transporter evidence into per-metabolite transport costs.

The localisation assignment (:func:`predict_localization`, :func:`assign_compartments`) charges a flat
cost per inter-compartment transport it must add. Applied blindly that is *indiscriminate*: it drops
real, functionally essential transporters as readily as spurious ones, because the cost ignores whether
a transporter actually exists (see :doc:`/studies/carvefungi_milp_benchmark`). This module makes the
cost **evidence-aware**: a transport is cheap when a transporter gene supports it (right substrate,
right membrane) and pays the full prior otherwise.

    ``transport_cost(metabolite) = base_cost * (1 - evidence(metabolite))``

:func:`evidence_aware_transport_cost` returns exactly the ``{metabolite_base: cost}`` mapping both
assignment functions already accept as their ``transport_cost`` argument — so no MILP change is needed.

Evidence is carrier-general (any transporter family/membrane) and organism-agnostic (sequence-derived;
the only per-organism input is the proteome). :func:`annotate_proteome` runs the bundled ``hmmsearch``
(Pfam transporter families) + ``diamond`` (TCDB) back-ends against a proteome FASTA;
:func:`annotate_transporters` takes a pre-computed table instead. Substrate matching is *coarse-first*
(family → coarse class); the per-substrate ChEBI layer is a later increment. See
:doc:`/reference/transport_evidence_scoring` for the full design.
"""
from __future__ import annotations

import re
import subprocess
import tempfile
from collections import defaultdict
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

import cobra
import pandas as pd

from raven_toolbox.binaries import resolve_binary
from raven_toolbox.data import ensure_data_file
from raven_toolbox.localization.transporter_tables import PFAM_TRANSPORTERS, TC_FAMILY_CLASS

__all__ = [
    "TransporterAnnotation",
    "annotate_proteome",
    "annotate_transporters",
    "evidence_aware_transport_cost",
]


@dataclass(frozen=True)
class TransporterAnnotation:
    """Per-gene transporter evidence, from any source (Pfam/hmmer, TCDB/diamond, orthology, or a
    hand-curated table). ``confidence`` is a 0..1 strength; ``substrate_classes`` are coarse shared
    classes (e.g. ``"amino_acid"``, ``"sugar"``, ``"organic_acid"``) used to match a transporter to a
    metabolite it can plausibly carry."""

    gene: str
    confidence: float = 0.0
    families: tuple[str, ...] = ()
    substrate_classes: frozenset[str] = field(default_factory=frozenset)
    mechanism: str | None = None  # "uniport" | "symport" | "antiport" | None


def _default_base(m: cobra.Metabolite) -> str:
    """Compartment-agnostic metabolite key: strip a trailing ``_<compartment>`` from the id."""
    if m.compartment and m.id.endswith(f"_{m.compartment}"):
        return m.id[: -(len(m.compartment) + 1)]
    return m.id


def annotate_transporters(
    table: pd.DataFrame,
    *,
    gene_col: str = "gene",
    confidence_col: str = "confidence",
    families_col: str | None = "families",
    substrate_col: str | None = "substrate_classes",
    mechanism_col: str | None = "mechanism",
    sep: str = ";",
) -> dict[str, TransporterAnnotation]:
    """Parse a per-gene transporter-annotation table into :class:`TransporterAnnotation` objects.

    This is the **bring-your-own** path: the table can come from any tool (eggNOG-mapper, InterProScan,
    a web service) or, in a later increment, from the bundled ``hmmsearch``/``diamond`` back-ends. List
    columns (``families``, ``substrate_classes``) may be real lists or ``sep``-joined strings. A gene
    appearing on several rows keeps its highest-confidence annotation and the union of families /
    substrate classes.
    """

    def _as_set(v) -> frozenset[str]:
        if v is None or (isinstance(v, float) and pd.isna(v)):
            return frozenset()
        if isinstance(v, str):
            return frozenset(x.strip() for x in v.split(sep) if x.strip())
        return frozenset(str(x).strip() for x in v if str(x).strip())

    out: dict[str, TransporterAnnotation] = {}
    for row in table.to_dict("records"):
        gene = str(row[gene_col])
        conf = float(row.get(confidence_col, 0.0) or 0.0)
        fams = _as_set(row.get(families_col)) if families_col else frozenset()
        subs = _as_set(row.get(substrate_col)) if substrate_col else frozenset()
        mech = row.get(mechanism_col) if mechanism_col else None
        mech = None if mech is None or (isinstance(mech, float) and pd.isna(mech)) else str(mech)
        prev = out.get(gene)
        if prev is None or conf >= prev.confidence:
            fams = tuple(sorted(fams | set(prev.families))) if prev else tuple(sorted(fams))
            subs = subs | prev.substrate_classes if prev else subs
            out[gene] = TransporterAnnotation(gene, max(conf, prev.confidence if prev else 0.0),
                                              fams, subs, mech or (prev.mechanism if prev else None))
        else:  # keep prev confidence but accumulate families/substrates
            out[gene] = TransporterAnnotation(
                gene, prev.confidence, tuple(sorted(set(prev.families) | fams)),
                prev.substrate_classes | subs, prev.mechanism)
    return out


def evidence_aware_transport_cost(
    model: cobra.Model,
    annotation: Mapping[str, TransporterAnnotation],
    gene_compartments: Mapping[str, Iterable[str]],
    *,
    substrate_of: Callable[[cobra.Metabolite], Iterable[str]] | None = None,
    base_cost: float = 0.5,
    base_metabolite: Callable[[cobra.Metabolite], str] | None = None,
) -> dict[str, float]:
    """Per-metabolite transport cost from transporter evidence, ready to pass as ``transport_cost``.

    For each (compartment-agnostic) metabolite *m*::

        evidence(m) = max over transporter genes g of g.confidence, restricted to genes that
                      (substrate) share a substrate class with m  [if ``substrate_of`` is given], and
                      (membrane)  are localised to a compartment where m occurs [if the gene has a
                                  predicted compartment].
        cost(m)     = base_cost * (1 - evidence(m))

    Parameters
    ----------
    gene_compartments:
        gene -> predicted compartment ids (the reliable DeepLoc *compartment* calls). A carrier at
        compartment *X* is taken to support transports across *X*'s boundary.
    substrate_of:
        metabolite -> its substrate class(es). Matching a transporter to a metabolite needs this; when
        omitted the match is compartment-only (coarse — every carrier at the right membrane counts).
    base_cost:
        the flat cost for an unsupported transport (recovers today's constant behaviour when no gene
        supports the metabolite).

    Returns every metabolite base -> cost (unsupported metabolites map to ``base_cost``), so the result
    is a self-contained ``transport_cost`` mapping.
    """
    base_of = base_metabolite or _default_base
    comps_of: dict[str, set[str]] = defaultdict(set)
    rep: dict[str, cobra.Metabolite] = {}
    for m in model.metabolites:
        b = base_of(m)
        if m.compartment:
            comps_of[b].add(m.compartment)
        rep.setdefault(b, m)

    carriers: list[tuple[set[str], frozenset[str], float]] = []
    for gene, ann in annotation.items():
        if ann.confidence <= 0:
            continue
        carriers.append((set(gene_compartments.get(gene, ())), ann.substrate_classes, ann.confidence))

    costs: dict[str, float] = {}
    for b, bcomps in comps_of.items():
        classes = frozenset(substrate_of(rep[b])) if substrate_of else None
        best = 0.0
        for gcomps, subs, conf in carriers:
            if classes is not None and not (classes & subs):
                continue  # substrate mismatch
            if gcomps and not (gcomps & bcomps):
                continue  # gene's membrane does not border a compartment this metabolite is in
            best = max(best, conf)
        costs[b] = base_cost * (1.0 - best)
    return costs


# --------------------------------------------------------------------------- annotation back-end
_TC_FAMILY = re.compile(r"(\d+\.[A-Z]\.\d+)\.\d+\.\d+")  # TC number -> family prefix (first 3 levels)


def _hmmsearch_families(proteome: Path, hmm_db: Path, threads: int) -> dict[str, tuple[set[str], float]]:
    """hmmsearch (transporter Pfam HMMs vs proteome) -> gene -> (Pfam accessions hit, best bitscore)."""
    with tempfile.TemporaryDirectory() as tmp:
        tbl = Path(tmp) / "hits.tbl"
        subprocess.run([resolve_binary("hmmsearch"), "--cut_ga", "--noali", "--cpu", str(threads),
                        "--tblout", str(tbl), str(hmm_db), str(proteome)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        out: dict[str, tuple[set[str], float]] = {}
        for line in tbl.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            gene, pfam_acc, score = fields[0], fields[3].split(".")[0], float(fields[5])
            fams, best = out.get(gene, (set(), 0.0))
            fams.add(pfam_acc)
            out[gene] = (fams, max(best, score))
        return out


def _diamond_tc(proteome: Path, tcdb_db: Path, threads: int) -> dict[str, tuple[set[str], float]]:
    """diamond blastp (proteome vs TCDB) -> gene -> (TC family prefixes hit, best fractional identity)."""
    with tempfile.TemporaryDirectory() as tmp:
        m8 = Path(tmp) / "hits.m8"
        subprocess.run([resolve_binary("diamond"), "blastp", "-q", str(proteome), "-d", str(tcdb_db),
                        "-o", str(m8), "--outfmt", "6", "qseqid", "sseqid", "pident", "-e", "1e-5",
                        "--max-target-seqs", "3", "-p", str(threads), "--quiet"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        out: dict[str, tuple[set[str], float]] = {}
        for line in m8.read_text().splitlines():
            if not line.strip():
                continue
            gene, sseqid, pident = line.split("\t")[:3]
            m = _TC_FAMILY.search(sseqid)
            if m is None:
                continue
            fams, best = out.get(gene, (set(), 0.0))
            fams.add(m.group(1))
            out[gene] = (fams, max(best, float(pident) / 100.0))
        return out


def annotate_proteome(
    proteome: str | Path,
    *,
    hmm_db: str | Path | None = None,
    tcdb_db: str | Path | None = None,
    threads: int = 1,
) -> dict[str, TransporterAnnotation]:
    """Annotate a proteome's transporters from sequence, via the bundled ``hmmsearch`` + ``diamond``.

    Scans the proteome FASTA (headers = gene ids, as written by :func:`prepare_deeploc_input`) against
    the transporter Pfam HMM database (``hmmsearch --cut_ga``) and the TCDB DIAMOND database
    (``diamond blastp``); both databases auto-download from raven-data if not supplied. Pfam families
    and TCDB TC-families are mapped to coarse substrate classes via
    :mod:`raven_toolbox.localization.transporter_tables`. Returns gene → :class:`TransporterAnnotation`
    (this coarse-first increment leaves ``substrate_chebi``/mechanism empty; the per-substrate ChEBI
    layer is a later increment). Genes with no transporter hit are absent from the result, so this
    doubles as a ``substrate_of``-independent membrane/family screen.
    """
    proteome = Path(proteome)
    hmm_db = Path(hmm_db) if hmm_db else ensure_data_file("transporters", "transporter_pfam.hmm")
    tcdb_db = Path(tcdb_db) if tcdb_db else ensure_data_file("transporters", "tcdb.dmnd")
    pfam = _hmmsearch_families(proteome, hmm_db, threads)
    tcdb = _diamond_tc(proteome, tcdb_db, threads)

    out: dict[str, TransporterAnnotation] = {}
    for gene in set(pfam) | set(tcdb):
        fams, pf_score = pfam.get(gene, (set(), 0.0))
        tc_fams, tc_id = tcdb.get(gene, (set(), 0.0))
        classes: set[str] = set()
        for acc in fams:
            if acc in PFAM_TRANSPORTERS:
                classes |= PFAM_TRANSPORTERS[acc][1]
        for tc in tc_fams:
            classes |= TC_FAMILY_CLASS.get(tc, frozenset())
        # Coarse confidence in [0,1]: a --cut_ga Pfam hit is curated-significant (floor 0.5, scaled by
        # bitscore, saturating ~200); a TCDB hit scales with sequence identity. Take the stronger. Tunable.
        conf = max(min(1.0, 0.5 + pf_score / 400.0) if fams else 0.0, tc_id)
        out[gene] = TransporterAnnotation(
            gene=gene, confidence=round(conf, 4), families=tuple(sorted(fams | tc_fams)),
            substrate_classes=frozenset(classes))
    return out
