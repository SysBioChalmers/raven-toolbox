"""Build a draft model for a KEGG species from the reference artefacts (step 3b.4).

Ports the **organism-ID** path of RAVEN ``getKEGGModelForOrganism`` (the branch
taken when no FASTA file is given). For an organism already annotated in KEGG it
needs no homology search: take the organism's gene↔KO assignments, map KO→reaction
against the gene-free reference model, OR-join the organism's genes into each
reaction's GPR, and keep the reactions that end up with genes (plus spontaneous
reactions, optionally). The HMM/FASTA path is step 3b.5.

Consumes the 3b.2 artefacts: the gene-free reference ``cobra.Model`` plus the
``ko_reaction``, ``organism_gene_ko`` and ``rxn_flags`` tables. The KO→reaction
mapping is taken from the ``ko_reaction`` table (a lossless published artefact)
rather than from the reference model's annotations, so it does not depend on KEGG
annotations surviving an SBML round-trip.
"""
from __future__ import annotations

from pathlib import Path

import cobra
import pandas as pd

from ravengem.reconstruction.kegg.parse import read_kegg_table

_NOTE = "Included by get_kegg_model_for_organism (no HMMs)"
_DOMAINS = {"eukaryotes", "prokaryotes"}


def _flag_set(rxn_flags: pd.DataFrame | None, column: str) -> set[str]:
    """Reaction ids whose ``column`` flag is truthy (handles bool or TSV strings)."""
    if rxn_flags is None or column not in rxn_flags:
        return set()
    mask = rxn_flags[column].map(lambda v: str(v).strip().lower() in ("true", "1"))
    return set(rxn_flags.loc[mask, "reaction"])


def get_kegg_model_for_organism(
    organism_id: str,
    reference_model: cobra.Model,
    ko_reaction: pd.DataFrame,
    organism_gene_ko: pd.DataFrame,
    *,
    rxn_flags: pd.DataFrame | None = None,
    keep_spontaneous: bool = True,
    keep_undefined_stoich: bool = True,
    keep_incomplete: bool = True,
    keep_general: bool = False,
) -> cobra.Model:
    """Reconstruct a draft model for a KEGG species from its KO annotations.

    Parameters
    ----------
    organism_id
        Three/four-letter KEGG organism code (e.g. ``"eco"``). Matched
        case-insensitively against the ``organism`` column.
    reference_model
        The gene-free KEGG reference model (from :func:`build_reference_model`).
    ko_reaction, organism_gene_ko, rxn_flags
        The relational tables from :func:`build_kegg_tables` (or read back with
        :func:`read_kegg_table`).
    keep_spontaneous, keep_undefined_stoich, keep_incomplete, keep_general
        Quality filters (RAVEN's ``keep*``). A reaction flagged in ``rxn_flags``
        is dropped unless its keep flag is set; this takes precedence over having
        genes. Spontaneous reactions are additionally kept *without* genes when
        ``keep_spontaneous`` is true.

    Returns
    -------
    cobra.Model
        A copy of the reference restricted to the organism's reactions, with GPRs
        built and ``kegg.genes`` annotations on the genes.
    """
    org = organism_id.lower()
    if org in _DOMAINS:
        raise NotImplementedError(
            "Domain-wide models ('eukaryotes'/'prokaryotes') need the "
            "phylogenetic-distance table from getPhylDist (step 3b.5); not yet "
            "implemented. Pass a species code such as 'eco'."
        )
    known = set(organism_gene_ko["organism"].str.lower())
    if org not in known:
        raise ValueError(
            f"Organism '{organism_id}' has no genes in organism_gene_ko. "
            f"Provide a KEGG species code present in the table."
        )

    # reaction -> set of KOs (from the lossless table, not model annotations).
    rxn_to_kos: dict[str, set[str]] = {}
    for ko, rid in zip(ko_reaction["ko"], ko_reaction["reaction"], strict=True):
        rxn_to_kos.setdefault(rid, set()).add(ko)

    # KO -> this organism's genes.
    sub = organism_gene_ko[organism_gene_ko["organism"].str.lower() == org]
    ko_to_genes: dict[str, list[str]] = {}
    for ko, gene in zip(sub["ko"], sub["gene"], strict=True):
        ko_to_genes.setdefault(ko, []).append(gene)

    spontaneous = _flag_set(rxn_flags, "spontaneous")
    drop_if = {
        "undefined_stoich": (keep_undefined_stoich, _flag_set(rxn_flags, "undefined_stoich")),
        "incomplete": (keep_incomplete, _flag_set(rxn_flags, "incomplete")),
        "general": (keep_general, _flag_set(rxn_flags, "general")),
    }

    gpr_map: dict[str, list[str]] = {}
    spontaneous_kept: set[str] = set()
    for rxn in reference_model.reactions:
        rid = rxn.id
        # Quality filters first: a filtered-out reaction is dropped even if it
        # would have genes (matches RAVEN's load-time pruning).
        if any(not keep_flag and rid in flagged for keep_flag, flagged in drop_if.values()):
            continue
        genes = sorted({g for ko in rxn_to_kos.get(rid, ()) for g in ko_to_genes.get(ko, ())})
        if genes:
            gpr_map[rid] = genes
        elif rid in spontaneous and keep_spontaneous:
            spontaneous_kept.add(rid)

    keep = set(gpr_map) | spontaneous_kept
    model = reference_model.copy()
    model.id = organism_id
    model.name = f"Generated by get_kegg_model_for_organism for {organism_id}"
    model.remove_reactions(
        [r for r in model.reactions if r.id not in keep], remove_orphans=True
    )
    for rid, genes in gpr_map.items():
        model.reactions.get_by_id(rid).gene_reaction_rule = " or ".join(genes)
    for rid in keep:
        model.reactions.get_by_id(rid).notes["note"] = _NOTE
    for gene in model.genes:
        gene.annotation["kegg.genes"] = f"{org}:{gene.id}"
    return model


def get_kegg_model_for_organism_from_artefacts(
    organism_id: str, artefact_dir: str | Path, **kwargs
) -> cobra.Model:
    """Load the published 3b.2 artefacts from ``artefact_dir`` and build the model.

    Reads ``reference_model.xml`` and the ``ko_reaction``/``organism_gene_ko``/
    ``rxn_flags`` gzipped-TSV tables, then calls :func:`get_kegg_model_for_organism`.
    """
    artefact_dir = Path(artefact_dir)
    reference_model = cobra.io.read_sbml_model(str(artefact_dir / "reference_model.xml"))
    ko_reaction = read_kegg_table(artefact_dir / "ko_reaction.tsv.gz")
    organism_gene_ko = read_kegg_table(artefact_dir / "organism_gene_ko.tsv.gz")
    rxn_flags = read_kegg_table(artefact_dir / "rxn_flags.tsv.gz")
    return get_kegg_model_for_organism(
        organism_id,
        reference_model,
        ko_reaction,
        organism_gene_ko,
        rxn_flags=rxn_flags,
        **kwargs,
    )
