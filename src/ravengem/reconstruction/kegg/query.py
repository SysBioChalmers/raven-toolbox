"""De-novo KEGG draft from a proteome FASTA via HMM search (step 3b.5).

Ports the FASTA/HMM branch of RAVEN ``getKEGGModelForOrganism``: search a query
proteome against the KO profile-HMM library (3b.3), assign genes to KOs using the
score cut-off and the two score-ratio filters, then build the draft model with the
shared assembler. For organisms not in KEGG.

Improvement over RAVEN: one ``hmmscan`` against the single ``hmmpress``-ed library
(K7) replaces RAVEN's per-KO ``hmmsearch`` loop. Phylogenetic-distance subsampling
is **not** used — our prebuilt prok90/euk90 libraries already fix the sequence set,
so picking the right domain library (not per-organism distance weighting) is the
relevant choice.

The scoring/assignment logic (:func:`assign_kos`, :func:`parse_hmmscan_tblout`) is
pure and unit-tested; running the search needs HMMER (``hmmscan``).
"""
from __future__ import annotations

import math
import subprocess
import tempfile
from pathlib import Path

import cobra
import pandas as pd

from ravengem.binaries import resolve_binary
from ravengem.io.yaml import read_yaml_model
from ravengem.reconstruction.kegg.assemble import assemble_model_from_ko_genes
from ravengem.reconstruction.kegg.parse import read_kegg_table

_NOTE = "Included by get_kegg_model_from_sequences (using HMMs)"
_MIN_EVALUE = 1e-250  # floor for a reported E-value of 0, to keep logs finite


def run_hmmscan(
    fasta: str | Path,
    library: str | Path,
    *,
    threads: int = 1,
    hmmscan: str | Path | None = None,
) -> str:
    """Run ``hmmscan`` of ``fasta`` against the pressed ``library``; return tblout text."""
    exe = resolve_binary("hmmscan", binary=hmmscan)
    with tempfile.TemporaryDirectory() as tmp:
        tbl = Path(tmp) / "hits.tbl"
        cmd = [exe, "--cpu", str(threads), "--tblout", str(tbl), str(library), str(fasta)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"hmmscan failed:\n{(proc.stderr or '').strip()}")
        return tbl.read_text()


def parse_hmmscan_tblout(text: str) -> pd.DataFrame:
    """Parse ``hmmscan --tblout`` text into a ``[ko, gene, evalue]`` table.

    In ``hmmscan`` the HMM database is the *target*, so column 1 (target name) is
    the KO, column 3 (query name) is the proteome gene, and column 5 is the
    full-sequence E-value.
    """
    rows = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 5:
            continue
        rows.append((fields[0], fields[2], float(fields[4])))
    return pd.DataFrame(rows, columns=["ko", "gene", "evalue"])


def assign_kos(
    hits: pd.DataFrame,
    *,
    cutoff: float = 1e-30,
    min_score_ratio_ko: float = 0.3,
    min_score_ratio_g: float = 0.9,
) -> dict[str, list[str]]:
    """Assign genes to KOs from HMM hits, applying the cut-off and ratio filters.

    Ports RAVEN's three steps on the KO×gene E-value matrix:

    1. keep hits with ``evalue <= cutoff``;
    2. **min_score_ratio_ko** — within a KO, drop genes whose
       ``log(evalue)/log(best_evalue_in_KO) < min_score_ratio_ko`` (prune weak
       members of a KO);
    3. **min_score_ratio_g** — within a gene, drop KOs whose
       ``log(evalue)/log(best_evalue_for_gene) < min_score_ratio_g`` (stop a gene
       that clearly belongs to one KO leaking into weaker ones).

    Smaller E-value = better; since all kept values are ``< 1`` their logs are
    negative, so the best (smallest) hit gives ratio 1 and weaker hits give a
    smaller positive ratio.

    Default calibration (see IMPROVEMENTS K15). Cross-validated against the true
    KEGG gene→KO annotation of four organisms spanning the prok/euk libraries and
    the well-/lesser-studied axis (*S. cerevisiae*, *Cyanidioschyzon merolae*,
    *E. coli*, *Mycoplasma genitalium*): real annotations score
    overwhelmingly (median E ≈ 1e-100…1e-155) while spurious hits pile up at
    ≈1e-8, so the two are separated by ~20 orders of magnitude. RAVEN's
    ``1e-50`` sits inside the *true* tail and silently drops real but divergent
    hits — costing 16% gene→KO recall on the divergent minimal genome
    (*M. genitalium*) for no noise-rejection benefit (noise is far weaker). The
    default is therefore loosened to **1e-30** (recovers that tail; still ~22
    orders above the noise floor), with the precision work moved to
    **min_score_ratio_g = 0.9** — the *effective* precision lever (it resolves
    multi-KO genes). ``min_score_ratio_ko`` proved empirically inert across all
    four organisms (identical output at 0.0/0.3/0.5) and is kept only for RAVEN
    parity.
    """
    # The ratio filters compare log(evalue)/log(best_evalue); when best == 1.0
    # the denominator is 0 → ZeroDivisionError. The default cutoff (1e-30) keeps
    # us safely away, but a caller-passed cutoff ≥ 1 is ambiguous and would
    # crash later. Reject it up front with a clear message.
    if cutoff >= 1:
        raise ValueError(
            f"cutoff must be < 1 (smaller E-value = better hit); got {cutoff!r}."
        )

    # Best (smallest) E-value per (ko, gene), filtered at the cut-off.
    mat: dict[str, dict[str, float]] = {}
    for ko, gene, evalue in zip(hits["ko"], hits["gene"], hits["evalue"], strict=True):
        if evalue > cutoff:
            continue
        e = evalue if evalue > 0 else _MIN_EVALUE
        per_ko = mat.setdefault(ko, {})
        if gene not in per_ko or e < per_ko[gene]:
            per_ko[gene] = e

    # Step 2: prune weak genes within each KO.
    for ko, genes in mat.items():
        log_best = math.log(min(genes.values()))
        mat[ko] = {
            g: e for g, e in genes.items() if math.log(e) / log_best >= min_score_ratio_ko
        }

    # Step 3: prune weak KOs within each gene (over the survivors of step 2).
    gene_kos: dict[str, dict[str, float]] = {}
    for ko, genes in mat.items():
        for g, e in genes.items():
            gene_kos.setdefault(g, {})[ko] = e
    dropped: set[tuple[str, str]] = set()
    for g, kos in gene_kos.items():
        log_best = math.log(min(kos.values()))
        for ko, e in kos.items():
            if math.log(e) / log_best < min_score_ratio_g:
                dropped.add((ko, g))

    result: dict[str, list[str]] = {}
    for ko, genes in mat.items():
        kept = sorted(g for g in genes if (ko, g) not in dropped)
        if kept:
            result[ko] = kept
    return result


def get_kegg_model_from_sequences(
    fasta: str | Path,
    reference_model: cobra.Model,
    ko_reaction: pd.DataFrame,
    library: str | Path,
    *,
    rxn_flags: pd.DataFrame | None = None,
    model_id: str | None = None,
    cutoff: float = 1e-30,
    min_score_ratio_ko: float = 0.3,
    min_score_ratio_g: float = 0.9,
    keep_spontaneous: bool = True,
    keep_undefined_stoich: bool = True,
    keep_incomplete: bool = True,
    keep_general: bool = False,
    threads: int = 1,
    hmmscan: str | Path | None = None,
) -> cobra.Model:
    """Reconstruct a draft model for a proteome by HMM-searching the KO library.

    Searches ``fasta`` against the pressed ``library`` (3b.3), assigns KOs
    (:func:`assign_kos`), and assembles the model against ``reference_model`` /
    ``ko_reaction``. Genes are the query proteome's identifiers.
    """
    hits = parse_hmmscan_tblout(run_hmmscan(fasta, library, threads=threads, hmmscan=hmmscan))
    ko_to_genes = assign_kos(
        hits,
        cutoff=cutoff,
        min_score_ratio_ko=min_score_ratio_ko,
        min_score_ratio_g=min_score_ratio_g,
    )
    model, _ = assemble_model_from_ko_genes(
        reference_model,
        ko_reaction,
        ko_to_genes,
        rxn_flags=rxn_flags,
        keep_spontaneous=keep_spontaneous,
        keep_undefined_stoich=keep_undefined_stoich,
        keep_incomplete=keep_incomplete,
        keep_general=keep_general,
        model_id=model_id,
        note=_NOTE,
    )
    return model


def get_kegg_model_from_sequences_with_artefacts(
    fasta: str | Path,
    artefact_dir: str | Path | None = None,
    library: str | Path | None = None,
    *,
    domain: str = "prokaryotes",
    version: str | None = None,
    **kwargs,
) -> cobra.Model:
    """Load reference model + tables from ``artefact_dir`` and run the HMM query.

    If ``artefact_dir`` / ``library`` are ``None`` they are fetched/cached via
    :func:`ravengem.data.ensure_kegg_data` / :func:`ravengem.data.ensure_kegg_hmm_library`
    (``domain`` selects the prok/euk library; ``version`` the release).
    """
    if artefact_dir is None or library is None:
        from ravengem.data import ensure_kegg_data, ensure_kegg_hmm_library

        if artefact_dir is None:
            artefact_dir = ensure_kegg_data(version=version)
        if library is None:
            library = ensure_kegg_hmm_library(domain, version=version)
    artefact_dir = Path(artefact_dir)
    reference_model = read_yaml_model(artefact_dir / "reference_model.yml.gz")
    ko_reaction = read_kegg_table(artefact_dir / "ko_reaction.tsv.gz")
    rxn_flags = read_kegg_table(artefact_dir / "rxn_flags.tsv.gz")
    return get_kegg_model_from_sequences(
        fasta, reference_model, ko_reaction, library, rxn_flags=rxn_flags, **kwargs
    )
