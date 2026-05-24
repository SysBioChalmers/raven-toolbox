"""Build per-KO HMM libraries from KEGG sequences (step 3b.3, maintainer-side).

Ports RAVEN ``constructMultiFasta`` plus the clustering/alignment/training stages
of ``getKEGGModelForOrganism``. Run once per KEGG release to produce the
``prok90`` / ``euk90`` HMM libraries that the de-novo query path (3b.5) searches.

Per KO, within one domain (prokaryote / eukaryote):

1. **Multi-FASTA** — gather the member genes' sequences from ``genes.pep``
   (:func:`build_ko_fastas`).
2. **CD-HIT** — dereplicate near-identical sequences (default 90 % identity).
3. **MAFFT** — multiple-sequence alignment (``--auto --anysymbol``).
4. **hmmbuild** — train the profile HMM.

Finally the per-KO HMMs are concatenated and ``hmmpress``-ed into a single
searchable library (an improvement over RAVEN's per-KO ``hmmsearch``: one
``hmmscan`` against a pressed database replaces thousands of invocations).

The pure parts (FASTA indexing/grouping, command construction, CD-HIT ``-n``
choice) are unit-tested; running the binaries needs HMMER/MAFFT/CD-HIT, located
via :func:`ravengem.binaries.resolve_binary`.
"""
from __future__ import annotations

import functools
import logging
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import pandas as pd

from ravengem.binaries import resolve_binary
from ravengem.reconstruction.kegg.taxonomy import organisms_in_domain

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# Step 1 — per-KO multi-FASTA (constructMultiFasta)
# --------------------------------------------------------------------------- #
def _full_id(organism: str, gene: str) -> str:
    """The genes.pep header key for a gene, i.e. ``organism:gene``."""
    return f"{organism}:{gene}"


def _index_fasta(path: str | Path, wanted: set[str]) -> dict[str, tuple[int, int]]:
    """Map each wanted record id to its ``(start, end)`` byte span in ``path``.

    The record id is the first whitespace-delimited token of the ``>`` header.
    One streaming pass; only wanted ids are kept (memory stays small).
    """
    index: dict[str, tuple[int, int]] = {}
    cur_id: str | None = None
    cur_start = 0
    pos = 0
    with open(path, "rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if cur_id is not None and cur_id in wanted:
                    index[cur_id] = (cur_start, pos)
                cur_id = line[1:].split(None, 1)[0].decode()
                cur_start = pos
            pos += len(line)
    if cur_id is not None and cur_id in wanted:
        index[cur_id] = (cur_start, pos)
    return index


def build_ko_fastas(
    organism_gene_ko: pd.DataFrame,
    genes_pep: str | Path,
    out_dir: str | Path,
    *,
    organisms: set[str] | None = None,
) -> dict[str, Path]:
    """Write one ``<KO>.fa`` per KO with its member genes' sequences.

    Port of RAVEN ``constructMultiFasta``, but with a stdlib offset index instead
    of the Java-hashtable byte scan. ``organisms`` restricts to a domain's
    organism codes (for the prok/euk split). Empty KOs are skipped (no file).
    Returns ``{ko: path}`` for the files written.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = organism_gene_ko
    if organisms is not None:
        rows = rows[rows["organism"].isin(organisms)]

    ko_to_ids: dict[str, list[str]] = {}
    wanted: set[str] = set()
    for organism, gene, ko in zip(rows["organism"], rows["gene"], rows["ko"], strict=True):
        fid = _full_id(organism, gene)
        ko_to_ids.setdefault(ko, []).append(fid)
        wanted.add(fid)

    index = _index_fasta(genes_pep, wanted)

    written: dict[str, Path] = {}
    with open(genes_pep, "rb") as src:
        for ko, ids in ko_to_ids.items():
            present = sorted({i for i in ids if i in index})
            if not present:
                continue
            path = out_dir / f"{ko}.fa"
            with open(path, "wb") as out:
                for fid in present:
                    start, end = index[fid]
                    src.seek(start)
                    out.write(src.read(end - start))
            written[ko] = path
    return written


# --------------------------------------------------------------------------- #
# Steps 2-4 — cluster, align, train (one KO)
# --------------------------------------------------------------------------- #
def _cdhit_word_size(seq_identity: float) -> str:
    """CD-HIT ``-n`` word size for a given identity threshold (per CD-HIT guide)."""
    if not 0.4 < seq_identity <= 1.0:
        raise ValueError("seq_identity must be in (0.4, 1.0] (or -1 to skip CD-HIT).")
    if seq_identity > 0.7:
        return "5"
    if seq_identity > 0.6:
        return "4"
    if seq_identity > 0.5:
        return "3"
    return "2"


def _count_sequences(fasta: Path) -> int:
    with open(fasta, "rb") as fh:
        return sum(1 for line in fh if line.startswith(b">"))


def _fasta_stats(fasta: Path) -> tuple[int, int]:
    """Return ``(sequence_count, total_residues)`` in one pass."""
    n = residues = 0
    with open(fasta, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                n += 1
            else:
                residues += len(line.strip())
    return n, residues


def _subsample_fasta(src: Path, dst: Path, n: int) -> int:
    """Write ``n`` evenly-spaced records of ``src`` to ``dst`` (caps huge KOs).

    Even spacing (rather than the first ``n``) spreads the kept sequences across
    the file, which is ordered by organism — a phylogeny-agnostic stand-in for
    RAVEN's ``nSequences`` most-related subsampling. Returns the count written.
    """
    records: list[list[str]] = []
    with open(src) as fh:
        for line in fh:
            if line.startswith(">"):
                records.append([])
            if records:
                records[-1].append(line)
    if len(records) <= n:
        shutil.copyfile(src, dst)
        return len(records)
    step = len(records) / n
    with open(dst, "w") as out:
        for i in range(n):
            out.write("".join(records[int(i * step)]))
    return n


def _cdhit_cmd(cdhit: str, inp: Path, out: Path, seq_identity: float, threads: int) -> list[str]:
    return [
        cdhit, "-i", str(inp), "-o", str(out),
        "-c", str(seq_identity), "-n", _cdhit_word_size(seq_identity),
        "-M", "2000", "-T", str(threads),
    ]


# MAFFT uses fast progressive FFT-NS-2 until an alignment's total *residues*
# (count x length — what actually drives memory) would push it past the RAM
# available to MAFFT, then switches to memory-light PartTree (which keeps all
# sequences; only the guide tree is approximated). The residue budget is derived
# from system memory (:func:`_auto_residue_budget`), calibrated on real KEGG 118:
# FFT-NS-2 used ~2 GB of MAFFT RSS per 1 M residues (K00901: 2.55 M residues
# needed ~5 GB; the same alignment was 0.69 GB under PartTree).
_FFTNS2_BYTES_PER_RESIDUE = 2000
# RAM not available to MAFFT (OS + WSL2 + Python); on WSL2 the reported VM total
# overcounts usable RAM by ~1.7-2.5 GB, so reserve the high end.
_MAFFT_MEMORY_OVERHEAD = int(2.5 * 1024**3)
_MEMORY_SAFETY = 0.7  # leave headroom; never budget MAFFT to the brink
_DEFAULT_RESIDUE_BUDGET = 1_000_000  # fallback when total memory can't be detected
_LOW_MEMORY_BYTES = 16 * 1024**3  # below this, warn that the budget is conservative


def _total_memory_bytes() -> int | None:
    try:
        return os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    except (AttributeError, ValueError, OSError):
        return None


@functools.lru_cache(maxsize=1)
def _auto_residue_budget() -> int:
    """Residue budget for the FFT-NS-2 -> PartTree cutover, from available RAM.

    Computed and logged once. Warns when memory is limited: the budget is then
    conservative, so more KOs use the approximate PartTree alignment.
    """
    total = _total_memory_bytes()
    if total is None:
        logger.warning(
            "Could not detect system memory; using default MAFFT residue budget %d. "
            "Pass parttree_residues to set it explicitly.", _DEFAULT_RESIDUE_BUDGET,
        )
        return _DEFAULT_RESIDUE_BUDGET
    mafft_mem = max(total - _MAFFT_MEMORY_OVERHEAD, 512 * 1024**2)
    budget = int(_MEMORY_SAFETY * mafft_mem / _FFTNS2_BYTES_PER_RESIDUE)
    logger.info(
        "MAFFT residue budget %d auto-set from %.1f GB total RAM (~%.1f GB for MAFFT)",
        budget, total / 1024**3, mafft_mem / 1024**3,
    )
    if total < _LOW_MEMORY_BYTES:
        logger.warning(
            "Limited memory (%.1f GB total): MAFFT residue budget set conservatively to "
            "%d, so more KOs will use the approximate PartTree alignment. On a machine "
            "with more RAM, raise parttree_residues (or add memory) to keep more KOs on "
            "the more accurate FFT-NS-2.", total / 1024**3, budget,
        )
    return budget


def _mafft_cmd(
    mafft: str, inp: Path, threads: int, *, fast: bool = True, parttree: bool = False
) -> list[str]:
    """Build the MAFFT command.

    ``fast`` selects FFT-NS-2 (``--retree 2 --maxiterate 0``) — fast progressive
    alignment, the right trade-off for building profile HMMs — instead of
    ``--auto`` (which picks slow iterative refinement on medium/large inputs).
    ``parttree`` adds MAFFT's PartTree approximation for very large inputs.
    """
    cmd = [mafft]
    if parttree:
        cmd += ["--retree", "2", "--parttree"]
    elif fast:
        cmd += ["--retree", "2", "--maxiterate", "0"]
    else:
        cmd += ["--auto"]
    cmd += ["--anysymbol", "--thread", str(threads), str(inp)]
    return cmd


def _hmmbuild_cmd(
    hmmbuild: str, out_hmm: Path, aligned: Path, threads: int, name: str | None = None
) -> list[str]:
    cmd = [hmmbuild, "--cpu", str(threads)]
    if name:  # name the profile after its KO so hmmscan targets are KO ids
        cmd += ["-n", name]
    cmd += [str(out_hmm), str(aligned)]
    return cmd


def _run(cmd: list[str], *, stdout_path: Path | None = None) -> str:
    """Run a command; optionally redirect stdout to a file. Raises on failure."""
    if stdout_path is not None:
        with open(stdout_path, "w") as out:
            proc = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
        stderr = proc.stderr or ""
    else:
        proc = subprocess.run(cmd, capture_output=True, text=True)
        stderr = proc.stderr or ""
    if proc.returncode != 0:
        raise RuntimeError(f"{Path(cmd[0]).name} failed:\n{stderr.strip()}")
    return stderr


def _staged_run(
    cmd: list[str], *, label: str, stage: str, verbose: bool, stdout_path: Path | None = None
) -> str:
    """Run a stage's command, logging which tool runs for which KO (when verbose).

    At INFO: ``[KO] stage: tool …`` before, and the elapsed time after. The tool's
    own stderr (MAFFT/CD-HIT/hmmbuild progress) is logged at DEBUG.
    """
    if verbose:
        logger.info("[%s] %s: running %s", label, stage, Path(cmd[0]).name)
    start = time.perf_counter()
    stderr = _run(cmd, stdout_path=stdout_path)
    if verbose:
        logger.info("[%s] %s: done in %.1fs", label, stage, time.perf_counter() - start)
        if stderr.strip():
            logger.debug("[%s] %s output:\n%s", label, stage, stderr.strip())
    return stderr


def build_ko_hmm(
    ko_fasta: str | Path,
    out_hmm: str | Path,
    *,
    seq_identity: float = 0.9,
    max_sequences: int | None = None,
    parttree_residues: int | None = None,
    threads: int = 1,
    fast: bool = True,
    verbose: bool = False,
    cdhit: str | Path | None = None,
    mafft: str | Path | None = None,
    hmmbuild: str | Path | None = None,
) -> Path:
    """Cluster, align and train a profile HMM for one KO's multi-FASTA.

    Single-sequence KOs skip CD-HIT/MAFFT (a lone sequence is its own alignment).
    ``seq_identity=-1`` skips CD-HIT. ``max_sequences`` caps the post-CD-HIT set
    (evenly subsampled) to bound MAFFT time/memory on over-represented KOs.
    ``fast`` uses MAFFT FFT-NS-2 (fast progressive) rather than ``--auto``'s slow
    iterative refinement. MAFFT switches to memory-light PartTree once the
    alignment exceeds ``parttree_residues`` total residues; when that is ``None``
    the budget is derived from available RAM (:func:`_auto_residue_budget`).
    ``verbose`` logs (via the ``logging`` module, INFO/DEBUG) which tool is running
    for this KO, sequence counts at each stage, timings, and the tools' own
    output. Returns ``out_hmm``.
    """
    ko_fasta = Path(ko_fasta)
    out_hmm = Path(out_hmm)
    label = out_hmm.stem
    out_hmm.parent.mkdir(parents=True, exist_ok=True)
    n = _count_sequences(ko_fasta)
    if n == 0:
        raise ValueError(f"{ko_fasta} contains no sequences.")
    if verbose:
        logger.info("[%s] start: %d sequences", label, n)

    hmmbuild = resolve_binary("hmmbuild", binary=hmmbuild)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        if n == 1:
            if verbose:
                logger.info("[%s] single sequence: skipping CD-HIT/MAFFT", label)
            aligned = ko_fasta  # trivially aligned
        else:
            clustered = ko_fasta
            if seq_identity != -1:
                clustered = tmp / "clustered.fa"
                _staged_run(
                    _cdhit_cmd(
                        resolve_binary("cd-hit", binary=cdhit), ko_fasta, clustered,
                        seq_identity, threads,
                    ),
                    label=label, stage=f"CD-HIT ({seq_identity})", verbose=verbose,
                )
            n_clustered = _count_sequences(clustered)
            if verbose and seq_identity != -1:
                logger.info("[%s] CD-HIT: %d -> %d sequences", label, n, n_clustered)
            if max_sequences and n_clustered > max_sequences:
                capped = tmp / "capped.fa"
                _subsample_fasta(clustered, capped, max_sequences)
                if verbose:
                    logger.info("[%s] capped %d -> %d sequences", label, n_clustered, max_sequences)
                clustered = capped
            aligned = tmp / "aligned.fa"
            n_clustered, residues = _fasta_stats(clustered)
            if n_clustered == 1:
                if verbose:
                    logger.info("[%s] one sequence after CD-HIT: skipping MAFFT", label)
                shutil.copyfile(clustered, aligned)  # MAFFT can't align a single seq
            else:
                # Memory tracks residues (count x length): use PartTree once the
                # alignment would exceed the residue budget (auto from RAM if unset).
                budget = _auto_residue_budget() if parttree_residues is None else parttree_residues
                parttree = residues > budget
                _staged_run(
                    _mafft_cmd(
                        resolve_binary("mafft", binary=mafft), clustered, threads,
                        fast=fast, parttree=parttree,
                    ),
                    label=label,
                    stage=f"MAFFT {'PartTree' if parttree else 'FFT-NS-2' if fast else 'auto'} "
                    f"({n_clustered} seqs, {residues} residues)",
                    verbose=verbose,
                    stdout_path=aligned,
                )
        _staged_run(
            _hmmbuild_cmd(hmmbuild, out_hmm, aligned, threads, name=label),
            label=label, stage="hmmbuild", verbose=verbose,
        )
    if verbose:
        logger.info("[%s] complete -> %s", label, out_hmm)
    return out_hmm


# --------------------------------------------------------------------------- #
# Orchestration — a full domain library
# --------------------------------------------------------------------------- #
def build_hmm_library(
    organism_gene_ko: pd.DataFrame,
    genes_pep: str | Path,
    taxonomy: str | Path,
    out_dir: str | Path,
    *,
    domain: str,
    seq_identity: float = 0.9,
    max_sequences: int | None = None,
    parttree_residues: int | None = None,
    threads: int = 1,
    fast: bool = True,
    verbose: bool = False,
    press: bool = True,
    cdhit: str | Path | None = None,
    mafft: str | Path | None = None,
    hmmbuild: str | Path | None = None,
    hmmpress: str | Path | None = None,
) -> dict[str, Path | list[Path]]:
    """Build a domain (``"prokaryotes"``/``"eukaryotes"``) HMM library.

    Restricts genes to the domain's organisms (from ``taxonomy``), builds a
    multi-FASTA and a profile HMM per KO under ``out_dir``, and (if ``press``)
    concatenates them into ``out_dir/library.hmm`` and ``hmmpress``-es it for fast
    ``hmmscan`` querying. Returns ``{"hmms": [...], "library": path | None}``.

    Heavy and binary-dependent — intended for the maintainer, run once per KEGG
    release. Skips KOs that already have an ``.hmm`` (resumable).
    """
    out_dir = Path(out_dir)
    fasta_dir = out_dir / "fasta"
    hmm_dir = out_dir / "hmms"
    hmm_dir.mkdir(parents=True, exist_ok=True)

    organisms = organisms_in_domain(taxonomy, domain)
    if not organisms:
        raise ValueError(f"No organisms found for domain {domain!r} in {taxonomy}.")

    ko_fastas = build_ko_fastas(organism_gene_ko, genes_pep, fasta_dir, organisms=organisms)

    hmms: list[Path] = []
    for ko, fasta in ko_fastas.items():
        out_hmm = hmm_dir / f"{ko}.hmm"
        if not out_hmm.exists():
            build_ko_hmm(
                fasta, out_hmm, seq_identity=seq_identity, max_sequences=max_sequences,
                parttree_residues=parttree_residues, threads=threads, fast=fast,
                verbose=verbose, cdhit=cdhit, mafft=mafft, hmmbuild=hmmbuild,
            )
        hmms.append(out_hmm)

    library: Path | None = None
    if press and hmms:
        library = out_dir / "library.hmm"
        with open(library, "wb") as out:
            for hmm in sorted(hmms):
                with open(hmm, "rb") as fh:
                    shutil.copyfileobj(fh, out)
        _run([resolve_binary("hmmpress", binary=hmmpress), "-f", str(library)])

    return {"hmms": hmms, "library": library}
