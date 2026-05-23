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

import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from ravengem.binaries import resolve_binary
from ravengem.reconstruction.kegg.taxonomy import organisms_in_domain


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


def _cdhit_cmd(cdhit: str, inp: Path, out: Path, seq_identity: float, threads: int) -> list[str]:
    return [
        cdhit, "-i", str(inp), "-o", str(out),
        "-c", str(seq_identity), "-n", _cdhit_word_size(seq_identity),
        "-M", "2000", "-T", str(threads),
    ]


def _mafft_cmd(mafft: str, inp: Path, threads: int) -> list[str]:
    return [mafft, "--auto", "--anysymbol", "--thread", str(threads), str(inp)]


def _hmmbuild_cmd(hmmbuild: str, out_hmm: Path, aligned: Path, threads: int) -> list[str]:
    return [hmmbuild, "--cpu", str(threads), str(out_hmm), str(aligned)]


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


def build_ko_hmm(
    ko_fasta: str | Path,
    out_hmm: str | Path,
    *,
    seq_identity: float = 0.9,
    threads: int = 1,
    cdhit: str | Path | None = None,
    mafft: str | Path | None = None,
    hmmbuild: str | Path | None = None,
) -> Path:
    """Cluster, align and train a profile HMM for one KO's multi-FASTA.

    Single-sequence KOs skip CD-HIT/MAFFT (a lone sequence is its own alignment).
    ``seq_identity=-1`` skips CD-HIT. Returns the written ``out_hmm`` path.
    """
    ko_fasta = Path(ko_fasta)
    out_hmm = Path(out_hmm)
    out_hmm.parent.mkdir(parents=True, exist_ok=True)
    n = _count_sequences(ko_fasta)
    if n == 0:
        raise ValueError(f"{ko_fasta} contains no sequences.")

    hmmbuild = resolve_binary("hmmbuild", binary=hmmbuild)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        if n == 1:
            aligned = ko_fasta  # trivially aligned
        else:
            clustered = ko_fasta
            if seq_identity != -1:
                clustered = tmp / "clustered.fa"
                _run(_cdhit_cmd(
                    resolve_binary("cd-hit", binary=cdhit), ko_fasta, clustered,
                    seq_identity, threads,
                ))
            aligned = tmp / "aligned.fa"
            if _count_sequences(clustered) == 1:
                shutil.copyfile(clustered, aligned)  # MAFFT can't align a single seq
            else:
                _run(
                    _mafft_cmd(resolve_binary("mafft", binary=mafft), clustered, threads),
                    stdout_path=aligned,
                )
        _run(_hmmbuild_cmd(hmmbuild, out_hmm, aligned, threads))
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
    threads: int = 1,
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
                fasta, out_hmm, seq_identity=seq_identity, threads=threads,
                cdhit=cdhit, mafft=mafft, hmmbuild=hmmbuild,
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
