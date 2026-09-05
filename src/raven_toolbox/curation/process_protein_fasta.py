"""Rename protein FASTA headers from a gene-mapping table — port of RAVEN's
``processProteinFastaFile``.

Reads and writes FASTA by hand rather than via a new dependency (e.g.
Biopython), matching ``readFasta.m``/``writeFasta.m``'s own reasoning for
being hand-rolled: something this simple doesn't need one.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

__all__ = ["process_protein_fasta_file"]


def _read_fasta(path: Path) -> list[dict[str, str]]:
    """Matches readFasta.m: header without the leading '>', trimmed;
    sequence with all line breaks removed; blank lines ignored."""
    records: list[dict[str, str]] = []
    header: str | None = None
    sequence: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append({"header": header, "sequence": "".join(sequence)})
                header = line[1:].strip()
                sequence = []
            else:
                sequence.append(line)
        if header is not None:
            records.append({"header": header, "sequence": "".join(sequence)})
    return records


def _write_fasta(path: Path, records: list[dict[str, str]], line_width: int = 80) -> None:
    """Matches writeFasta.m: wraps sequences at line_width residues per line."""
    with open(path, "w", encoding="utf-8") as fh:
        for record in records:
            fh.write(f">{record['header']}\n")
            seq = record["sequence"]
            for i in range(0, len(seq), line_width):
                fh.write(seq[i:i + line_width] + "\n")


def process_protein_fasta_file(
    faa_file: str | Path,
    gene_table: pd.DataFrame | str | Path,
    header_col: str,
    *,
    output_dir: str | Path = ".",
) -> Path:
    """Rename a protein FASTA's headers using a gene mapping table.

    Replaces each sequence header with the value from ``header_col``,
    matched via the ``GenBank_protein`` accession present in the original
    FASTA header (its first whitespace-delimited token). Sequences whose
    accession is not found in ``gene_table`` keep their original header.

    Parameters
    ----------
    faa_file:
        Path to the protein FASTA file (``.faa``) to process.
    gene_table:
        Either a DataFrame, or a path to a tab-delimited file to load one
        from (e.g. from :func:`raven_toolbox.curation.get_gene_data`). Must
        contain at least the columns ``"GenBank_protein"`` and ``header_col``.
    header_col:
        Column whose values replace each FASTA header (e.g. ``"locus_tag"``).
    output_dir:
        Directory the processed FASTA is saved in. The output file name is
        the original base name with ``"_processed"`` appended before the
        extension.

    Returns
    -------
    Path
        Path to the written processed FASTA file.

    Raises
    ------
    FileNotFoundError
        If ``faa_file`` doesn't exist.
    ValueError
        If ``gene_table`` is missing ``"GenBank_protein"`` or ``header_col``.

    Notes
    -----
    A row is skipped only if its ``GenBank_protein`` is empty — not if
    ``header_col`` is, matching ``processProteinFastaFile.m`` exactly (an
    accession mapped to an empty ``header_col`` value still renames that
    sequence's header to an empty string).
    """
    faa_file = Path(faa_file)
    if not faa_file.is_file():
        raise FileNotFoundError(f"FASTA file not found: {str(faa_file)!r}.")

    if isinstance(gene_table, (str, Path)):
        gene_table = pd.read_csv(gene_table, sep="\t")

    for col in ("GenBank_protein", header_col):
        if col not in gene_table.columns:
            raise ValueError(f"gene_table is missing required column: {col!r}.")

    protein_map: dict[str, str] = {}
    for protein_id, header_value in zip(
        gene_table["GenBank_protein"], gene_table[header_col], strict=True
    ):
        if pd.isna(protein_id) or str(protein_id).strip() == "":
            continue
        protein_map[str(protein_id)] = "" if pd.isna(header_value) else str(header_value)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    processed_faa_file = output_dir / f"{faa_file.stem}_processed{faa_file.suffix}"

    records = _read_fasta(faa_file)
    matched = 0
    for record in records:
        tokens = record["header"].split()
        protein_acc = tokens[0] if tokens else ""
        if protein_acc in protein_map:
            record["header"] = protein_map[protein_acc]
            matched += 1

    _write_fasta(processed_faa_file, records)

    print(f"Processed FASTA written to: {processed_faa_file}")
    print(f"  Total sequences  : {len(records)}")
    print(f"  Renamed          : {matched}")
    print(f"  Kept (no match)  : {len(records) - matched}")

    return processed_faa_file
