#!/usr/bin/env python
"""Driver for the KEGG arm of the homology cut-off study.

Reproduces docs/studies/homology_cutoff_calibration.md: fetch proteomes keyed by
KEGG gene id, align the template against each target once, then score every
threshold combination against KEGG's own orthology (do the two genes share a KO).

    python scripts/homology_cutoff_kegg.py fetch    --out work/
    python scripts/homology_cutoff_kegg.py align    --out work/
    python scripts/homology_cutoff_kegg.py score    --out work/ \\
        --gene-ko kegg118_organism_gene_ko.tsv.gz

The alignment is the only slow step (~6-13 min per pair) and is cached, so
re-scoring with a different loss function costs nothing.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import json
import pathlib
import sys
import time
import urllib.request

TEMPLATE = "sce"

# KEGG organism code -> UniProt reference proteome. The code names a *strain*, so
# these pairings are checked by the KEGG cross-reference prefix, not assumed.
PROTEOMES = {
    "sce": "UP000002311",   # S. cerevisiae S288C -- template
    "kla": "UP000000598",   # close
    "yli": "UP000001300",   # medium
    "ani": "UP000000560",   # distant
    "eco": "UP000000625",   # very distant
}

# OMA is keyed by taxon id, except where only its own genome code resolves.
OMA_TEMPLATE = 559292
OMA_TARGETS = {"kla": 284590, "yli": 284591, "ani": 227321, "eco": "ECOLI"}

BASE = {"max_evalue": 1e-30, "min_align_len": 200, "min_identity": 40}
GRID = {
    "min_identity": [20, 25, 30, 35, 40, 45, 50, 60],
    "min_align_len": [50, 100, 150, 200, 300],
    "max_evalue": [1e-100, 1e-50, 1e-30, 1e-10, 1e-4],
}


def _stream(url: str) -> str:
    with urllib.request.urlopen(url, timeout=900) as fh:
        return fh.read().decode("utf-8")


def cmd_fetch(args: argparse.Namespace) -> int:
    """Download each proteome and relabel its headers with KEGG gene ids."""
    out = args.out / "proteomes"
    out.mkdir(parents=True, exist_ok=True)

    for code, upid in PROTEOMES.items():
        target = out / f"{code}.faa"
        if target.exists():
            print(f"{code}: present, skipping")
            continue

        mapping = {}
        rows = _stream(
            f"https://rest.uniprot.org/uniprotkb/stream?query=proteome:{upid}"
            f"&format=tsv&fields=accession,xref_kegg"
        ).splitlines()[1:]
        for line in rows:
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            for entry in filter(None, (p.strip() for p in parts[1].split(";"))):
                prefix, _, gene = entry.partition(":")
                # Exact code only: K-12 entries carry both eco: and ecj:, and
                # "whichever is commonest" would pick a different genome.
                if prefix == code:
                    mapping[parts[0]] = gene
                    break

        fasta = _stream(
            f"https://rest.uniprot.org/uniprotkb/stream?query=proteome:{upid}&format=fasta"
        )
        kept, dropped, lines, emit = 0, 0, [], False
        for line in fasta.splitlines():
            if line.startswith(">"):
                parts = line[1:].split("|")
                accession = parts[1] if len(parts) > 2 else parts[0]
                gene = mapping.get(accession)
                emit = gene is not None
                if emit:
                    lines.append(f">{gene}")
                    kept += 1
                else:
                    dropped += 1
            elif emit:
                lines.append(line)

        target.write_text("\n".join(lines) + "\n", encoding="utf-8")
        coverage = kept / (kept + dropped) if kept + dropped else 0.0
        print(f"{code}: {kept} kept, {dropped} unmapped (coverage {coverage:.1%})")
    return 0


def cmd_align(args: argparse.Namespace) -> int:
    """One bidirectional alignment per (template, target) pair, cached."""
    from raven_toolbox.reconstruction.homology import run_blast, run_diamond

    prot = args.out / "proteomes"
    prefix = "hits" if args.aligner == "blast" else f"hits_{args.aligner}"
    for code in PROTEOMES:
        if code == TEMPLATE:
            continue
        target = args.out / f"{prefix}_{code}.csv"
        if target.exists():
            print(f"{code}: present, skipping")
            continue
        started = time.perf_counter()
        if args.aligner == "blast":
            hits = run_blast(
                code, prot / f"{code}.faa", [TEMPLATE], [prot / f"{TEMPLATE}.faa"],
                evalue=args.collection_evalue, threads=args.threads,
            )
        else:
            hits = run_diamond(
                code, prot / f"{code}.faa", [TEMPLATE], [prot / f"{TEMPLATE}.faa"],
                threads=args.threads,
            )
        hits.to_csv(target, index=False)
        print(f"{code}: {len(hits)} hits in {time.perf_counter() - started:.0f}s")
    return 0


def _oma_truth(code: str, target_map: dict, source_map: dict, cache_dir: pathlib.Path):
    """Ortholog pairs from OMA, as KEGG gene names.

    OMA infers orthologs independently of BLAST, so it answers the objection that
    KEGG's annotations partly reflect the method under test. Its identifiers are
    UniProt, hence the mapping through both sides.
    """
    import re

    cache = cache_dir / f"oma_pairs_{code}.json"
    if cache.is_file():
        raw = json.loads(cache.read_text(encoding="utf-8"))
    else:
        raw, url = [], (
            f"https://omabrowser.org/api/pairs/{OMA_TARGETS[code]}/{OMA_TEMPLATE}/?page_size=1000"
        )
        while url:
            request = urllib.request.Request(
                url, headers={"Accept": "application/json", "User-Agent": "raven-toolbox/1.0"}
            )
            with urllib.request.urlopen(request, timeout=300) as fh:
                page, headers = json.load(fh), fh.headers
            raw += [
                (e["entry_1"]["canonicalid"], e["entry_2"]["canonicalid"]) for e in page
            ]
            match = re.search(r'<([^>]+)>;\s*rel="next"', headers.get("Link", ""))
            url = match.group(1) if match else None
            time.sleep(0.5)
        cache.write_text(json.dumps(raw), encoding="utf-8")

    truth, unmapped = set(), 0
    for target_id, template_id in raw:
        tg, sg = target_map.get(target_id), source_map.get(template_id)
        if tg and sg:
            truth.add((sg, tg))
        else:
            unmapped += 1
    return truth, len(raw), unmapped


def _load_gene_ko(path: pathlib.Path, wanted: set[str]) -> dict:
    opener = gzip.open if path.suffix == ".gz" else open
    per: dict[str, dict[str, set]] = collections.defaultdict(lambda: collections.defaultdict(set))
    with opener(path, "rt", encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            org, gene, ko = line.rstrip("\n").split("\t")
            if org in wanted:
                per[org][gene].add(ko)
    return per


def _uniprot_keys(code: str) -> dict[str, str]:
    """Every UniProt identifier for this organism -> its KEGG gene name.

    Both the accession and the entry name are included, because OMA uses one or
    the other depending on the genome.
    """
    url = (
        f"https://rest.uniprot.org/uniprotkb/stream?query=proteome:{PROTEOMES[code]}"
        f"&format=tsv&fields=accession,id,xref_kegg"
    )
    mapping = {}
    for line in _stream(url).splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        for entry in filter(None, (p.strip() for p in parts[2].split(";"))):
            prefix, _, gene = entry.partition(":")
            if prefix == code:
                mapping[parts[0]] = gene
                mapping[parts[1]] = gene
                break
    return mapping


def cmd_score(args: argparse.Namespace) -> int:
    """Score every threshold combination against the chosen reference."""
    import pandas as pd

    gene_ko = _load_gene_ko(args.gene_ko, set(PROTEOMES))
    source_ko = gene_ko[TEMPLATE]
    source_keys = _uniprot_keys(TEMPLATE) if args.reference == "oma" else {}
    prefix = "hits" if args.aligner == "blast" else f"hits_{args.aligner}"

    results = {}
    for code in PROTEOMES:
        if code == TEMPLATE:
            continue
        path = args.out / f"{prefix}_{code}.csv"
        if not path.is_file():
            print(f"{code}: no hit table, run 'align' first")
            continue

        hits = pd.read_csv(path)
        target_ko = gene_ko[code]

        by_ko: dict[str, list[str]] = {}
        for gene, kos in source_ko.items():
            for ko in kos:
                by_ko.setdefault(ko, []).append(gene)
        truth = {
            (source_gene, gene)
            for gene, kos in target_ko.items()
            for ko in kos
            for source_gene in by_ko.get(ko, ())
        }
        judged_source, judged_target = set(source_ko), set(target_ko)

        if args.reference == "oma":
            target_keys = _uniprot_keys(code)
            truth, n_pairs, unmapped = _oma_truth(
                code, target_keys, source_keys, args.out
            )
            # OMA covers whole proteomes, so every gene with a KEGG name counts.
            judged_source = set(source_keys.values())
            judged_target = set(target_keys.values())
            print(f"{code}: {n_pairs} OMA pairs, {len(truth)} usable "
                  f"({unmapped} unmapped on one side)")

        combos = [dict(BASE)]
        for name, values in GRID.items():
            for value in values:
                combo = {**BASE, name: value}
                if combo not in combos:
                    combos.append(combo)

        rows = []
        for combo in combos:
            kept = hits[
                (hits.evalue <= combo["max_evalue"])
                & (hits.align_len >= combo["min_align_len"])
                & (hits.identity >= combo["min_identity"])
            ]
            called, unjudgeable = set(), 0
            for a, b in zip(kept.from_gene, kept.to_gene, strict=True):
                if a in judged_source and b in judged_target:
                    called.add((a, b))
                elif b in judged_source and a in judged_target:
                    called.add((b, a))
                else:
                    # The reference says nothing about one of these genes:
                    # unjudgeable, not wrong.
                    unjudgeable += 1
            tp, fp = len(called & truth), len(called - truth)
            fn = len(truth - called)
            precision = tp / (tp + fp) if tp + fp else 0.0
            recall = tp / (tp + fn) if tp + fn else 0.0
            # beta < 1 weights precision above recall. For reconstruction that is
            # the honest weighting: a missed reaction can be recovered by
            # gap-filling, while a wrongly transferred one pollutes the model and
            # its gene associations silently.
            b2 = args.beta ** 2
            score = (
                (1 + b2) * precision * recall / (b2 * precision + recall)
                if precision and recall else 0.0
            )
            rows.append({**combo, "judged": len(called), "unjudgeable": unjudgeable,
                         "tp": tp, "fp": fp, "fn": fn, "precision": round(precision, 4),
                         "recall": round(recall, 4), "beta": args.beta,
                         "fbeta": round(score, 4)})

        results[code] = {"truth_pairs": len(truth), "n_hits": len(hits), "rows": rows}
        base = next(r for r in rows if all(r[k] == v for k, v in BASE.items()))
        best = max(rows, key=lambda r: r["fbeta"])
        label = f"F{args.beta:g}"
        print(f"{code}: default {label}={base['fbeta']:.3f} (P={base['precision']:.3f} "
              f"R={base['recall']:.3f})  best {label}={best['fbeta']:.3f} at ide="
              f"{best['min_identity']} len={best['min_align_len']}")

    (args.out / "hitlevel_results.json").write_text(
        json.dumps(results, indent=2), encoding="utf-8"
    )
    print(f"wrote {args.out / 'hitlevel_results.json'}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=pathlib.Path, default=pathlib.Path("work"),
                        help="working directory for proteomes, hit tables and results")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("fetch", help="download proteomes, relabelled with KEGG gene ids")

    align = sub.add_parser("align", help="one bidirectional alignment per pair")
    align.add_argument("--aligner", choices=("blast", "diamond"), default="blast",
                       help="DIAMOND is 10-20x faster and reaches the same optimum")
    align.add_argument("--threads", type=int, default=4)
    align.add_argument("--collection-evalue", type=float, default=1e-4,
                       help="hit-collection cutoff (RAVEN's getBlast hardcodes 1e-4)")

    score = sub.add_parser("score", help="score thresholds against KO sharing")
    score.add_argument("--gene-ko", type=pathlib.Path, required=True,
                       help="KEGG organism_gene_ko table (.tsv or .tsv.gz)")
    score.add_argument("--aligner", choices=("blast", "diamond"), default="blast")
    score.add_argument("--reference", choices=("kegg", "oma"), default="kegg",
                       help="which source decides whether two genes are counterparts")
    score.add_argument("--beta", type=float, default=0.5,
                       help="F-beta weighting; below 1 favours precision, which is the "
                            "agreed weighting for reconstruction (default: 0.5)")

    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    return {"fetch": cmd_fetch, "align": cmd_align, "score": cmd_score}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
