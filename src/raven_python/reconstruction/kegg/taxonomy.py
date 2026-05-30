"""Parse the KEGG ``taxonomy`` file into per-organism category lineages.

Ports the file-reading half of RAVEN ``getPhylDist`` (the distance-matrix half is
step 3b.5). The ``taxonomy`` file is an indented tree: ``#``-prefixed lines name a
category, the number of leading ``#`` giving its depth; organism lines are
tab-separated ``T-number<tab>org_id<tab>name<tab>...``. Each organism inherits the
stack of categories above it, the first of which is its domain (``Prokaryotes`` /
``Eukaryotes``).

Used by 3b.3 to split genes into the prok/euk HMM libraries, and (later) by 3b.5
for phylogenetic distances.
"""
from __future__ import annotations

import warnings
from pathlib import Path


def parse_taxonomy(path: str | Path) -> dict[str, list[str]]:
    """Return ``{organism_id: [category, ...]}`` from outermost to innermost."""
    org_categories: dict[str, list[str]] = {}
    stack: list[str] = []
    skipped_level_warned = False
    with open(path, encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("#"):
                depth = len(line) - len(line.lstrip("#"))
                name = line[depth:].strip()
                if depth - 1 > len(stack):
                    # Depth-skip (e.g. ## then ####): the original `stack[:depth-1]`
                    # silently produced a too-short lineage. Pad with explicit
                    # blanks so downstream slices stay aligned; warn once.
                    if not skipped_level_warned:
                        warnings.warn(
                            f"{path}: taxonomy depth skips a level near line {line_no} "
                            f"({'#' * depth} {name!r} appeared with stack {stack!r}); "
                            "padding the missing levels with '' (later occurrences silenced).",
                            stacklevel=2,
                        )
                        skipped_level_warned = True
                    stack = stack + [""] * (depth - 1 - len(stack))
                else:
                    stack = stack[: depth - 1]
                stack.append(name)
            else:
                fields = line.split("\t") if "\t" in line else line.split()
                if len(fields) < 2:
                    continue
                org_categories[fields[1].strip()] = list(stack)
    return org_categories


def organism_domains(path: str | Path) -> dict[str, str]:
    """Return ``{organism_id: domain}`` (the outermost category)."""
    return {org: cats[0] for org, cats in parse_taxonomy(path).items() if cats}


def organisms_in_domain(path: str | Path, domain: str) -> set[str]:
    """Organism ids whose outermost category matches ``domain`` (case-insensitive).

    Accepts a prefix, so ``"prok"`` matches ``"Prokaryotes"``.
    """
    needle = domain.lower()
    return {
        org
        for org, dom in organism_domains(path).items()
        if dom.lower().startswith(needle) or needle.startswith(dom.lower())
    }
