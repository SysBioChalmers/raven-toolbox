#!/usr/bin/env python
"""Curated-GEM arm of the homology cut-off study -- and the test that retires it.

Reconstructs H. polymorpha from yeast-GEM and rhto-GEM and compares against the
curated hanpo-GEM. Kept in the repository not because the comparison is usable
but because the check that showed it is *not* usable is worth being able to
repeat -- on this pair, and on any other model somebody proposes as ground truth.

    python scripts/homology_cutoff_curated.py --hanpo-gem /path/to/hanpo-GEM

The curated model's draft was built by getModelFromHomology at
(1e-30, 150, 35) -- see its code/reconstructionProtocol.m. The ladder therefore
walks through those settings and reports what fraction of newly admitted
reactions land in the curated model. A cliff at the build settings means the
reference is echoing its own construction; a smooth decline would mean the
agreement is real.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

BUILD = {"max_evalue": 1e-30, "min_align_len": 150, "min_identity": 35}
DEFAULTS = {"max_evalue": 1e-30, "min_align_len": 200, "min_identity": 40}

LADDER = [
    {"max_evalue": 1e-30, "min_align_len": 200, "min_identity": 50},
    {"max_evalue": 1e-30, "min_align_len": 200, "min_identity": 45},
    DEFAULTS,
    {"max_evalue": 1e-30, "min_align_len": 175, "min_identity": 37},
    BUILD,
    {"max_evalue": 1e-30, "min_align_len": 140, "min_identity": 33},
    {"max_evalue": 1e-30, "min_align_len": 125, "min_identity": 30},
    {"max_evalue": 1e-30, "min_align_len": 100, "min_identity": 25},
]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hanpo-gem", type=pathlib.Path, required=True,
                        help="a hanpo-GEM checkout (supplies proteomes, templates and truth)")
    parser.add_argument("--hits", type=pathlib.Path,
                        help="cached hit table; computed and saved here if absent")
    parser.add_argument("--out", type=pathlib.Path, default=pathlib.Path("work"))
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    import cobra
    import pandas as pd

    from raven_toolbox.io import read_yaml_model
    from raven_toolbox.reconstruction.homology import get_model_from_homology, run_blast

    genomes = args.hanpo_gem / "data" / "genomes"
    templates = args.hanpo_gem / "data" / "templateModels"

    hits_path = args.hits or args.out / "hits_hanpo.csv"
    if hits_path.is_file():
        hits = pd.read_csv(hits_path)
    else:
        hits = run_blast(
            "hanpo", genomes / "hanpo.faa", ["sce", "rhto"],
            [genomes / "sce.faa", genomes / "rhto.faa"], evalue=1e-4, threads=4,
        )
        hits.to_csv(hits_path, index=False)
    print(f"{len(hits)} hits")

    sce = cobra.io.read_sbml_model(str(templates / "yeastGEM.xml"))
    sce.id = "sce"
    rhto = cobra.io.read_sbml_model(str(templates / "rhto.xml"))
    rhto.id = "rhto"
    curated = {r.id for r in read_yaml_model(args.hanpo_gem / "model" / "hanpo-GEM.yml").reactions}

    reachable = curated & ({r.id for r in sce.reactions} | {r.id for r in rhto.reactions})
    print(f"curated {len(curated)} reactions, {len(reachable)} reachable "
          f"-> ceiling {len(reachable) / len(curated):.3f}")

    rows, previous = [], None
    for combo in LADDER:
        draft = get_model_from_homology(
            [sce, rhto], hits, "hanpo",
            preferred_order=["sce", "rhto"], strictness=1, **combo,
        ).model
        current = {r.id for r in draft.reactions}
        row = {**combo, "n_draft": len(current), "in_curated": len(current & curated)}
        if previous is not None:
            gained = current - previous
            row["gained"] = len(gained)
            row["gained_in_curated"] = len(gained & curated)
            row["hit_rate"] = len(gained & curated) / len(gained) if gained else None
        rows.append(row)
        previous = current

        note = "  <- build settings" if combo == BUILD else ""
        rate = f"{row['hit_rate']:.3f}" if row.get("hit_rate") is not None else "baseline"
        print(f"len={combo['min_align_len']:<4d} ide={combo['min_identity']:<3d} "
              f"draft={row['n_draft']:<5d} newly-admitted-in-curated={rate}{note}")

    (args.out / "contamination.json").write_text(
        json.dumps({"curated": len(curated), "ladder": rows}, indent=2), encoding="utf-8"
    )
    print(f"wrote {args.out / 'contamination.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
