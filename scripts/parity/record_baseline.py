"""Record the set-level baseline the parity suite compares against.

Set-level checks (tier 2 in ``tests/parity/README.md``) cannot assert a single
right answer: extraction is a MILP and equally-good optima exist. What they can
do is detect *drift* -- today's result differing from the result that was
inspected and accepted.

This script records that accepted result. Run it when a change is expected to
move the extraction, inspect the diff the test reports, and commit the new
baseline with a note in the PR saying why it moved.

    python scripts/parity/record_baseline.py            # uses $RAVEN_ROOT
    python scripts/parity/record_baseline.py --raven-root /path/to/RAVEN

The baseline records which implementation produced it. It is currently seeded
from raven-toolbox itself, which makes the check a regression guard rather than
a cross-language one; when ``scripts/parity/generate_oracles.m`` grows an
extraction oracle, the same comparison runs against MATLAB RAVEN's answer and
the ``source`` field changes to say so.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import platform
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
BASELINE_DIR = REPO_ROOT / "tests" / "data" / "parity" / "baselines"

# The fixture is RAVEN's own small yeast model, read from a RAVEN checkout
# rather than vendored: RAVEN is GPL and this package is MIT.
FIXTURE = ("tutorial", "smallYeast.yml")


def reaction_scores(model) -> dict[str, float]:
    """Deterministic scores that force a real decision.

    A third of the reactions are clearly worth keeping, a third clearly not,
    and a third marginal -- so the extraction has to choose rather than keep
    everything.
    """
    return {
        rxn.id: (10.0 if i % 3 == 0 else (-5.0 if i % 3 == 1 else 1.0))
        for i, rxn in enumerate(model.reactions)
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raven-root",
        default=os.environ.get("RAVEN_ROOT"),
        help="a MATLAB RAVEN checkout (default: $RAVEN_ROOT)",
    )
    args = parser.parse_args(argv)

    if not args.raven_root:
        parser.error(
            "no RAVEN checkout given. Clone RAVEN and pass --raven-root, or set "
            "RAVEN_ROOT: git clone --depth 1 -b develop3 "
            "https://github.com/SysBioChalmers/RAVEN"
        )

    import cobra

    import raven_toolbox
    from raven_toolbox.init import run_init
    from raven_toolbox.io import read_yaml_model

    fixture = Path(args.raven_root).joinpath(*FIXTURE)
    if not fixture.is_file():
        parser.error(f"fixture not found: {fixture}")

    model = read_yaml_model(fixture)
    scores = reaction_scores(model)
    result = run_init(model, scores)
    kept = sorted(r.id for r in result.model.reactions)

    solver = cobra.Configuration().solver.__name__.split(".")[-1]
    baseline = {
        "source": "raven-toolbox",
        "description": (
            "Reactions kept by run_init on RAVEN's tutorial/smallYeast.yml with "
            "the scores from scripts/parity/record_baseline.py. Seeded from "
            "raven-toolbox, so this detects drift in this package rather than "
            "disagreement with MATLAB RAVEN."
        ),
        "fixture": "/".join(FIXTURE),
        "model_reactions": len(model.reactions),
        "kept_reactions": kept,
        "recorded": {
            "date": dt.date.today().isoformat(),
            "raven_toolbox": getattr(raven_toolbox, "__version__", "unknown"),
            "cobra": cobra.__version__,
            "solver": solver,
            "python": platform.python_version(),
        },
    }

    BASELINE_DIR.mkdir(parents=True, exist_ok=True)
    target = BASELINE_DIR / "init_smallyeast.json"
    target.write_text(json.dumps(baseline, indent=2) + "\n", encoding="utf-8")
    print(f"recorded {len(kept)}/{len(model.reactions)} reactions to {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
