"""Tier 2: the extracted model must not drift unnoticed.

Extraction is a MILP, so there is no single right answer to assert -- equally
good optima exist and a different set is not automatically a wrong set. What can
be asserted is that today's result matches the result that was last inspected
and accepted, recorded by ``scripts/parity/record_baseline.py``.

The baseline is currently seeded from raven-toolbox itself, so this is a
regression guard rather than a cross-language check; the ``source`` field says
which, and the test reports it on failure so nobody mistakes one for the other.
When an extraction oracle is generated from MATLAB, the same comparison runs
against RAVEN's answer.

Exact equality is asserted rather than an overlap band because it was measured:
on this fixture GLPK and Gurobi return the same 13 reactions, and each is
identical across repeated runs. A difference therefore means something in this
package changed, not that the solver picked another optimum. Should that stop
being true on a larger fixture, the right move is a band with a measured floor
-- not a loosened threshold here.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from raven_toolbox.init import run_init
from raven_toolbox.io import read_yaml_model

pytestmark = pytest.mark.parity

BASELINE = (
    Path(__file__).resolve().parent.parent
    / "data"
    / "parity"
    / "baselines"
    / "init_smallyeast.json"
)


def _scores(model) -> dict[str, float]:
    """Must match scripts/parity/record_baseline.py, or the comparison is void."""
    return {
        rxn.id: (10.0 if i % 3 == 0 else (-5.0 if i % 3 == 1 else 1.0))
        for i, rxn in enumerate(model.reactions)
    }


@pytest.fixture(scope="module")
def baseline() -> dict:
    if not BASELINE.is_file():
        pytest.skip(
            f"no baseline at {BASELINE}. Record one with "
            f"python scripts/parity/record_baseline.py"
        )
    return json.loads(BASELINE.read_text(encoding="utf-8"))


def test_extraction_matches_the_recorded_baseline(raven_root, baseline):
    fixture = raven_root.joinpath(*baseline["fixture"].split("/"))
    if not fixture.is_file():
        pytest.skip(f"fixture {baseline['fixture']} not in this RAVEN checkout")

    model = read_yaml_model(fixture)
    assert len(model.reactions) == baseline["model_reactions"], (
        "the fixture itself changed upstream, so the baseline describes a "
        "different model -- re-record it rather than adjusting this test"
    )

    kept = sorted(r.id for r in run_init(model, _scores(model)).model.reactions)
    expected = baseline["kept_reactions"]

    if kept == expected:
        return

    got, want = set(kept), set(expected)
    overlap = len(got & want) / len(got | want)
    pytest.fail(
        f"extraction drifted from the recorded baseline "
        f"(source: {baseline['source']}, recorded {baseline['recorded']['date']} "
        f"with {baseline['recorded']['solver']}).\n"
        f"  Jaccard: {overlap:.3f} ({len(got & want)} shared of {len(got | want)})\n"
        f"  added:   {sorted(got - want)}\n"
        f"  dropped: {sorted(want - got)}\n"
        f"If the change is intended, re-record with "
        f"python scripts/parity/record_baseline.py and say in the PR why it moved."
    )


def test_baseline_describes_a_real_decision(baseline):
    """A baseline that kept everything (or nothing) would assert nothing."""
    kept = len(baseline["kept_reactions"])
    total = baseline["model_reactions"]
    assert 0 < kept < total, (
        f"the baseline keeps {kept} of {total} reactions, so the extraction is "
        f"not choosing -- fix the fixture or the scores, not this assertion"
    )


def test_solvers_agree_on_the_extraction(raven_root, baseline):
    """The extracted set must not depend on which solver found it.

    This is the measurement the exact-equality assertion above rests on, kept
    as a test so it stays true. If a future model makes the solvers disagree,
    this fails first and says so -- which is the signal to move that fixture to
    an overlap band with a measured floor, rather than to loosen the assertion
    above and lose the ability to detect real drift.

    Solvers that are not installed, or that refuse the model (the size-limited
    licence bundled with pip-installed gurobipy will refuse a large one), are
    skipped individually rather than failing the test.
    """
    from cobra.util.solver import solvers as available

    fixture = raven_root.joinpath(*baseline["fixture"].split("/"))
    if not fixture.is_file():
        pytest.skip(f"fixture {baseline['fixture']} not in this RAVEN checkout")

    candidates = [name for name in ("glpk", "gurobi", "cplex") if name in available]
    results: dict[str, list[str]] = {}
    refused: dict[str, str] = {}

    for name in candidates:
        model = read_yaml_model(fixture)
        try:
            model.solver = name
            results[name] = sorted(
                r.id for r in run_init(model, _scores(model)).model.reactions
            )
        except Exception as exc:  # noqa: BLE001 - reported, not swallowed
            refused[name] = f"{type(exc).__name__}: {exc}"

    if len(results) < 2:
        pytest.skip(
            f"need two usable solvers to compare; ran {sorted(results)}, "
            f"refused {refused or 'none'}"
        )

    distinct = {tuple(v) for v in results.values()}
    assert len(distinct) == 1, (
        "solvers disagree on the extracted set, so exact equality is no longer "
        "the right assertion for this fixture:\n"
        + "\n".join(f"  {name}: {len(kept)} reactions" for name, kept in results.items())
        + f"\n  symmetric difference: "
        f"{sorted(set(results[candidates[0]]) ^ set(results[candidates[1]]))}"
    )
