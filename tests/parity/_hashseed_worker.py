"""Run one gap-fill and print a digest of the result. Not a test module.

Invoked as a subprocess by ``test_determinism.py`` under different
``PYTHONHASHSEED`` values. Python randomises string hashing per process, so a
set iterated to build constraint rows, or a dict deciding a tie-break, produces
a stable answer *within* one process and a different one in the next. An
in-process repetition cannot see that; two processes with different seeds can.

Prints one JSON object on the last line so the caller can parse it out of
whatever the solver has written to stdout.
"""
from __future__ import annotations

import hashlib
import json
import sys


def build_models():
    import cobra

    def model_with(model_id, reactions, objective=None):
        model = cobra.Model(model_id)
        mets = {
            mid: cobra.Metabolite(mid, compartment="c")
            for mid in ("a", "b", "c", "d")
        }
        model.add_metabolites(list(mets.values()))
        for rid, stoich, lb, ub in reactions:
            rxn = cobra.Reaction(rid, lower_bound=lb, upper_bound=ub)
            rxn.add_metabolites({mets[m]: v for m, v in stoich.items()})
            model.add_reactions([rxn])
        if objective:
            model.objective = objective
        return model

    draft = model_with(
        "draft",
        [
            ("EX_a", {"a": 1}, -10.0, 1000.0),
            ("R1", {"a": -1, "b": 1}, 0.0, 1000.0),
            ("EX_d", {"d": -1}, 0.0, 1000.0),
        ],
        objective="EX_d",
    )
    # Two symmetric routes plus a two-step alternative: nothing in the problem
    # distinguishes them, so any hash-order leak shows up as a different pick.
    template = model_with(
        "template",
        [
            ("T_bd_1", {"b": -1, "d": 1}, 0.0, 1000.0),
            ("T_bd_2", {"b": -1, "d": 1}, 0.0, 1000.0),
            ("T_bc", {"b": -1, "c": 1}, 0.0, 1000.0),
            ("T_cd", {"c": -1, "d": 1}, 0.0, 1000.0),
        ],
    )
    return draft, template


def main() -> int:
    import cobra

    from raven_toolbox.gapfilling import connect_blocked_reactions

    # Single-process: the gap-filler runs FVA, which otherwise starts a worker
    # pool per invocation. On a four-reaction model that costs far more than it
    # saves, and three of these workers running pools at once makes the test
    # slower than everything else in the suite combined.
    cobra.Configuration().processes = 1

    draft, template = build_models()
    result = connect_blocked_reactions(draft, [template])
    payload = {
        "added": sorted(result.added_reactions),
        "connected": sorted(result.newly_connected),
    }
    payload["digest"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode("utf-8")
    ).hexdigest()[:16]
    print("HASHSEED_RESULT " + json.dumps(payload))
    return 0


if __name__ == "__main__":
    sys.exit(main())
