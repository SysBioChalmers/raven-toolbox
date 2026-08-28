"""Fail unless the Gurobi licence in effect will accept a genome-scale model.

The licence bundled with pip-installed gurobipy stops at 2000 variables and
2000 constraints -- ample for the small fixtures, useless for a real GEM. A
nightly job that silently fell back to it would skip the genome-scale checks
and still report success, which is exactly the kind of empty green this harness
exists to prevent. So the licence is proven before the tests run.
"""
from __future__ import annotations

import sys

SIZE = 20_000  # comfortably past the limited licence, comfortably under a GEM


def main() -> int:
    try:
        import gurobipy as gp
    except ImportError as exc:
        print(f"gurobipy is not installed: {exc}", file=sys.stderr)
        return 1

    try:
        with gp.Env() as env, gp.Model(env=env) as model:
            model.setParam("OutputFlag", 0)
            variables = model.addVars(SIZE, vtype=gp.GRB.BINARY)
            model.addConstr(variables.sum() >= 1)
            model.optimize()
            print(f"licence OK: solved a {model.NumVars}-variable MIP")
    except Exception as exc:  # noqa: BLE001 - the message is the whole point
        print(
            f"Gurobi refused a {SIZE}-variable problem: {exc}\n"
            f"The size-limited licence bundled with gurobipy cannot run the "
            f"genome-scale checks; a WLS licence must be in effect.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
