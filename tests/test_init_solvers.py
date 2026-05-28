"""Cross-solver smoke tests for the (f)tINIT MILP path.

The clean-data calibration and robustness studies were run on Gurobi; the tractability
choices (big-M=100, MIP gap, time limits) and the Gurobi-specific param plumbing
(``opt.problem.Params.MIPGap``) only matter if those choices also work on the *other*
MILP backends real users have. These tests assert that each available MILP-capable
optlang interface produces the same reaction-set verdict as Gurobi on the toy models the
unit tests use — so a regression in solver portability fails CI instead of being found
months later on a user's machine.

Solvers tested: every MILP-capable cobra/optlang interface that imports in this env
(Gurobi, HiGHS via ``hybrid``, GLPK). Missing ones are skipped automatically. Genome-scale
behaviour is measured separately by ``scripts/analyze_init_solvers.py`` (manual benchmark).
"""
from __future__ import annotations

import importlib

import cobra
import pytest

from ravengem.init import ftinit, prep_init_model, run_ftinit, run_init
from ravengem.tasks import Task, check_tasks

# Detect which MILP-capable optlang interfaces are installed; skip the rest.
_INTERFACES = {"gurobi": "gurobi_interface", "hybrid": "hybrid_interface", "glpk": "glpk_interface"}
_AVAILABLE = [name for name, mod in _INTERFACES.items()
              if importlib.util.find_spec(f"optlang.{mod}") is not None]


@pytest.fixture(params=_AVAILABLE)
def solver(request):
    """One installed MILP solver per parameter value."""
    return request.param


# ----------------------------------------------------------------------- toy fixtures

def _met(mid, comp="c"):
    return cobra.Metabolite(mid, name=mid.split("_")[0], compartment=comp)


def _toy_init_model() -> cobra.Model:
    """EX_A → A → B → C → D (r1, r2 good; r3 bad). Same network as test_init.py."""
    def rxn(rid, lb, ub, mets):
        r = cobra.Reaction(rid, lower_bound=lb, upper_bound=ub)
        r.add_metabolites(mets)
        return r
    m = cobra.Model("toy")
    A, B, C, D = (_met(x) for x in ("A_c", "B_c", "C_c", "D_c"))
    m.add_metabolites([A, B, C, D])
    m.add_reactions([rxn("EX_A", -1000, 1000, {A: -1}),
                     rxn("r1", 0, 1000, {A: -1, B: 1}),
                     rxn("r2", 0, 1000, {B: -1, C: 1}),
                     rxn("r3", 0, 1000, {C: -1, D: 1})])
    return m


def _toy_ftinit_model() -> cobra.Model:
    """Small flux-consistent network for ftINIT: A→B, B→C, parallel A→C (negative-score)."""
    def rxn(rid, lb, ub, mets):
        r = cobra.Reaction(rid, lower_bound=lb, upper_bound=ub)
        r.add_metabolites(mets)
        return r
    m = cobra.Model("ftoy")
    A, B, C = (_met(x) for x in ("A_c", "B_c", "C_c"))
    m.add_metabolites([A, B, C])
    m.add_reactions([rxn("EX_A", -1000, 0, {A: -1}),
                     rxn("EX_C", 0, 1000, {C: -1}),
                     rxn("r1", 0, 1000, {A: -1, B: 1}),
                     rxn("r2", 0, 1000, {B: -1, C: 1}),
                     rxn("rbad", 0, 1000, {A: -1, C: 1})])
    return m


# --------------------------------------------------------------------- tests

def test_run_init_same_verdict(solver):
    """tINIT MILP on a small network drops the negative-score reaction with any solver."""
    m = _toy_init_model()
    m.solver = solver
    res = run_init(m, {"r1": 1.0, "r2": 1.0, "r3": -1.0}, prod_weight=0.0, allow_excretion=True)
    assert "r3" in res.deleted_reactions
    assert sorted(set(r.id for r in res.model.reactions)) == ["EX_A", "r1", "r2"]


def test_run_ftinit_same_verdict(solver):
    """ftINIT MILP picks the same on-set across solvers on a small network."""
    m = _toy_ftinit_model()
    m.solver = solver
    res = run_ftinit(m, {"r1": 1.0, "r2": 1.0, "rbad": -1.0}, allow_excretion=True)
    assert "rbad" not in res.on_reactions
    assert {"r1", "r2"}.issubset(res.on_reactions)


def test_check_tasks_works_per_solver(solver):
    """check_tasks (one slim_optimize per task) works with each solver."""
    m = _toy_ftinit_model()
    m.solver = solver
    task = Task(id="make_c", inputs=[("A[c]", 0.0, 1000.0)], outputs=[("C[c]", 1.0, 1.0)])
    results = check_tasks(m, [task])
    assert results[0].passed


def test_ftinit_pipeline_with_tasks(solver):
    """The full ftinit() pipeline (prep + staged MILP + gap-fill) runs with each solver."""
    m = _toy_ftinit_model()
    m.solver = solver
    task = Task(id="make_c", inputs=[("A[c]", 0.0, 1000.0)], outputs=[("C[c]", 1.0, 1.0)])
    prep = prep_init_model(m, [task])
    out = ftinit(prep, {"r1": 1.0, "r2": 1.0, "rbad": -1.0}, series="1+1")
    # Functional: the target task remains satisfiable in the extracted model.
    assert check_tasks(out, [task])[0].passed


def test_solver_param_plumbing(solver):
    """mip_gap / time_limit reach the solver without raising (graceful per backend)."""
    m = _toy_ftinit_model()
    m.solver = solver
    # Tight time limit + loose gap on a trivial problem; just verify the call returns.
    res = run_ftinit(m, {"r1": 1.0, "rbad": -1.0}, allow_excretion=True,
                     mip_gap=0.05, time_limit=60)
    assert res.objective is not None
