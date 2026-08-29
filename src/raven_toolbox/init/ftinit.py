"""The ftINIT MILP — the faster staged variant of INIT.

ftINIT keeps tINIT's objective — pick the reaction subset best matching expression
scores while staying flux-consistent — but with a cheaper MILP encoding that is the
reason it is *fast*: a **positive-score reaction needs no binary**. Because the
objective *maximises* ``Σ score·y`` with ``score > 0``, the optimiser pushes its
continuous indicator ``y ∈ [0,1]`` to 1, and the gate ``net_flux ≥ force_on·y`` only
lets ``y`` reach 1 if the reaction can actually carry flux. Only *negative*-score
reactions need a true ``{0,1}`` binary (their indicator would otherwise sit at 0 for
free). This roughly halves the integer count — the dominant MILP cost.

Reaction categories (RAVEN's six), by score sign × reversibility:

* **score 0** — left in the model, *not* in the problem: a free flux variable that can
  carry flux for connectivity but is neither scored nor removable.
* **positive, irreversible** — continuous ``y∈[0,1]``; ``v ≥ force_on·y``. No binary.
* **positive, reversible** — split ``v = v⁺ − v⁻``; continuous ``y``; a single
  direction binary keeps one of ``v⁺/v⁻`` at 0 (no fwd/back loop faking "on");
  ``v⁺+v⁻ ≥ force_on·y``.
* **negative, irreversible** — binary ``x∈{0,1}``; ``v ≤ ub·x``.
* **negative, reversible** — split; binary ``x``; ``v⁺+v⁻ ≤ cap·x``.
* **essential** — forced on (``v ≥ force_on_ess``); no indicator. Assumed already
  oriented irreversible in its forced direction (``prepINITModel`` does this).

Objective: **maximise** ``Σ score·indicator``. Unlike classic INIT
(:func:`raven_toolbox.init.run_init`), ftINIT does **not** reward production of every
metabolite — ``prod_weight`` applies only to metabolomics-detected metabolites (not
yet implemented; passing a non-empty ``metabolomics`` argument raises
``NotImplementedError``). Connectivity comes solely from the flux gates plus any
essential reactions. ``allow_excretion`` relaxes ``S·v = 0`` to ``≥ 0``; ``rem_pos_rev``
drops positive reversible reactions from the problem (used in the staging schedule).

Needs a MILP solver (cobra's configured optlang solver; only Gurobi is fully viable at
genome scale — see ``docs/init_solver_benchmark.md``). Magic numbers
(``force_on``/``force_on_ess`` = 0.1, ``big_m`` = 100) are exposed and scale-dependent;
calibration tables are in ``docs/init_param_calibration.md``. ``big_m`` caps a *scored*
reaction's flux in its on/off (direction) constraint — using a fixed 100 rather than
the reaction's ±1000 bound keeps the LP relaxation tight (what makes the genome-scale
MILP tractable). Free / essential reactions keep their real bounds.

⚠️ **Loops.** The MILP has *no* loopless constraint: an internal
thermodynamically-infeasible cycle is flux-consistent (``S·v = 0``), so if its
reactions carry positive net score the optimiser will "include" them with no real
exchange flux. RAVEN tolerates this — loop-free models come from the staged pipeline
+ exchange handling, and at genome scale real exchange reactions make such cycles not
score-optimal. A loopless option could be layered on later if needed.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field

import cobra
from cobra.exceptions import OptimizationError
from optlang.interface import Variable  # untyped (resolves to Any); used for container hints
from optlang.symbolics import Real, add, mul

from raven_toolbox.init.genes import remove_low_score_genes
from raven_toolbox.init.merge import group_rxn_scores
from raven_toolbox.init.steps import get_init_steps
from raven_toolbox.init.taskfill import fill_tasks

_FORCE_ON = 0.1  # min flux for a reaction to count as "on" (RAVEN forceOnLim)
_BIG_M = 100.0   # indicator/direction big-M cap on a *scored* reaction's flux (RAVEN's 100)
_STRICT_ABS_GAP = 0.05  # absolute gap for the opt-in strict mode (below the 0.1 score granularity)
_EXTRACT_SEED = 1234  # RAVEN optimizeProb Seed; a module constant so a determinism probe
                      # can vary it (see the ftinit-determinism study on raven-docs) without
                      # editing code


def _dbg(msg: str) -> None:
    """Print a diagnostic line to stderr when FTINIT_DEBUG is set (off by default)."""
    import os
    import sys
    if os.environ.get("FTINIT_DEBUG"):
        print(msg, file=sys.stderr, flush=True)


@dataclass
class FtInitResult:
    """Result of :func:`run_ftinit`."""

    model: cobra.Model
    kept_reactions: list[str]
    deleted_reactions: list[str]
    fluxes: dict[str, float]
    objective: float
    on_reactions: set[str] = field(default_factory=set)  # scored reactions turned on (indicator)
    achieved_gap: float | None = None  # solver's final relative MIP gap (RAVEN's multi-run driver)


def run_ftinit(
    model: cobra.Model,
    rxn_scores: Mapping[str, float] | None = None,
    *,
    essential_rxns: Iterable[str] | None = None,
    essential_directions: Mapping[str, int] | None = None,
    essential_force: Mapping[str, float] | None = None,
    allow_excretion: bool = False,
    rem_pos_rev: bool = False,
    ignore_mets: Iterable[str] = (),
    force_on: float = _FORCE_ON,
    force_on_ess: float = _FORCE_ON,
    big_m: float = _BIG_M,
    mip_gap: float | None = None,
    mip_gap_abs: float | None = None,
    time_limit: float | None = None,
    strict_abs_gap: float | None = None,
    canonical: bool = False,
) -> FtInitResult:
    """Run the single-step ftINIT MILP and return the extracted model.

    ``rxn_scores`` maps reaction id → score (default 0 → reaction left free in the
    model, not scored or removable). ``essential_rxns`` are forced to carry flux
    (≥ ``force_on_ess``, overridable per reaction via ``essential_force``);
    ``essential_directions`` maps an essential reaction id to ``+1`` (forward) or
    ``-1`` (reverse) for the forced direction (default forward). ``ignore_mets`` are
    metabolite **names** whose mass balance is dropped (RAVEN's per-step "simple
    metabolite" removal, e.g. H2O/H+). See the module docstring for the formulation.
    This is the single-step variant; the staged schedule
    (:func:`raven_toolbox.init.ftinit`) calls it per step.

    ``mip_gap`` / ``time_limit``: the default ``None`` uses the solver's own
    defaults (Gurobi: MIPGap≈1e-4, no time cap), which is fine for a single step
    — solve time here is dominated by model construction, so a tight gap is
    nearly free (measured on genome-scale Human-GEM, see the `INIT parameter
    calibration study
    <https://github.com/edkerk/raven-docs/blob/main/docs/parameter-tuning/studies/init-param-calibration.md>`_
    on raven-docs). For the genome-scale staged
    schedule (:func:`raven_toolbox.init.ftinit`), an essential-forced step can
    run away without a cap (one severely-degraded case ran >75 min unbounded);
    set ``time_limit`` ≈120-600 s/step, and loosen ``mip_gap`` to ``0.01`` (or
    ``0.005``) for a ~37% speedup at ≥0.99 Jaccard vs. the tight-gap model.

    ``strict_abs_gap``, when set, proves the solve to this fixed absolute objective
    gap instead of a relative ``mip_gap``/``mip_gap_abs`` — see the ``strict_gap``
    parameter of :func:`ftinit` for when to use it.

    ``canonical`` (opt-in, default off to preserve RAVEN parity) resolves the MILP's
    degeneracy deterministically: after the score optimum is found, a lexicographic
    phase 2 holds the objective and minimises the id-ordered count of "on" reactions, so
    the kept set is the unique sparsest optimum, independent of solver seed/version,
    instead of an arbitrary tie-break. See :func:`_canonicalize`.
    """
    scores = dict(rxn_scores or {})
    essential = set(essential_rxns or [])
    directions = dict(essential_directions or {})
    essential_force = dict(essential_force or {})
    ignore_met_names = set(ignore_mets)
    prob = model.problem
    opt = prob.Model()

    variables: list = []
    constraints: list = []
    flux_terms: dict[str, list[tuple[Variable, float]]] = {}  # rxn id -> [(var, sign)]
    indicators: dict[str, tuple[Variable, float]] = {}  # rxn id -> (indicator var, score)
    free_or_essential: set[str] = set()               # kept regardless of an indicator

    def add_constraint(expr, **kw):
        constraints.append(prob.Constraint(expr, **kw))

    for rxn in model.reactions:
        rid = rxn.id
        lb, ub = rxn.lower_bound, rxn.upper_bound
        score = float(scores.get(rid, 0.0))
        if rem_pos_rev and score > 0 and lb < 0 < ub:
            score = 0.0  # staging step 1: positive reversibles dropped from the problem

        if rid in essential:
            # Forced to carry flux in its forced direction (default forward); respect a
            # stricter native bound if the model already forces more flux. The forced
            # magnitude may be set per reaction (RAVEN's min(0.99·|prev flux|, 0.1), so
            # a reaction is never forced above what it carried before).
            force = essential_force.get(rid, force_on_ess) if essential_force else force_on_ess
            if directions.get(rid, 1) >= 0:
                forced = min(force, ub)  # clamp to capacity so we never make lb > ub
                v = prob.Variable(f"v_{rid}", lb=max(forced, lb, 0.0), ub=ub)
            else:  # reverse: flux ≤ -force
                forced = min(force, -lb)
                v = prob.Variable(f"v_{rid}", lb=lb, ub=min(-forced, ub))
            variables.append(v)
            flux_terms[rid] = [(v, 1.0)]
            free_or_essential.add(rid)
            continue

        if score == 0.0:  # free: carries flux for connectivity, not scored/removable
            v = prob.Variable(f"v_{rid}", lb=lb, ub=ub)
            variables.append(v)
            flux_terms[rid] = [(v, 1.0)]
            free_or_essential.add(rid)
            continue

        reversible = lb < 0 < ub
        if reversible:
            vp = prob.Variable(f"vp_{rid}", lb=0.0, ub=ub)
            vn = prob.Variable(f"vn_{rid}", lb=0.0, ub=-lb)
            variables += [vp, vn]
            flux_terms[rid] = [(vp, 1.0), (vn, -1.0)]
            total = vp + vn  # |flux| (one of vp/vn pinned to 0 below), used by the gates
        else:  # single-direction: keep the model's own [lb, ub] (incl. any forced lb>0)
            v = prob.Variable(f"v_{rid}", lb=lb, ub=ub)
            variables.append(v)
            flux_terms[rid] = [(v, 1.0)]
            total = v if ub > 0 else -v  # magnitude for a single-direction reaction

        if score > 0:
            y = prob.Variable(f"y_{rid}", lb=0.0, ub=1.0)  # continuous indicator, no binary
            variables.append(y)
            indicators[rid] = (y, score)
            add_constraint(total - force_on * y, lb=0.0, name=f"on_{rid}")  # y=1 ⇒ |flux| ≥ force_on
            if reversible:  # one direction binary stops a fwd/back loop faking "on"
                b = prob.Variable(f"b_{rid}", type="binary")
                variables.append(b)
                add_constraint(vp - big_m * b, ub=0.0, name=f"dirp_{rid}")          # vp ≤ M·b
                add_constraint(vn + big_m * b, ub=big_m, name=f"dirn_{rid}")        # vn ≤ M·(1-b)
        else:  # score < 0
            x = prob.Variable(f"x_{rid}", type="binary")
            variables.append(x)
            indicators[rid] = (x, score)
            add_constraint(total - big_m * x, ub=0.0, name=f"off_{rid}")  # flux>0 ⇒ x=1

    # Steady state S·v {== 0 | >= 0}; ignored metabolites are left unbalanced.
    # NOTE (allow_excretion sign): we relax to S·v >= 0 (net production / excretion
    # allowed), matching the constraint's name and classic INIT (RAVEN runINIT). RAVEN's
    # *ftINIT* builder instead uses csense 'L' → S·v <= 0 (net consumption), which reads as
    # inconsistent with its own "allowExcretion" intent. We deliberately keep the
    # intent-faithful S·v >= 0: allow_excretion is unused in the default '1+1' schedule
    # (both steps have allow_met_secr=False) and only fires in '2+1'/'2+0' step 1, so this
    # does not affect the gene-essentiality pipeline. Revisit only if exact '2+1' parity
    # with RAVEN is required.
    # Build each metabolite's balance as a *flat* list of (coeff·sign)·var terms and sum
    # it with optlang.symbolics.add. Python's builtin sum re-canonicalises a growing
    # sympy expression at every step (O(n²)); for hub metabolites that appear in ~10³
    # reactions that is minutes per constraint. add() builds the sum in one pass.
    met_terms: dict = {m: [] for m in model.metabolites if m.name not in ignore_met_names}
    for rxn in model.reactions:
        terms = flux_terms[rxn.id]
        for met, coeff in rxn.metabolites.items():
            bucket = met_terms.get(met)
            if bucket is None:
                continue
            for var, sign in terms:
                bucket.append(mul([Real(coeff * sign), var]))
    for termlist in met_terms.values():
        if termlist:
            add_constraint(add(termlist), lb=0.0, ub=None if allow_excretion else 0.0)

    opt.add(variables + constraints)
    obj_expr = add([mul([Real(score), ind]) for ind, score in indicators.values()])
    opt.objective = prob.Objective(obj_expr, direction="max")
    try:  # Gurobi-specific; harmless if the backend differs. Match RAVEN's optimizeProb
        # defaults exactly, because the ftINIT MILP is highly degenerate and the chosen
        # incumbent (hence which reactions are kept) depends on these:
        #   * Threads=1 — RAVEN forces single-threaded solving; multi-threaded Gurobi
        #     picks among equal optima non-deterministically and can even report the MILP
        #     infeasible (RAVEN issue #607). This is the dominant reproducibility lever.
        #   * Presolve=2, FeasibilityTol/OptimalityTol/IntFeasTol=1e-9 — RAVEN's
        #     optimizeProb defaults; they steer which optimal vertex a degenerate MILP
        #     lands on and how binaries round at the 0.5 on/off cut.
        #   * Seed (``_EXTRACT_SEED``=1234) — fixed seed so tie-breaking is reproducible.
        opt.problem.Params.Threads = 1
        opt.problem.Params.Presolve = 2
        opt.problem.Params.FeasibilityTol = 1e-9
        opt.problem.Params.OptimalityTol = 1e-9
        opt.problem.Params.IntFeasTol = 1e-9
        opt.problem.Params.Seed = _EXTRACT_SEED
    except Exception:  # noqa: BLE001
        pass
    if strict_abs_gap is not None:
        # Strict mode: prove the optimum to a fixed *absolute* gap below the reaction-score
        # granularity (scores are nudged to |score| ≥ 0.1), so the kept-set objective is
        # proven optimal at any objective scale. A relative gap is meaningless on the
        # near-zero-objective final step; this also bypasses the loose escalation below.
        try:  # Gurobi-specific; harmless if the backend differs
            opt.problem.Params.MIPGap = 0.0
            opt.problem.Params.MIPGapAbs = strict_abs_gap
        except Exception:  # noqa: BLE001
            pass
    elif mip_gap is not None:
        try:  # Gurobi-specific; harmless if the backend differs
            opt.problem.Params.MIPGap = mip_gap
        except Exception:  # noqa: BLE001
            pass

    if strict_abs_gap is None and mip_gap_abs is not None:
        # RAVEN's multi-run gap strategy (ftINIT.m). The final staged step has a
        # near-zero objective (it mostly removes small negative-score reactions), so a
        # fixed *relative* gap becomes an almost-zero absolute gap the solver cannot
        # reach in the time limit — it stops with extra bypass reactions on, giving a
        # suboptimal model. A fixed *absolute* gap (Gurobi MIPGapAbs) is unsafe: on a
        # small-objective model the solver stops at the first garbage incumbent within
        # it. Instead: a quick estimate solve, then re-solve (warm-started, same model)
        # with a *relative* gap scaled to the objective — max(mip_gap, mip_gap_abs/|obj|)
        # = RAVEN's AbsMIPGap/|objVal|, safe at any scale.
        if time_limit is not None:
            opt.configuration.timeout = max(1, min(int(time_limit) // 10, 30))
        opt.optimize()
        obj = abs(opt.objective.value) if opt.objective.value is not None else 0.0
        _dbg(f"[ftinit] estimate solve: obj={opt.objective.value} status={opt.status}")
        if obj > 0:
            eff_gap = min(max(mip_gap or 0.0, mip_gap_abs / obj), 1.0)
            _dbg(f"[ftinit] eff_gap = max({mip_gap}, {mip_gap_abs}/{obj:.1f}) = {eff_gap:.5f}")
            try:
                opt.problem.Params.MIPGap = eff_gap
            except Exception:  # noqa: BLE001
                pass

    if time_limit is not None:
        opt.configuration.timeout = int(time_limit)
    opt.optimize()
    try:
        _achieved = opt.problem.MIPGap
    except Exception:  # noqa: BLE001
        _achieved = None
    _dbg(f"[ftinit] final solve: obj={opt.objective.value} status={opt.status} "
         f"achieved_gap={_achieved}")
    # Accept a near-optimal incumbent (when a MIP gap / time limit is set), as RAVEN does,
    # but only if the solver actually holds one to read.
    if opt.status not in ("optimal", "feasible", "suboptimal", "time_limit") \
            or not _has_solution(opt):
        raise OptimizationError(
            f"ftINIT MILP produced no usable solution (status: {opt.status}); "
            "increase time_limit or disable strict_gap."
        )

    # Report the *primary* (score) objective; a canonical phase 2 replaces the objective
    # in place, so capture it before that.
    primary_obj = float(opt.objective.value) if opt.objective.value is not None else 0.0

    # RAVEN: a reaction is "on" iff its indicator ≥ 0.5 (positive indicators are
    # continuous and can land fractionally when a reaction can carry only tiny flux).
    def _read_solution():
        on = {rid for rid, (ind, _) in indicators.items() if (ind.primal or 0.0) >= 0.5}
        fluxes = {rid: sum(sign * (var.primal or 0.0) for var, sign in terms)
                  for rid, terms in flux_terms.items()}
        return on, fluxes

    on, fluxes = _read_solution()  # the primary optimum
    # canonical is best-effort: keep its result only if phase 2 actually converged.
    if canonical and indicators and _canonicalize(opt, prob, obj_expr, indicators,
                                                  primary_obj, time_limit):
        on, fluxes = _read_solution()

    kept = free_or_essential | on
    deleted = [r.id for r in model.reactions if r.id not in kept]

    out = model.copy()
    out.remove_reactions(deleted, remove_orphans=True)
    return FtInitResult(out, sorted(kept), sorted(deleted), fluxes,
                        primary_obj, on_reactions=on,
                        achieved_gap=_achieved)


def _has_solution(opt) -> bool:
    """Whether the solver currently holds a readable primal solution.

    A MILP can finish with an accepted status (notably ``time_limit``) yet no incumbent,
    so reading ``.primal`` would raise. On Gurobi this is the solution count; on other
    backends the status is taken as authoritative.
    """
    try:
        return opt.problem.SolCount > 0
    except Exception:  # noqa: BLE001 - non-Gurobi backend
        return True


def _canonicalize(opt, prob, obj_expr, indicators, primary, time_limit) -> bool:
    """Pin ftINIT's degenerate optimum to a single canonical solution (in place on ``opt``).

    The ftINIT MILP is highly degenerate: many reaction subsets reach the same score
    optimum and the solver returns an arbitrary one — reproducible for a fixed solver
    build + seed, but fragile to the Gurobi version or platform (a changed tie-break
    moves ~1-2% of the kept reactions, which flips a handful of gene-essentiality calls).

    This runs a lexicographic phase 2, holding the score objective at its optimum with a
    floor constraint: first minimise the count of kept removable reactions (the sparsest
    optimum), then, among the sparsest, minimise their summed id rank (prefer lower ids).
    The result is a stable, near-unique optimum independent of seed/solver version. ``opt``
    is left holding the phase-2 solution, which the caller reads for the "on" set and
    fluxes; the reported objective stays the phase-1 ``primary`` value.

    Only the negative-score reactions carry a true 0/1 "keep" binary (the positive
    indicators are continuous and pinned near 1 by the score objective), so the two phases
    run over those: their count and id-sum are integers, provable with a cheap absolute gap
    below 1 rather than the near-full proof a tiny relative gap would need. This
    canonicalises the removable-reaction choices — the genuine seed-fragile degeneracy; a
    residual flux-distribution degeneracy (which only feeds the next step's small
    ``ess_force`` clamp) is left to the solver.

    Returns ``True`` if phase 2 produced a usable solution (the caller then reads the
    canonical "on" set and fluxes), ``False`` if it did not converge to an incumbent (the
    caller keeps the phase-1 optimum). Canonicalisation is best-effort: it never fails the
    extraction.
    """
    binaries = {rid: ind for rid, (ind, score) in indicators.items() if score < 0}
    if not binaries:  # only positive / free reactions: nothing removable to canonicalise
        return _has_solution(opt)
    tol = max(abs(primary) * 1e-7, 1e-7)
    opt.add(prob.Constraint(obj_expr, lb=primary - tol, name="_canon_obj_floor"))
    count = add([mul([Real(1.0), ind]) for ind in binaries.values()])

    def _phase(objective) -> bool:
        opt.objective = objective
        try:  # integer objective: an absolute gap < 1 proves the optimum cheaply.
            opt.problem.Params.MIPGap = 0.0
            opt.problem.Params.MIPGapAbs = 0.4
        except Exception:  # noqa: BLE001 - GLPK solves exactly; harmless
            pass
        if time_limit is not None:
            opt.configuration.timeout = int(time_limit)
        opt.optimize()
        return (opt.status in ("optimal", "feasible", "suboptimal", "time_limit")
                and _has_solution(opt))

    # Phase 2a — parsimony: the fewest kept removable reactions.
    if not _phase(prob.Objective(count, direction="min")):
        return False
    kmin = opt.objective.value or 0.0
    # Phase 2b — among the sparsest, prefer lower reaction ids (deterministic tie-break).
    opt.add(prob.Constraint(count, ub=kmin + 0.5, name="_canon_count_cap"))
    ranks = {rid: i for i, rid in enumerate(sorted(binaries))}
    idsum = add([mul([Real(float(1 + ranks[rid])), ind]) for rid, ind in binaries.items()])
    return _phase(prob.Objective(idsum, direction="min"))


def _nudge_scores(rxn_scores: Mapping[str, float]) -> dict[str, float]:
    """Push tiny reaction scores off zero (RAVEN ``ftINIT.m:160-161``).

    Scores in ``(-0.1, 0]`` become ``-0.1`` and in ``(0, 0.1)`` become ``0.1`` — near-zero
    scores cause the MILP numerical trouble. Exactly-zero scores become ``-0.1`` (scored,
    removable) here; genuinely ignore-masked reactions are set back to 0 afterwards by
    :func:`group_rxn_scores`.
    """
    nudged: dict[str, float] = {}
    for rid, s in rxn_scores.items():
        if -0.1 < s <= 0:
            s = -0.1
        elif 0 < s < 0.1:
            s = 0.1
        nudged[rid] = s
    return nudged


def _solve_step(
    min_model, scores, step, *, essential, directions, ess_force, force_on, big_m,
    mip_gap, mip_gap_abs, time_limit, strict_gap=False, canonical=False,
) -> FtInitResult:
    """Solve one ftINIT step, following RAVEN's multi-run gap-escalation schedule.

    RAVEN (``ftINIT.m:227-283``) runs a step up to ``len(step.milp_runs)`` times: run 1 at
    a flat relative gap; each later run at ``max(mip_gap, abs_gap/|obj|)`` — the absolute
    gap made relative, which is safe when the step's objective is near zero — breaking as
    soon as the previous run's achieved gap already meets the next (looser) target. A
    caller ``time_limit`` caps each run's own limit. With no schedule (``step.milp_runs``
    empty, e.g. the ``'full'`` series) this is a single solve at the caller's gap.
    """
    def _run(mg, mga, tl, strict=None):
        return run_ftinit(
            min_model, scores, essential_rxns=essential, essential_directions=directions,
            essential_force=ess_force, allow_excretion=step.allow_met_secr,
            rem_pos_rev=step.pos_rev_off, ignore_mets=step.mets_to_ignore,
            force_on=force_on, force_on_ess=force_on, big_m=big_m,
            mip_gap=mg, mip_gap_abs=mga, time_limit=tl,
            strict_abs_gap=strict, canonical=canonical,
        )

    # Strict mode: one solve per step proven to a fixed absolute gap (below the 0.1 score
    # granularity), bypassing RAVEN's loose relative escalation — whose final
    # near-zero-objective run otherwise accepts an arbitrary within-gap incumbent. Trades
    # runtime for a stable, well-defined optimum (pairs naturally with ``canonical``).
    if strict_gap:
        return _run(None, None, time_limit, strict=_STRICT_ABS_GAP)

    if not step.milp_runs:
        return _run(mip_gap, mip_gap_abs, time_limit)

    res: FtInitResult | None = None
    last_obj: float | None = None
    achieved: float | None = None
    for i, run in enumerate(step.milp_runs):
        target = run["mip_gap"]
        tl = run["time_limit"] if time_limit is None else min(run["time_limit"], time_limit)
        if i > 0:
            # RAVEN: near-zero objective → abs_gap/|obj| blows up → accept any incumbent.
            target = (min(max(run["mip_gap"], run["abs_gap"] / abs(last_obj)), 1.0)
                      if last_obj else 1.0)
            if achieved is not None and achieved <= target:
                break
        res = _run(target, None, tl)
        last_obj, achieved = res.objective, res.achieved_gap
    assert res is not None  # milp_runs is non-empty here
    return res


def ftinit(
    prep,
    rxn_scores: Mapping[str, float],
    *,
    gene_scores: Mapping[str, float] | None = None,
    series: str = "1+1",
    steps=None,
    fill_gaps: bool = True,
    metabolomics: Iterable[str] | None = None,
    force_on: float = _FORCE_ON,
    big_m: float = _BIG_M,
    mip_gap: float | None = None,
    mip_gap_abs: float | None = 10.0,
    time_limit: float | None = None,
    strict_gap: bool = False,
    canonical: bool = False,
) -> cobra.Model:
    """Run the full ftINIT pipeline on prepData and return the context-specific model.

    ``prep`` is a :class:`raven_toolbox.init.PrepData`. ``rxn_scores`` maps **original**
    reaction id → score (e.g. from :func:`score_reactions_from_genes` on the template).
    ``series`` selects which staged schedule to run — see :func:`get_init_steps` for
    what each option changes (default ``'1+1'``); pass ``steps`` directly to use a
    custom schedule instead. Each step regroups scores under its ``ignore_mask``, fixes
    the reactions turned on by earlier steps as essential (in their flux direction), and
    solves :func:`run_ftinit` on the merged model. Reactions never turned on (and not
    essential or left-in) are removed from the reference model; exchange reactions are
    always kept (RAVEN re-adds them).

    If ``fill_gaps`` and ``prep`` carries tasks, reactions are added back so every task
    is feasible (:func:`raven_toolbox.init.fill_tasks`). If ``gene_scores`` is given,
    negative-scoring genes are pruned from the GPRs at the end
    (:func:`raven_toolbox.init.remove_low_score_genes`).

    Essential reactions are forced to carry ``force_on`` (default 0.1) of flux in the
    forced direction. On genome-scale models a stricter regime is needed (the previous
    step's actual carried flux instead of a flat 0.1) — exposed via per-reaction
    ``essential_force`` on :func:`run_ftinit`.

    ``metabolomics`` (a list of detected metabolite names to reward producing) is
    **not yet implemented**: the linear merge eliminates degree-2 detected metabolites,
    so it needs a producer-group-mapping + negative-producer force-flux block — the
    most intricate MILP piece, for the least-used input. Passing a non-empty value
    raises ``NotImplementedError``.

    ``mip_gap``/``time_limit`` are forwarded to each :func:`run_ftinit` solve. On
    genome-scale models they are essential for tractability — see
    ``docs/init_param_calibration.md`` for the calibration table.

    ``strict_gap`` and ``canonical`` (both opt-in, default off → exact RAVEN behaviour)
    reduce the arbitrariness of the degenerate MILP's tie-break, yielding a more
    parsimonious and more *reproducible* extracted model. They do **not** make the model
    biologically more accurate: the alternative optima they choose among are equally
    consistent with the expression data, so this only pins *which* optimum is returned. It
    reduces the fragility of the MILP, not its correctness:

      * ``strict_gap`` replaces the loose relative-gap escalation with a single solve per
        step proven to a fixed *absolute* gap (below the 0.1 reaction-score granularity),
        so the kept-set objective is optimal at any scale — not an arbitrary within-gap
        incumbent. Slower, but removes the largest source of tie-break drift.
      * ``canonical`` adds a lexicographic phase 2 that selects the unique sparsest (then
        lowest-id) optimum, so the degenerate choice is pinned rather than left to the
        solver — applied both to each extraction step and to the task gap-fill. Best used
        together with ``strict_gap`` (a well-defined primary optimum to canonicalise).

    Caveats worth carrying into any workflow that uses these: they reduce run-to-run and
    platform fragility but do *not* guarantee reproducibility across Gurobi versions
    (proving the genome-scale optimum is intractable, so ``strict_gap`` may fall back to an
    incumbent), and a downstream metric such as gene essentiality can even shift or worsen,
    because the sparser ``canonical`` model is more sensitive to the residual (mostly
    transport) degeneracy. For reproducible gene essentiality specifically, pin the solver
    stack (raven-toolbox commit + ``gurobipy`` version) rather than relying on these flags.
    Rule of thumb: for many models the baseline is fine; reach for these when one or a few
    stable, parsimonious model artifacts are wanted.
    """
    if metabolomics:
        raise NotImplementedError(
            "metabolomics production-bonus is not yet implemented."
        )
    steps = steps if steps is not None else get_init_steps(series)
    min_model, group_of = prep.min_model, prep.group_of

    # RAVEN nudges tiny reaction scores off zero once, before per-step grouping.
    rxn_scores = _nudge_scores(rxn_scores)

    turned_on: dict[str, float] = {}   # merged reaction id -> flux (accumulated)
    left_in: set[str] = set()          # merged reactions with score 0 in the last step
    # RAVEN seeds every reaction's "carried flux" at force_on (0.1) and updates it after
    # each step; an essential reaction is forced at min(0.99·|carried flux|, force_on), so
    # it is never forced above what it last carried (ftINIT.m:172,248) — this applies to
    # the permanent (prep) essentials too, not only reactions turned on by a prior step.
    flux_of: dict[str, float] = {r.id: force_on for r in min_model.reactions}
    for step in steps:
        to_zero = prep.masks.ignored(step.ignore_mask)
        scores = group_rxn_scores(min_model, rxn_scores, prep.orig_rxn_ids,
                                  prep.group_ids, to_zero)
        essential = set(prep.essential_rxns)  # pre-oriented forward (default direction)
        directions: dict[str, int] = {}
        if step.how_to_use_prev == "essential":
            for rid, flux in turned_on.items():
                essential.add(rid)
                directions[rid] = 1 if flux >= 0 else -1
        ess_force = {rid: min(abs(flux_of.get(rid, force_on)) * 0.99, force_on)
                     for rid in essential}
        res = _solve_step(
            min_model, scores, step, essential=essential, directions=directions,
            ess_force=ess_force, force_on=force_on, big_m=big_m,
            mip_gap=mip_gap, mip_gap_abs=mip_gap_abs, time_limit=time_limit,
            strict_gap=strict_gap, canonical=canonical,
        )
        for rid in res.on_reactions:
            turned_on[rid] = res.fluxes[rid]
        flux_of.update(res.fluxes)  # carry this step's fluxes into the next
        left_in = {rid for rid, s in scores.items() if s == 0.0}

    # Merged reactions to keep: turned on + permanently essential + left-in (score 0).
    kept_min = set(turned_on) | set(prep.essential_rxns) | left_in
    deleted_min = [r.id for r in min_model.reactions if r.id not in kept_min]

    # Map deleted merged reactions back to all originals in their groups.
    removed_groups = {group_of[rid] for rid in deleted_min if group_of[rid] != 0}
    to_remove = {o for o in prep.orig_rxn_ids if group_of[o] and group_of[o] in removed_groups}
    to_remove |= {rid for rid in deleted_min if group_of[rid] == 0}  # unmerged
    # Keep the surviving originals plus all exchange reactions (always re-added).
    final_kept = (set(prep.orig_rxn_ids) - to_remove) | prep.masks.exchange

    out = prep.ref_model.copy()
    out.remove_reactions([r.id for r in out.reactions if r.id not in final_kept],
                         remove_orphans=True)

    if fill_gaps and prep.tasks:  # add reactions back so every task is feasible
        # The gap-fill MILP is its own problem (RAVEN ftINITFillGaps); it uses RAVEN's
        # fixed per-task 300 s limit and seed, not the main extraction's time_limit.
        out = fill_tasks(out, prep.ref_model, prep.tasks, rxn_scores=rxn_scores,
                         canonical=canonical).model
    if gene_scores is not None:   # prune negative-scoring genes from the GPRs
        out, _ = remove_low_score_genes(out, gene_scores)
    return out
