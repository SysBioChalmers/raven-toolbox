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


@dataclass
class FtInitResult:
    """Result of :func:`run_ftinit`."""

    model: cobra.Model
    kept_reactions: list[str]
    deleted_reactions: list[str]
    fluxes: dict[str, float]
    objective: float
    on_reactions: set[str] = field(default_factory=set)  # scored reactions turned on (indicator)


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
    time_limit: float | None = None,
) -> FtInitResult:
    """Run the single-step ftINIT MILP and return the extracted model.

    ``rxn_scores`` maps reaction id → score (default 0 → reaction left free in the
    model, not scored or removable). ``essential_rxns`` are forced to carry flux
    (≥ ``force_on_ess``); ``essential_directions`` maps an essential reaction id to
    ``+1`` (forward) or ``-1`` (reverse) for the forced direction (default forward).
    ``ignore_mets`` are metabolite **names** whose mass balance is dropped (RAVEN's
    per-step "simple metabolite" removal, e.g. H2O/H+). See the module docstring for
    the formulation. This is the single-step variant; the staged schedule
    (:func:`raven_toolbox.init.ftinit`) calls it per step.
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
    opt.objective = prob.Objective(
        add([mul([Real(score), ind]) for ind, score in indicators.values()]), direction="max"
    )
    if time_limit is not None:
        opt.configuration.timeout = int(time_limit)
    if mip_gap is not None:
        try:  # Gurobi-specific; harmless if the backend differs
            opt.problem.Params.MIPGap = mip_gap
        except Exception:  # noqa: BLE001
            pass
    opt.optimize()
    # Accept a near-optimal incumbent (when a MIP gap / time limit is set), as RAVEN does.
    if opt.status not in ("optimal", "feasible", "suboptimal", "time_limit"):
        raise OptimizationError(f"ftINIT MILP did not solve (status: {opt.status}).")

    # RAVEN: a reaction is "on" iff its indicator ≥ 0.5 (positive indicators are
    # continuous and can land fractionally when a reaction can carry only tiny flux).
    on = {rid for rid, (ind, _) in indicators.items() if (ind.primal or 0.0) >= 0.5}
    kept = free_or_essential | on
    deleted = [r.id for r in model.reactions if r.id not in kept]
    fluxes: dict[str, float] = {
        rid: sum(sign * (var.primal or 0.0) for var, sign in terms)
        for rid, terms in flux_terms.items()
    }

    out = model.copy()
    out.remove_reactions(deleted, remove_orphans=True)
    return FtInitResult(out, sorted(kept), sorted(deleted), fluxes,
                        float(opt.objective.value), on_reactions=on)


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
    time_limit: float | None = None,
) -> cobra.Model:
    """Run the full ftINIT pipeline on prepData and return the context-specific model.

    ``prep`` is a :class:`raven_toolbox.init.PrepData`. ``rxn_scores`` maps **original**
    reaction id → score (e.g. from :func:`score_reactions_from_genes` on the template).
    Each step (:func:`raven_toolbox.init.get_init_steps`) regroups scores under its
    ``ignore_mask``, fixes the reactions turned on by earlier steps as essential (in
    their flux direction), and solves :func:`run_ftinit` on the merged model. Reactions
    never turned on (and not essential or left-in) are removed from the reference model;
    exchange reactions are always kept (RAVEN re-adds them).

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
    """
    if metabolomics:
        raise NotImplementedError(
            "metabolomics production-bonus is not yet implemented."
        )
    steps = steps if steps is not None else get_init_steps(series)
    min_model, group_of = prep.min_model, prep.group_of

    turned_on: dict[str, float] = {}   # merged reaction id -> flux (accumulated)
    left_in: set[str] = set()          # merged reactions with score 0 in the last step
    for step in steps:
        to_zero = prep.masks.ignored(step.ignore_mask)
        scores = group_rxn_scores(min_model, rxn_scores, prep.orig_rxn_ids,
                                  prep.group_ids, to_zero)
        essential = set(prep.essential_rxns)  # pre-oriented forward (default direction)
        directions: dict[str, int] = {}
        ess_force: dict[str, float] = {}
        if step.how_to_use_prev == "essential":
            for rid, flux in turned_on.items():
                essential.add(rid)
                directions[rid] = 1 if flux >= 0 else -1
                # never force more flux than the reaction carried before (RAVEN)
                ess_force[rid] = min(abs(flux) * 0.99, force_on)
        res = run_ftinit(
            min_model, scores, essential_rxns=essential, essential_directions=directions,
            essential_force=ess_force, allow_excretion=step.allow_met_secr,
            rem_pos_rev=step.pos_rev_off, ignore_mets=step.mets_to_ignore,
            force_on=force_on, force_on_ess=force_on, big_m=big_m,
            mip_gap=mip_gap, time_limit=time_limit,
        )
        for rid in res.on_reactions:
            turned_on[rid] = res.fluxes[rid]
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
        out = fill_tasks(out, prep.ref_model, prep.tasks, rxn_scores=rxn_scores,
                         mip_gap=mip_gap, time_limit=time_limit).model
    if gene_scores is not None:   # prune negative-scoring genes from the GPRs
        out, _ = remove_low_score_genes(out, gene_scores)
    return out
