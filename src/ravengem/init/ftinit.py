"""The ftINIT MILP (port of RAVEN ``ftINITInternalAlg``) — Phase 4d.3.

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

Objective: **maximise** ``Σ score·indicator``. Unlike classic INIT (our
:func:`ravengem.init.run_init`), ftINIT does **not** reward production of every
metabolite — ``prodWeight`` applies only to metabolomics-detected metabolites, which
are deferred to 4d.6; here connectivity comes solely from the flux gates plus any
essential reactions. ``allow_excretion`` relaxes ``S·v = 0`` to ``≥ 0``; ``rem_pos_rev``
drops positive reversible reactions from the problem (used in staging, 4d.3b).

Needs a MILP solver (cobra's configured optlang solver). Magic numbers
(``force_on``/``force_on_ess`` = 0.1, big-M = each reaction's own bound) are exposed;
they are scale-dependent and calibrated in 4d.7 (see docs/ftinit_review_and_plan.md).

⚠️ **Loops.** Like RAVEN's MILP, this has *no* loopless constraint: an internal
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

_FORCE_ON = 0.1  # min flux for a reaction to count as "on" (RAVEN forceOnLim)


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
    allow_excretion: bool = False,
    rem_pos_rev: bool = False,
    ignore_mets: Iterable[str] = (),
    force_on: float = _FORCE_ON,
    force_on_ess: float = _FORCE_ON,
) -> FtInitResult:
    """Run the single-step ftINIT MILP and return the extracted model.

    ``rxn_scores`` maps reaction id → score (default 0 → reaction left free in the
    model, not scored or removable). ``essential_rxns`` are forced to carry flux
    (≥ ``force_on_ess``); ``essential_directions`` maps an essential reaction id to
    ``+1`` (forward) or ``-1`` (reverse) for the forced direction (default forward).
    ``ignore_mets`` are metabolite **names** whose mass balance is dropped (RAVEN's
    per-step "simple metabolite" removal, e.g. H2O/H+). See the module docstring for
    the formulation. This is the single-step variant; the staged schedule
    (:func:`ravengem.init.ftinit`) calls it per step.
    """
    scores = dict(rxn_scores or {})
    essential = set(essential_rxns or [])
    directions = dict(essential_directions or {})
    ignore_met_names = set(ignore_mets)
    prob = model.problem
    opt = prob.Model()

    variables: list = []
    constraints: list = []
    flux_terms: dict[str, list[tuple[object, float]]] = {}  # rxn id -> [(var, sign)]
    indicators: dict[str, tuple[object, float]] = {}  # rxn id -> (indicator var, score)
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
            # stricter native bound if the model already forces more flux.
            if directions.get(rid, 1) >= 0:
                v = prob.Variable(f"v_{rid}", lb=max(force_on_ess, lb), ub=ub)
            else:  # reverse: flux ≤ -force_on_ess
                v = prob.Variable(f"v_{rid}", lb=lb, ub=min(-force_on_ess, ub))
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
                add_constraint(vp - ub * b, ub=0.0, name=f"dirp_{rid}")        # vp ≤ ub·b
                add_constraint(vn - lb * b, ub=-lb, name=f"dirn_{rid}")        # vn ≤ -lb·(1-b)
        else:  # score < 0
            x = prob.Variable(f"x_{rid}", type="binary")
            variables.append(x)
            indicators[rid] = (x, score)
            cap = (ub - lb) if reversible else (ub if ub > 0 else -lb)
            add_constraint(total - cap * x, ub=0.0, name=f"off_{rid}")  # flux>0 ⇒ x=1

    def net(rid):  # net flux expression vp - vn (or v), the single source of truth
        return sum(sign * var for var, sign in flux_terms[rid])

    # Steady state S·v {== 0 | >= 0}; ignored metabolites are left unbalanced.
    for met in model.metabolites:
        if met.name in ignore_met_names:
            continue
        expr = sum(coeff * net(r.id) for r in met.reactions if (coeff := r.metabolites[met]))
        if expr != 0:
            add_constraint(expr, lb=0.0, ub=None if allow_excretion else 0.0)

    opt.add(variables + constraints)
    opt.objective = prob.Objective(
        sum(score * ind for ind, score in indicators.values()), direction="max"
    )
    opt.optimize()
    if opt.status != "optimal":
        raise RuntimeError(f"ftINIT MILP did not solve to optimality (status: {opt.status}).")

    # RAVEN: a reaction is "on" iff its indicator ≥ 0.5 (positive indicators are
    # continuous and can land fractionally when a reaction can carry only tiny flux).
    on = {rid for rid, (ind, _) in indicators.items() if (ind.primal or 0.0) >= 0.5}
    kept = free_or_essential | on
    deleted = [r.id for r in model.reactions if r.id not in kept]
    fluxes = {
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
    series: str = "1+1",
    steps=None,
    force_on: float = _FORCE_ON,
) -> cobra.Model:
    """Run the staged ftINIT pipeline on prepData and return the extracted model.

    ``prep`` is a :class:`ravengem.init.PrepData`. ``rxn_scores`` maps **original**
    reaction id → score (e.g. from :func:`score_reactions_from_genes` on the template).
    Each step (:func:`ravengem.init.get_init_steps`) regroups scores under its
    ``ignore_mask``, fixes the reactions turned on by earlier steps as essential (in
    their flux direction), and solves :func:`run_ftinit` on the merged model. Reactions
    never turned on (and not essential or left-in) are removed from the reference model;
    exchange reactions are always kept (RAVEN re-adds them).
    """
    from ravengem.init.merge import group_rxn_scores
    from ravengem.init.steps import get_init_steps

    steps = steps if steps is not None else get_init_steps(series)
    min_model, group_of = prep.min_model, prep.group_of

    turned_on: dict[str, float] = {}   # merged reaction id -> flux (accumulated)
    left_in: set[str] = set()          # merged reactions with score 0 in the last step
    for step in steps:
        to_zero = prep.masks.ignored(step.ignore_mask)
        scores = group_rxn_scores(min_model, rxn_scores, prep.orig_rxn_ids,
                                  prep.group_ids, to_zero)
        essential = set(prep.essential_rxns)
        directions = dict(prep.essential_directions)
        if step.how_to_use_prev == "essential":
            for rid, flux in turned_on.items():
                essential.add(rid)
                directions[rid] = 1 if flux >= 0 else -1
        res = run_ftinit(
            min_model, scores, essential_rxns=essential, essential_directions=directions,
            allow_excretion=step.allow_met_secr, rem_pos_rev=step.pos_rev_off,
            ignore_mets=step.mets_to_ignore, force_on=force_on, force_on_ess=force_on,
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
    return out
