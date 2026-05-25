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
from dataclasses import dataclass

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


def run_ftinit(
    model: cobra.Model,
    rxn_scores: Mapping[str, float] | None = None,
    *,
    essential_rxns: Iterable[str] | None = None,
    allow_excretion: bool = False,
    rem_pos_rev: bool = False,
    force_on: float = _FORCE_ON,
    force_on_ess: float = _FORCE_ON,
) -> FtInitResult:
    """Run the single-step ftINIT MILP and return the extracted model.

    ``rxn_scores`` maps reaction id → score (default 0 → reaction left free in the
    model, not scored or removable). ``essential_rxns`` are forced to carry flux
    (≥ ``force_on_ess``) and should already be oriented irreversibly. See the module
    docstring for the formulation. This is the ``'full'`` (single-step) variant;
    linear merging and the staged ``'1+1'`` schedule are layered on later (4d.2/4d.3b).
    """
    scores = dict(rxn_scores or {})
    essential = set(essential_rxns or [])
    prob = model.problem
    opt = prob.Model()

    variables: list = []
    constraints: list = []
    net_flux: dict[str, object] = {}                  # rxn id -> optlang expr (for S·v)
    flux_terms: dict[str, list[tuple[object, float]]] = {}  # rxn id -> [(var, sign)] for values
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
            # Forced on, oriented forward (prepINITModel makes essentials irreversible).
            v = prob.Variable(f"v_{rid}", lb=max(force_on_ess, max(lb, 0.0)), ub=ub)
            variables.append(v)
            net_flux[rid], flux_terms[rid] = v, [(v, 1.0)]
            free_or_essential.add(rid)
            continue

        if score == 0.0:  # free: carries flux for connectivity, not scored/removable
            v = prob.Variable(f"v_{rid}", lb=lb, ub=ub)
            variables.append(v)
            net_flux[rid], flux_terms[rid] = v, [(v, 1.0)]
            free_or_essential.add(rid)
            continue

        reversible = lb < 0 < ub
        if reversible:
            vp = prob.Variable(f"vp_{rid}", lb=0.0, ub=ub)
            vn = prob.Variable(f"vn_{rid}", lb=0.0, ub=-lb)
            variables += [vp, vn]
            net_flux[rid], flux_terms[rid] = vp - vn, [(vp, 1.0), (vn, -1.0)]
            total = vp + vn  # |flux|, used by the on/off gates
        else:
            v = prob.Variable(f"v_{rid}", lb=min(lb, 0.0), ub=max(ub, 0.0))
            variables.append(v)
            net_flux[rid], flux_terms[rid] = v, [(v, 1.0)]
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

    # Steady state S·v {== 0 | >= 0}.
    for met in model.metabolites:
        expr = sum(coeff * net_flux[r.id] for r in met.reactions if (coeff := r.metabolites[met]))
        if expr != 0:
            add_constraint(expr, lb=0.0, ub=None if allow_excretion else 0.0)

    opt.add(variables + constraints)
    opt.objective = prob.Objective(
        sum(score * ind for ind, score in indicators.values()), direction="max"
    )
    opt.optimize()
    if opt.status != "optimal":
        raise RuntimeError(f"ftINIT MILP did not solve to optimality (status: {opt.status}).")

    on = {rid for rid, (ind, _) in indicators.items() if (ind.primal or 0.0) > 0.5}
    kept = free_or_essential | on
    deleted = [r.id for r in model.reactions if r.id not in kept]
    fluxes = {
        rid: sum(sign * (var.primal or 0.0) for var, sign in terms)
        for rid, terms in flux_terms.items()
    }

    out = model.copy()
    out.remove_reactions(deleted, remove_orphans=True)
    return FtInitResult(out, sorted(kept), sorted(deleted), fluxes, float(opt.objective.value))
