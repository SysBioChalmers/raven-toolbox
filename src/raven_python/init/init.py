"""The INIT MILP — tINIT core.

INIT (Agren et al., PLoS Comput Biol 2012) extracts a context-specific model: keep a
flux-consistent subnetwork that maximises the summed score of *included* reactions
(positive score = evidence to keep, negative = evidence to remove), optionally
rewarding net production of metabolites.

Formulation:

* Reversible reactions are split into forward / reverse directed reactions (flux ≥ 0).
* Each non-essential directed reaction gets a binary ``x`` (included ⇔ ``x=1``) with
  ``eps·x ≤ v ≤ ub·x`` — included reactions must carry flux ≥ ``eps`` (connectivity),
  excluded ones carry none.
* Essential reactions (``essential_rxns``) are forced to carry flux (``v ≥ eps``) and
  skip the binary.
* ``no_rev_loops`` adds ``x_fwd + x_rev ≤ 1`` so a reversible reaction can't look
  "connected" via an internal forward/back loop.
* Steady state ``S·v = 0`` per metabolite; ``allow_excretion`` relaxes it to ``≥ 0``
  (net production allowed). With ``prod_weight > 0`` a per-metabolite sink
  ``s_m ∈ [0,1]`` is added and rewarded, giving a reason to include connectivity
  reactions.
* Objective: **maximise** ``Σ score·x + prod_weight·Σ s_m``.

Needs a MILP solver (cobra's configured optlang solver). On genome-scale problems,
Gurobi is the only backend that is fully usable today (see
``docs/init_solver_benchmark.md``).

**Parameter caveat — magic numbers are scale-dependent.** ``eps`` (the flux an
included reaction must carry, default 1.0) and ``prod_weight`` (default 0.5) only make
sense when reaction bounds are ~±1000 and scores are O(1); the right values depend on
the model's flux magnitudes and the score distribution. The upper gate uses each
reaction's own ``ub`` as the big-M by default (adapts to the model); pass ``big_m`` to
override with a fixed cap for a tighter LP relaxation. Calibration tables live in
``docs/init_param_calibration.md``.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import cobra
from cobra.exceptions import OptimizationError
from optlang.symbolics import Real, add, mul

_EPS = 1.0  # flux an included reaction must carry (RAVEN's fake-met unit)


@dataclass
class _Directed:
    """One directed reaction in the split (irreversible) problem."""

    key: str
    origin: str  # original reaction id
    coeffs: dict[str, float]  # met id -> stoichiometry (already sign-adjusted)
    ub: float
    score: float
    essential: bool


@dataclass
class InitResult:
    """Result of :func:`run_init`."""

    model: cobra.Model
    deleted_reactions: list[str]
    met_production: dict[str, bool]  # present-met name -> producible?
    objective: float


def _split_reactions(
    model: cobra.Model, scores: Mapping[str, float], essential: set[str]
) -> list[_Directed]:
    directed: list[_Directed] = []
    for rxn in model.reactions:
        score = float(scores.get(rxn.id, 0.0))
        coeffs = {m.id: c for m, c in rxn.metabolites.items()}
        rev_coeffs = {m: -c for m, c in coeffs.items()}
        if rxn.id in essential:
            # Force flux in a *single* direction (forward if it can run forward, else
            # reverse) — like an irreversible essential reaction. Emitting both halves
            # as essential would force fwd ≥ eps AND rev ≥ eps, i.e. a phantom
            # eps-magnitude self-loop that can starve out the real pathway.
            if rxn.upper_bound > 0:
                directed.append(_Directed(rxn.id, rxn.id, coeffs, rxn.upper_bound, score, True))
            else:
                directed.append(_Directed(f"{rxn.id}__rev", rxn.id, rev_coeffs,
                                          -rxn.lower_bound, score, True))
            continue
        if rxn.upper_bound > 0:
            directed.append(_Directed(rxn.id, rxn.id, coeffs, rxn.upper_bound, score, False))
        if rxn.lower_bound < 0:  # reverse direction as its own non-negative flux
            directed.append(
                _Directed(f"{rxn.id}__rev", rxn.id, rev_coeffs, -rxn.lower_bound, score, False)
            )
    return directed


def run_init(
    model: cobra.Model,
    rxn_scores: Mapping[str, float] | None = None,
    *,
    present_mets: Iterable[str] | None = None,
    essential_rxns: Iterable[str] | None = None,
    prod_weight: float = 0.5,
    allow_excretion: bool = False,
    no_rev_loops: bool = False,
    eps: float = _EPS,
    big_m: float | None = None,
    mip_gap: float | None = None,
    time_limit: float | None = None,
) -> InitResult:
    """Run the INIT MILP and return the extracted model.

    ``rxn_scores`` maps reaction id → score (default 0). ``essential_rxns`` must be
    kept (forced to carry flux). ``present_mets`` are metabolite *names* that the
    network should be able to produce; each is tested and reported in
    ``met_production``. See the module docstring for the formulation.

    Note on score 0 (classic INIT vs. ftINIT divergence): in classic INIT a
    reaction with score exactly 0 receives an include-indicator with **zero
    reward**, so the optimiser is free to drop it. This matches RAVEN's
    `runINIT` semantics. ftINIT inverts that — score-0 reactions stay in the
    model unless they actively hurt feasibility — so a score of exactly 0
    means *different things* in the two variants. If you want score-0
    reactions kept here, pass a small positive value (e.g. ``min_score`` from
    `gene_scores_from_expression`) instead of 0.
    """
    scores = dict(rxn_scores or {})
    essential = set(essential_rxns or [])
    present = list(present_mets or [])

    directed = _split_reactions(model, scores, essential)
    prob = model.problem
    opt = prob.Model()

    # Flux variables for every directed reaction.
    flux = {d.key: prob.Variable(f"v_{d.key}", lb=0.0, ub=d.ub) for d in directed}

    # Binary include-indicators for non-essential reactions; eps*x <= v <= ub*x.
    keep: dict[str, object] = {}
    gates = []
    for d in directed:
        if d.essential:
            flux[d.key].lb = max(eps, 0.0)  # forced to carry flux
            continue
        x = prob.Variable(f"x_{d.key}", type="binary")
        keep[d.key] = x
        cap = d.ub if big_m is None else big_m  # big-M: per-reaction bound (default) or fixed
        gates.append(prob.Constraint(flux[d.key] - cap * x, ub=0.0, name=f"ub_{d.key}"))
        gates.append(prob.Constraint(flux[d.key] - eps * x, lb=0.0, name=f"lb_{d.key}"))

    # no_rev_loops: at most one direction of a reversible reaction is included.
    by_origin: dict[str, list[str]] = {}
    for d in directed:
        by_origin.setdefault(d.origin, []).append(d.key)
    if no_rev_loops:
        for keys in by_origin.values():
            xs = [keep[k] for k in keys if k in keep]
            if len(xs) > 1:
                gates.append(prob.Constraint(sum(xs), ub=1.0, name=f"onedir_{keys[0]}"))

    # Steady-state constraints S·v (- sink) {==0 | >=0}, plus prod_weight sinks.
    # Accumulate each metabolite's terms by iterating reactions once (avoids the
    # O(mets·rxns) per-metabolite filter) and sum with optlang.symbolics.add — Python
    # sum() re-canonicalises a growing sympy expression each step (O(n²)), which is
    # minutes per hub metabolite at genome scale.
    met_terms: dict[str, list] = {met.id: [] for met in model.metabolites}
    for d in directed:
        v = flux[d.key]
        for mid, coeff in d.coeffs.items():
            met_terms[mid].append(mul([Real(coeff), v]))

    sinks: dict[str, object] = {}
    met_constraints: dict[str, object] = {}
    ub = None if allow_excretion else 0.0
    for met in model.metabolites:
        terms = met_terms[met.id]
        if prod_weight != 0:
            s = prob.Variable(f"s_{met.id}", lb=0.0, ub=1.0)
            sinks[met.id] = s
            terms = [*terms, mul([Real(-1.0), s])]  # net production drained into rewarded sink
        if terms:
            met_constraints[met.id] = prob.Constraint(add(terms), lb=0.0, ub=ub)

    opt.add(list(flux.values()) + list(keep.values()) + list(sinks.values())
            + gates + list(met_constraints.values()))

    objective = prob.Objective(
        add([mul([Real(d.score), keep[d.key]]) for d in directed if d.key in keep]
            + [mul([Real(prod_weight), s]) for s in sinks.values()]),
        direction="max",
    )
    opt.objective = objective

    met_production = _check_present_mets(prob, present, model, directed, allow_excretion)

    if time_limit is not None:
        opt.configuration.timeout = int(time_limit)
    if mip_gap is not None:
        try:  # Gurobi-specific; harmless if the backend differs
            opt.problem.Params.MIPGap = mip_gap
        except Exception:  # noqa: BLE001
            pass
    opt.optimize()
    # With a MIP gap / time limit set, accept a near-optimal incumbent (as RAVEN does).
    if opt.status not in ("optimal", "feasible", "suboptimal", "time_limit"):
        raise OptimizationError(f"INIT MILP did not solve (status: {opt.status}).")

    # A reaction is kept if any of its directed parts is essential or has x≈1.
    kept_origins = {d.origin for d in directed if d.essential}
    kept_origins |= {d.origin for d in directed if d.key in keep and (keep[d.key].primal or 0) > 0.5}
    deleted = [r.id for r in model.reactions if r.id not in kept_origins]

    out = model.copy()
    out.remove_reactions(deleted, remove_orphans=True)
    return InitResult(out, sorted(deleted), met_production, float(opt.objective.value))


def _check_present_mets(prob, present, model, directed, allow_excretion) -> dict[str, bool]:
    """Whether each present metabolite (by name) can be net-produced at all.

    A small LP per metabolite (no score/binary, so it's the LP relaxation, as RAVEN
    does): all reactions available, steady state, and a demand draining ≥1 unit of
    any compartment form of the metabolite — feasible ⇔ producible.
    """
    if not present:
        return {}
    name_to_ids: dict[str, list[str]] = {}
    for met in model.metabolites:
        name_to_ids.setdefault((met.name or met.id).upper(), []).append(met.id)

    result: dict[str, bool] = {}
    for name in present:
        ids = name_to_ids.get(name.upper())
        if not ids:
            result[name] = False
            continue
        lp = prob.Model()
        flux = {d.key: prob.Variable(f"v_{d.key}", lb=0.0, ub=d.ub) for d in directed}
        drains = {mid: prob.Variable(f"drain_{mid}", lb=0.0, ub=1e6) for mid in ids}
        terms: dict[str, list] = {met.id: [] for met in model.metabolites}
        for d in directed:
            v = flux[d.key]
            for mid, c in d.coeffs.items():
                terms[mid].append(mul([Real(c), v]))
        for mid in drains:
            terms[mid].append(mul([Real(-1.0), drains[mid]]))
        cons = [prob.Constraint(add(t), lb=0.0, ub=None if allow_excretion else 0.0)
                for t in terms.values() if t]
        require = prob.Constraint(add(list(drains.values())), lb=1.0, name="_require_production")
        lp.add(list(flux.values()) + list(drains.values()) + cons + [require])
        lp.objective = prob.Objective(prob.Variable("_zero", lb=0, ub=0), direction="max")
        lp.optimize()
        result[name] = lp.status == "optimal"
    return result
