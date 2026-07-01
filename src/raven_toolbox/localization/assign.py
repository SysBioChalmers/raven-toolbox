"""Functionality-constrained, multi-localization compartment assignment.

``assign_compartments`` places a set of reactions into subcellular compartments by a single
MILP that maximises agreement with soft per-compartment localisation scores **subject to the
compartmentalised network still producing biomass** (or reaching a growth floor). The biomass
constraint is what makes the result functional and connected, and it is what makes a reaction
land in a compartment *against its own top score* when its pathway needs it there — pathway
coherence is an emergent property of requiring flux, not a bespoke heuristic.

The MILP fuses two pieces:

* the assignment structure (binary ``x[r,c]`` reaction-in-compartment, ``y[g,c]`` gene-in-
  compartment, ``t[m,c]`` transport, with multi-localisation), and
* a **compartment-expanded flux** model (continuous fluxes over per-compartment metabolite
  nodes, Big-M-gated by the placement binaries, with ``v_biomass >= min_growth``).

Localisation evidence is an agnostic ``gene x compartment`` score table
(:class:`raven_toolbox.localization.scores.LocalizationScores`) — any predictor or database.

Built on raven-toolbox + cobra.
"""
from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field

import cobra
from cobra.util import linear_reaction_coefficients
from optlang.symbolics import Real, add, mul

from raven_toolbox.localization.scores import LocalizationScores

__all__ = ["AssignmentProposal", "assign_compartments", "apply_assignment"]


# --------------------------------------------------------------------------- helpers
def _base_met(m: cobra.Metabolite) -> str:
    if m.compartment and m.id.endswith(f"_{m.compartment}"):
        return m.id[: -(len(m.compartment) + 1)]
    return m.id


def _rxn_compartment(rxn: cobra.Reaction) -> str | None:
    comps = {m.compartment for m in rxn.metabolites if m.compartment}
    return next(iter(comps)) if len(comps) == 1 else None


def _tighten_integrality(opt, tol: float) -> None:
    """Best-effort tightening of the MILP integer-feasibility tolerance across solver backends.

    The Big-M flux gating is only sound if the integer tolerance is tight: otherwise ghost flux
    up to ``ub * tol`` leaks through a reaction whose placement binary is rounded to zero, which
    can spuriously satisfy the growth floor. Backends without the knob are left at their default.
    """
    if not tol or tol <= 0:
        return
    # optlang's solver-agnostic knob first (covers GLPK, the CI solver), then backend-specific.
    try:
        opt.configuration.tolerances.integrality = tol
        return
    except Exception:  # noqa: BLE001 — backend may not expose it via optlang
        pass
    prob = getattr(opt, "problem", None)
    if prob is None:
        return
    for setter in (
        lambda: setattr(prob.Params, "IntFeasTol", tol),                # Gurobi
        lambda: prob.parameters.mip.tolerances.integrality.set(tol),    # CPLEX
        lambda: prob.setOptionValue("mip_feasibility_tolerance", tol),  # HiGHS
    ):
        try:
            setter()
            return
        except Exception:  # noqa: BLE001 — wrong backend for this knob; try the next
            continue


@dataclass
class AssignmentProposal:
    """Outcome of :func:`assign_compartments`.

    Attributes
    ----------
    placements:
        ``reaction_id -> [compartments]`` the reaction is assigned to (one or more).
    gene_compartments:
        ``gene_id -> [compartments]``.
    added_transports:
        ``(base_metabolite_id, compartment)`` pairs needing a transport to/from
        ``default_compartment``.
    added_reactions:
        ids of gap-fill reactions pulled from the ``universal`` model to keep the
        compartmentalised network functional.
    unplaced_reactions:
        relocate-set reactions with no scored gene support (placed by function only).
    objective:
        MILP objective (localisation score − transport − multi-penalty − gap-fill cost).
    min_growth:
        the growth floor enforced.
    status:
        solver status (``"optimal"``, ``"infeasible"``, ...).
    """

    placements: dict[str, list[str]] = field(default_factory=dict)
    gene_compartments: dict[str, list[str]] = field(default_factory=dict)
    added_transports: list[tuple[str, str]] = field(default_factory=list)
    added_reactions: list[str] = field(default_factory=list)
    unplaced_reactions: list[str] = field(default_factory=list)
    objective: float = 0.0
    min_growth: float = 0.0
    status: str = "not_solved"


def assign_compartments(
    model: cobra.Model,
    scores: LocalizationScores,
    reactions_to_relocate: Iterable[str],
    *,
    default_compartment: str = "c",
    transport_cost: float | Mapping[str, float] = 0.5,
    multi_compartment_penalty: float = 0.5,
    transportable: Iterable[str] | None = None,
    base_metabolite: Callable[[cobra.Metabolite], str] | None = None,
    universal: cobra.Model | None = None,
    gapfill_cost: float | Mapping[str, float] = 1.0,
    biomass_reaction: str | None = None,
    min_growth: float | None = None,
    big_m: float = 1000.0,
    time_limit: float | None = None,
    mip_gap: float | None = None,
    integrality_tol: float = 1e-9,
    multi_localization: bool = False,
    eps_flux: float | None = None,
) -> AssignmentProposal:
    """Assign ``reactions_to_relocate`` to compartments while keeping the model functional.

    Parameters
    ----------
    model:
        Draft model (one or few compartments) with GPRs and a biomass objective.
    scores:
        Agnostic ``gene x compartment`` soft scores (higher = stronger evidence).
    reactions_to_relocate:
        Reaction ids to (re)place. Boundary reactions, existing multi-compartment reactions,
        and the biomass reaction are always pinned (filtered out even if listed).
    default_compartment:
        Compartment that transports route through (usually cytosol).
    transport_cost:
        Scalar or ``{base_met_id: cost}`` cost per added inter-compartment transport.
    multi_compartment_penalty:
        Cost per *extra* compartment a gene is placed in (its primary is free).
    transportable:
        Base metabolite ids that may receive inter-compartment transports. ``None`` = all.
        Restricting it lets compartment-confined metabolites force functionality-driven
        placement (a reaction must sit where its substrates actually are).
    base_metabolite:
        Callable mapping a metabolite to its **compartment-agnostic key** (the same chemical
        across compartments must map to one key). Default strips an ``_<compartment>`` id
        suffix (e.g. ``atp_c`` → ``atp``). Models that key the same species to different ids
        per compartment — e.g. yeast-GEM's ``s_####`` ids — need ``base_metabolite=lambda m:
        m.name`` (or an annotation-based key) so the compartments are unified.
    universal:
        Optional candidate-reaction database. Reactions whose id is not already in ``model``
        become gap-fill candidates: the MILP may switch them on (at ``gapfill_cost`` each) to
        keep the compartmentalised network functional. Because candidates earn no localisation
        score and only cost, they are added only when biomass feasibility requires them — the
        same parsimonious gap-fill pattern as cost-weighted gap-filling, fused into assignment.
    gapfill_cost:
        Scalar or ``{reaction_id: cost}`` cost per gap-fill reaction added from ``universal``.
    biomass_reaction:
        Reaction id whose flux must reach ``min_growth``. Default: the model objective.
    min_growth:
        Required biomass flux. Default: 10 % of the un-compartmentalised FBA optimum. If the
        draft cannot grow on its own (a true gap), pass ``min_growth`` explicitly and supply
        ``universal`` so the gap-fill can restore growth.
    big_m:
        Big-M for transport-flux gating.
    integrality_tol:
        Integer-feasibility tolerance applied to the solver (best-effort across backends). Must
        be tight enough that ``ub * integrality_tol`` is well below ``min_growth``; otherwise the
        Big-M flux gating leaks ghost flux through a reaction placed elsewhere and can certify a
        non-functional placement. Default ``1e-9``.
    multi_localization:
        Opt-in **reaction-level multi-localisation** (default ``False`` = mono-localisation). When
        ``True`` a reaction may be placed in several compartments at once (``Σ_c x[r,c] ≥ 1``), but
        every *extra* placement must **carry flux** (``flux-activity coupling``): a placement earns
        its localisation score only if it is *active* — i.e. its flux magnitude reaches ``eps_flux``
        in the same biomass-producing solution. Each reaction may keep at most one **inactive
        "home"** placement (scored by prior alone, exactly like a mono-localised reaction that is
        not needed for biomass), and that home is allowed **only if the reaction carries no flux at
        all**; the moment a reaction is used, *all* its placements must be flux-active. This makes
        the feature **sound**: a reaction can never collect a high-score compartment's reward
        through a *dead* duplicate (zero-flux placement). It is solver-independent — the soundness
        is enforced by constraints, not by which optimum a solver happens to pick. Genuine
        dual-targeting (a confined precursor that must be produced in two compartments) is recovered
        because both placements then carry flux and are rewarded. Cost: adds per movable reaction an
        activity binary per compartment (plus a reverse-direction binary for reversible reactions)
        and a home/used binary, so the MILP is materially larger at genome scale. The
        mono-localisation default is unchanged.

        **Residual caveat (internal cycles).** The activity threshold proves a placement carries
        ``≥ eps_flux`` somewhere in the biomass solution, but flux through an *internal/futile cycle*
        (a reaction running in a net-zero loop, e.g. a reversible reaction or one with a return leg)
        also clears the threshold without doing productive work — so such a placement could still
        harvest its score. Eliminating that needs loopless/thermodynamic (loop-law) constraints, the
        standard FBA internal-cycle problem, and is **deferred**. It is a far narrower gap than the
        zero-flux exploit it replaces (which scored placements carrying *no* flux at all).
        See ``docs/multi_localization_design.md``.
    eps_flux:
        Activity threshold for ``multi_localization``: an extra placement counts as *active* (and so
        may earn its score) only if its flux magnitude is ``≥ eps_flux`` in the biomass solution.
        ``None`` (default) sets it to ``10 * big_m * integrality_tol`` (``1e-5`` with the defaults) —
        just above the ghost-flux floor a tolerance-rounded binary could fake, and small enough to
        keep the activity *deadzone* ``(0, eps_flux)`` below the meaningful fluxes of genome-scale
        models. A larger fixed value (e.g. ``1e-4``) can make such models **infeasible**: a reaction
        forced to carry a tiny but nonzero flux in a compartment cannot then clear the threshold and,
        being *used*, is allowed no inactive home. Must exceed ``big_m * integrality_tol``. Ignored
        when ``multi_localization`` is ``False``.

    Returns
    -------
    AssignmentProposal
    """
    ghost_floor = big_m * integrality_tol
    if eps_flux is None:
        # Sit just above the ghost-flux floor (where a tolerance-rounded binary could fake activity)
        # — a 10x margin. Keeping eps_flux this small minimises the activity "deadzone" (0, eps): a
        # used reaction that must carry a tiny but nonzero flux in a compartment needs that flux to
        # clear eps, so an eps far above the floor (e.g. a fixed 1e-4) makes genome-scale models with
        # small fluxes — like yeast-GEM, biomass ~0.08 — infeasible. 10x the floor (1e-5 by default)
        # stays sound and feasible across yeast-GEM relocate sets.
        eps_flux = 10.0 * ghost_floor
    if multi_localization and eps_flux <= ghost_floor:
        raise ValueError(
            f"eps_flux={eps_flux} must exceed the ghost-flux floor big_m*integrality_tol="
            f"{ghost_floor}; otherwise a tolerance-rounded binary can fake activity."
        )
    compartments = sorted(set(model.compartments) | set(scores.compartments))
    if default_compartment not in compartments:
        raise ValueError(f"default_compartment={default_compartment!r} not in {compartments}")

    # ---- biomass reaction + growth floor --------------------------------------
    if biomass_reaction is None:
        obj = linear_reaction_coefficients(model)
        if not obj:
            raise ValueError("model has no objective; pass biomass_reaction=...")
        biomass_reaction = max(obj, key=lambda r: obj[r]).id
    if min_growth is None:
        base = model.slim_optimize(error_value=0.0)
        if not base or base <= 0:
            raise ValueError(
                "the draft model does not grow; cannot set a growth floor (pass min_growth=...)."
            )
        min_growth = 0.1 * base

    # ---- scope: movable vs pinned ---------------------------------------------
    to_relocate = set(reactions_to_relocate)
    movable: list[cobra.Reaction] = []
    for rid in sorted(to_relocate):
        r = model.reactions.get_by_id(rid)
        if r.boundary or _rxn_compartment(r) is None or r.id == biomass_reaction:
            continue
        movable.append(r)
    movable_ids = {r.id for r in movable}
    pinned = [r for r in model.reactions if r.id not in movable_ids]

    # ---- gap-fill candidates from the universal model -------------------------
    candidates: list[cobra.Reaction] = []
    if universal is not None:
        existing = set(model.reactions.list_attr("id"))
        candidates = [r for r in universal.reactions if r.id not in existing]

    def _gcost(rid: str) -> float:
        return float(gapfill_cost) if isinstance(gapfill_cost, (int, float)) \
            else float(gapfill_cost.get(rid, 1.0))

    # ---- genes in scope + score lookup ----------------------------------------
    score_df = scores.df
    genes_in_scope: set[str] = set()
    unplaced: list[str] = []
    for r in movable:
        gs = [g.id for g in r.genes]
        scored = [g for g in gs if g in score_df.index]
        if gs and not scored:
            unplaced.append(r.id)
        genes_in_scope.update(scored)

    def _score(g: str, c: str) -> float:
        if c not in score_df.columns or g not in score_df.index:
            return 0.0
        v = score_df.at[g, c]
        return 0.0 if v is None or (isinstance(v, float) and math.isnan(v)) else float(v)

    # ---- transportable base metabolites touched by movable reactions ----------
    base = base_metabolite if base_metabolite is not None else _base_met
    movable_base = {base(m) for r in movable for m in r.metabolites}
    transp = set(movable_base) if transportable is None else set(transportable) & movable_base

    # Stable integer index per base key. SBML ids are name-safe, but base keys may be metabolite
    # *names* (e.g. yeast-GEM) containing spaces/symbols that optlang forbids in variable and
    # constraint names — so names are built from this index, never the raw key.
    all_bases = {base(m) for r in (movable + pinned + candidates) for m in r.metabolites}
    bidx = {b: i for i, b in enumerate(sorted(all_bases))}

    def _tcost(b: str) -> float:
        return float(transport_cost) if isinstance(transport_cost, (int, float)) \
            else float(transport_cost.get(b, 0.5))

    # ---- build the MILP -------------------------------------------------------
    model.solver  # noqa: B018 — initialise the solver so model.problem is usable
    prob = model.problem
    opt = prob.Model()

    x = {(r.id, c): prob.Variable(f"x_{r.id}_{c}", type="binary")
         for r in movable for c in compartments}
    y = {(g, c): prob.Variable(f"y_{g}_{c}", type="binary")
         for g in genes_in_scope for c in compartments}
    t = {(b, c): prob.Variable(f"t_{bidx[b]}_{c}", type="binary")
         for b in transp for c in compartments if c != default_compartment}

    fmove = {(r.id, c): prob.Variable(
        f"vm_{r.id}_{c}", lb=min(r.lower_bound, 0.0), ub=max(r.upper_bound, 0.0))
        for r in movable for c in compartments}
    fpin = {r.id: prob.Variable(f"vp_{r.id}", lb=r.lower_bound, ub=r.upper_bound)
            for r in pinned}
    ftr = {(b, c): prob.Variable(f"vt_{bidx[b]}_{c}", lb=-big_m, ub=big_m) for (b, c) in t}
    # gap-fill candidate fluxes (vu) gated by z binaries
    vu = {r.id: prob.Variable(
        f"vu_{r.id}", lb=min(r.lower_bound, 0.0), ub=max(r.upper_bound, 0.0))
        for r in candidates}
    z = {r.id: prob.Variable(f"z_{r.id}", type="binary") for r in candidates}

    # ---- multi-localisation flux-activity coupling (opt-in) -------------------
    # A reaction earns a compartment's score only where it is *active* (carries >= eps_flux). Each
    # reaction may keep one inactive "home" (scored by prior, like a mono reaction not needed for
    # biomass), but only while it carries no flux anywhere; once it is used, every placement must be
    # active. This makes multi-localisation sound (no dead duplicate can harvest score) without
    # depending on which optimum the solver returns. aF/aR are the forward/reverse activity
    # binaries, h the home binary, u whether the reaction is used at all.
    aF: dict[tuple[str, str], object] = {}
    aR: dict[tuple[str, str], object] = {}
    h: dict[tuple[str, str], object] = {}
    used: dict[str, object] = {}
    if multi_localization:
        for r in movable:
            used[r.id] = prob.Variable(f"u_{r.id}", type="binary")
            for c in compartments:
                h[r.id, c] = prob.Variable(f"h_{r.id}_{c}", type="binary")
                if r.upper_bound > 0:
                    aF[r.id, c] = prob.Variable(f"aF_{r.id}_{c}", type="binary")
                if r.lower_bound < 0:
                    aR[r.id, c] = prob.Variable(f"aR_{r.id}_{c}", type="binary")

    def _active(rid, c):
        """Activity terms (forward + reverse) for placement (rid, c)."""
        return [v for v in (aF.get((rid, c)), aR.get((rid, c))) if v is not None]

    cons: list = []

    # mass balance per (base_metabolite, compartment) node: Σ coeff·flux = 0
    node_terms: dict[tuple[str, str], list] = {}

    def _add(node, var, coeff):
        node_terms.setdefault(node, []).append((float(coeff), var))

    for r in movable:
        for c in compartments:
            for m, coeff in r.metabolites.items():
                _add((base(m), c), fmove[r.id, c], coeff)
    for r in pinned:
        for m, coeff in r.metabolites.items():
            _add((base(m), m.compartment), fpin[r.id], coeff)
    for r in candidates:
        for m, coeff in r.metabolites.items():
            _add((base(m), m.compartment), vu[r.id], coeff)
    for (b, c), v in ftr.items():
        _add((b, default_compartment), v, -1.0)
        _add((b, c), v, 1.0)

    for node, terms in node_terms.items():
        cons.append(prob.Constraint(
            add([mul([Real(co), va]) for co, va in terms]),
            lb=0.0, ub=0.0, name=f"bal_{bidx[node[0]]}_{node[1]}"))

    # flux gating by placement binaries: lb·x ≤ v ≤ ub·x
    for r in movable:
        for c in compartments:
            xv, fv = x[r.id, c], fmove[r.id, c]
            cons.append(prob.Constraint(fv - mul([Real(r.upper_bound), xv]),
                                        ub=0.0, name=f"ub_{r.id}_{c}"))
            cons.append(prob.Constraint(fv - mul([Real(r.lower_bound), xv]),
                                        lb=0.0, name=f"lb_{r.id}_{c}"))
    for (b, c), v in ftr.items():
        tv = t[b, c]
        cons.append(prob.Constraint(v - mul([Real(big_m), tv]), ub=0.0, name=f"tub_{bidx[b]}_{c}"))
        cons.append(prob.Constraint(v + mul([Real(big_m), tv]), lb=0.0, name=f"tlb_{bidx[b]}_{c}"))
    # gap-fill candidate flux gating: lb·z ≤ vu ≤ ub·z (flux only if the candidate is added)
    for r in candidates:
        zv, fv = z[r.id], vu[r.id]
        cons.append(prob.Constraint(fv - mul([Real(r.upper_bound), zv]),
                                    ub=0.0, name=f"zub_{r.id}"))
        cons.append(prob.Constraint(fv - mul([Real(r.lower_bound), zv]),
                                    lb=0.0, name=f"zlb_{r.id}"))

    # assignment + gene coupling. By default each reaction goes to exactly one compartment
    # (mono-localisation): this prevents a reaction being "placed" in a high-score compartment
    # where it carries no flux just to harvest the score. With ``multi_localization=True`` a
    # reaction may occupy several compartments (Σx ≥ 1); the flux-activity coupling below then keeps
    # every extra placement honest (it must carry flux to earn its score).
    place_ub = None if multi_localization else 1.0
    for r in movable:
        cons.append(prob.Constraint(add([x[r.id, c] for c in compartments]),
                                    lb=1.0, ub=place_ub, name=f"place_{r.id}"))
        for g in {gg.id for gg in r.genes}:
            if g not in genes_in_scope:
                continue
            for c in compartments:
                cons.append(prob.Constraint(x[r.id, c] - y[g, c], ub=0.0,
                                            name=f"couple_{r.id}_{g}_{c}"))
    for g in genes_in_scope:
        cons.append(prob.Constraint(add([y[g, c] for c in compartments]),
                                    lb=1.0, name=f"gene1_{g}"))

    # A gene is "in" a compartment only if it catalyses a movable reaction placed there.
    # Without this, a gene could collect a compartment's localisation score for free (no
    # reaction there), which breaks the score/placement link and inflates multi-localisation.
    gene_rxns: dict[str, list[str]] = {g: [] for g in genes_in_scope}
    for r in movable:
        for g in {gg.id for gg in r.genes} & genes_in_scope:
            gene_rxns[g].append(r.id)
    for g in genes_in_scope:
        for c in compartments:
            expr = add([x[rid, c] for rid in gene_rxns[g]]) if gene_rxns[g] else Real(0.0)
            cons.append(prob.Constraint(y[g, c] - expr, ub=0.0, name=f"has_{g}_{c}"))

    # ---- multi-localisation flux-activity coupling -----------------------------
    # Make multi-localisation sound: a placement that carries no flux ("dead") cannot keep its
    # localisation score. For each movable reaction r and compartment c:
    #   * aF/aR are activity binaries that can be 1 only if the flux reaches +/- eps_flux
    #       eps*aF <= v[r,c]      (forward active => v >= eps)
    #       eps*aR <= -v[r,c]     (reverse active => v <= -eps)
    #   * a placement is materialised only as a home or an active placement
    #       x[r,c] <= h[r,c] + aF[r,c] + aR[r,c]      and  h,aF,aR <= x
    #   * the reaction is "used" iff active somewhere; the inactive home exists iff it is NOT used
    #       aF[r,c],aR[r,c] <= used[r] ;  used[r] <= Σ_c (aF+aR) ;  Σ_c h[r,c] == 1 - used[r]
    # Together: a used reaction has no inactive placement (every placement carries flux), an unused
    # reaction keeps exactly one inactive home (scored by prior, like a mono reaction). So a reaction
    # can never collect an extra compartment's score through a zero-flux duplicate — independent of
    # the solver. ``place_`` (Σx >= 1) keeps every reaction materialised somewhere.
    if multi_localization:
        for r in movable:
            # fmove[r,c] bounds (signed): lbm <= v <= ubm.
            lbm, ubm = min(r.lower_bound, 0.0), max(r.upper_bound, 0.0)
            acts = [v for c in compartments for v in _active(r.id, c)]
            for c in compartments:
                av = _active(r.id, c)
                # presence => home or active
                cons.append(prob.Constraint(x[r.id, c] - add([h[r.id, c]] + av),
                                            ub=0.0, name=f"present_{r.id}_{c}"))
                cons.append(prob.Constraint(h[r.id, c] - x[r.id, c], ub=0.0,
                                            name=f"home_le_x_{r.id}_{c}"))
                if (r.id, c) in aF:
                    cons.append(prob.Constraint(aF[r.id, c] - x[r.id, c], ub=0.0,
                                                name=f"aF_le_x_{r.id}_{c}"))
                    # aF=1 => v >= eps; aF=0 => v >= lbm (vacuous). A plain v >= eps*aF would
                    # force v >= 0 even when aF=0, wrongly forbidding reverse flux on a reversible
                    # reaction — so use the Big-M implication v - (eps-lbm)*aF >= lbm.
                    cons.append(prob.Constraint(
                        fmove[r.id, c] - mul([Real(eps_flux - lbm), aF[r.id, c]]),
                        lb=lbm, name=f"actF_{r.id}_{c}"))
                    cons.append(prob.Constraint(aF[r.id, c] - used[r.id], ub=0.0,
                                                name=f"usedF_{r.id}_{c}"))
                if (r.id, c) in aR:
                    cons.append(prob.Constraint(aR[r.id, c] - x[r.id, c], ub=0.0,
                                                name=f"aR_le_x_{r.id}_{c}"))
                    # aR=1 => v <= -eps; aR=0 => v <= ubm (vacuous): v + (ubm+eps)*aR <= ubm.
                    cons.append(prob.Constraint(
                        fmove[r.id, c] + mul([Real(ubm + eps_flux), aR[r.id, c]]),
                        ub=ubm, name=f"actR_{r.id}_{c}"))
                    cons.append(prob.Constraint(aR[r.id, c] - used[r.id], ub=0.0,
                                                name=f"usedR_{r.id}_{c}"))
            # used <= Σ activity  (used only if active somewhere)
            cons.append(prob.Constraint(used[r.id] - add(acts) if acts else used[r.id],
                                        ub=0.0, name=f"used_{r.id}"))
            # Σ home == 1 - used  (exactly one inactive home iff the reaction is unused)
            cons.append(prob.Constraint(add([h[r.id, c] for c in compartments]) + used[r.id],
                                        lb=1.0, ub=1.0, name=f"homecount_{r.id}"))

    # growth floor
    cons.append(prob.Constraint(mul([Real(1.0), fpin[biomass_reaction]]),
                                lb=min_growth, name="growth"))

    opt.add(list(x.values()) + list(y.values()) + list(t.values())
            + list(fmove.values()) + list(fpin.values()) + list(ftr.values())
            + list(vu.values()) + list(z.values())
            + list(aF.values()) + list(aR.values()) + list(h.values()) + list(used.values())
            + cons)

    # objective: maximise score − multi-penalty − transport cost − gap-fill cost
    obj_terms = []
    for g in genes_in_scope:
        for c in compartments:
            s = _score(g, c)
            if s:
                obj_terms.append(mul([Real(s), y[g, c]]))
    offset = 0.0
    if multi_compartment_penalty:
        for v in y.values():
            obj_terms.append(mul([Real(-multi_compartment_penalty), v]))
        offset = multi_compartment_penalty * len(genes_in_scope)
    for (b, _c), tv in t.items():
        obj_terms.append(mul([Real(-_tcost(b)), tv]))
    for rid, zv in z.items():
        obj_terms.append(mul([Real(-_gcost(rid)), zv]))

    opt.objective = prob.Objective(add(obj_terms) if obj_terms else Real(0.0), direction="max")
    if time_limit is not None:
        opt.configuration.timeout = int(time_limit)
    if mip_gap is not None:
        try:
            opt.problem.Params.MIPGap = mip_gap
        except Exception:  # noqa: BLE001
            pass

    # The Big-M flux gating (|v[r,c]| <= ub*x[r,c]) is only as tight as the solver's integer
    # tolerance: a placement binary rounded to ~tol still lets up to ub*tol of "ghost" flux pass
    # through a reaction placed elsewhere. With ub=1000 and a default tolerance of 1e-5 that is
    # 0.01 — enough to spuriously satisfy a small growth floor and certify a non-functional
    # placement. Tighten the integer-feasibility tolerance so the leak is negligible.
    _tighten_integrality(opt, integrality_tol)

    opt.optimize()
    status = opt.status

    proposal = AssignmentProposal(min_growth=min_growth, status=status,
                                  unplaced_reactions=unplaced)
    if status not in ("optimal", "feasible", "suboptimal", "time_limit"):
        return proposal

    for r in movable:
        proposal.placements[r.id] = [c for c in compartments
                                     if (x[r.id, c].primal or 0.0) >= 0.5]
    proposal.gene_compartments = {
        g: [c for c in compartments if (y[g, c].primal or 0.0) >= 0.5]
        for g in genes_in_scope}
    proposal.added_transports = [(b, c) for (b, c), tv in t.items()
                                 if (tv.primal or 0.0) >= 0.5]
    proposal.added_reactions = [rid for rid, zv in z.items()
                                if (zv.primal or 0.0) >= 0.5]
    proposal.objective = float((opt.objective.value or 0.0) + offset)
    return proposal


# --------------------------------------------------------------------------- apply
def apply_assignment(
    model: cobra.Model, proposal: AssignmentProposal, *, default_compartment: str = "c",
    base_metabolite: Callable[[cobra.Metabolite], str] | None = None,
    universal: cobra.Model | None = None,
) -> cobra.Model:
    """Build the compartmentalised model from a proposal (deep copy; original untouched).

    Pass the same ``base_metabolite`` and ``universal`` used in :func:`assign_compartments` so the
    same compartment-agnostic keying is used (existing per-compartment metabolites are reused
    rather than duplicated) and the chosen gap-fill reactions are added.
    """
    out = model.copy()
    base = base_metabolite if base_metabolite is not None else _base_met

    # (base key, compartment) -> existing metabolite, so a reaction moved into a compartment
    # reuses the species already there instead of creating a parallel duplicate.
    index: dict[tuple[str, str], cobra.Metabolite] = {}
    for m in out.metabolites:
        index.setdefault((base(m), m.compartment), m)

    def resolve(template: cobra.Metabolite, compartment: str) -> cobra.Metabolite:
        key = (base(template), compartment)
        met = index.get(key)
        if met is None:
            new_id = f"{template.id}__{compartment}"
            met = (out.metabolites.get_by_id(new_id) if new_id in out.metabolites
                   else cobra.Metabolite(new_id, name=template.name, compartment=compartment,
                                         formula=template.formula, charge=template.charge))
            if new_id not in out.metabolites:
                out.add_metabolites([met])
            index[key] = met
        return met

    for rid, comps in proposal.placements.items():
        if not comps:
            continue
        rxn = out.reactions.get_by_id(rid)
        _move_reaction(rxn, comps[0], resolve)
        for extra in comps[1:]:
            _duplicate_reaction(out, rxn, extra, resolve)
    for n, (b, c) in enumerate(proposal.added_transports):
        _add_transport(out, b, c, default_compartment, index, resolve, n)
    if universal is not None:
        for rid in proposal.added_reactions:
            _add_universal_reaction(out, universal.reactions.get_by_id(rid))
    return out


def _add_universal_reaction(model, urxn):
    if urxn.id in model.reactions:
        return
    new = cobra.Reaction(urxn.id, name=urxn.name,
                         lower_bound=urxn.lower_bound, upper_bound=urxn.upper_bound)
    model.add_reactions([new])
    mets = {}
    for m, coeff in urxn.metabolites.items():
        if m.id in model.metabolites:
            mm = model.metabolites.get_by_id(m.id)
        else:
            mm = cobra.Metabolite(m.id, name=m.name, compartment=m.compartment,
                                  formula=m.formula, charge=m.charge)
        mets[mm] = coeff
    new.add_metabolites(mets)
    new.gene_reaction_rule = urxn.gene_reaction_rule


def _move_reaction(rxn, compartment, resolve):
    old = list(rxn.metabolites.items())
    if all(m.compartment == compartment for m, _ in old):
        return
    rxn.subtract_metabolites(dict(old))
    rxn.add_metabolites({resolve(m, compartment): coeff for m, coeff in old})


def _duplicate_reaction(model, rxn, compartment, resolve):
    new_id = f"{rxn.id}_{compartment}"
    if new_id in model.reactions:
        return
    dup = cobra.Reaction(new_id, name=rxn.name,
                         lower_bound=rxn.lower_bound, upper_bound=rxn.upper_bound)
    model.add_reactions([dup])
    dup.add_metabolites({resolve(m, compartment): coeff for m, coeff in rxn.metabolites.items()})
    dup.gene_reaction_rule = rxn.gene_reaction_rule


def _add_transport(model, base_key, compartment, default_compartment, index, resolve, n):
    src = index.get((base_key, default_compartment))
    if src is None:
        return
    tr_id = f"tr_{n}_{compartment}"
    if tr_id in model.reactions:
        return
    dest = resolve(src, compartment)
    tr = cobra.Reaction(tr_id, name=f"transport {src.name} ({default_compartment}<->{compartment})",
                        lower_bound=-1000, upper_bound=1000)
    tr.add_metabolites({src: -1.0, dest: 1.0})
    model.add_reactions([tr])
