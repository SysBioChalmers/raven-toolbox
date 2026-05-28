"""Subcellular localisation by MILP (port of RAVEN ``predictLocalization``).

Assigns reactions to compartments by maximising per-gene localisation evidence minus
inter-compartment transport cost. **Critical differences from RAVEN's function**
(see [docs/localization_design.md](../../docs/localization_design.md)):

* Existing compartmentalisation is **respected by default** — only reactions explicitly
  flagged (``reactions.notes['localization']='uncertain'``) or passed in
  ``reactions_to_relocate`` are placed; everything else is pinned.
* Incomplete models are tolerated — no silent reaction removal for "metabolite not
  produced". Reactions with no scored genes are reported in ``unplaced_reactions``.
* Boundary reactions and existing inter-compartment transports are auto-pinned.
* A deterministic MILP (Gurobi / HiGHS / GLPK) replaces simulated annealing.
* ``apply=False`` returns a :class:`LocalizationProposal` (a diff) without mutating.
* Gene mono-localisation by default; ``multi_compartment_genes=True`` allows multiple
  compartments with a configurable per-extra-compartment penalty.

Limitations vs RAVEN to be aware of:

* Isozyme separation (RAVEN's ``expandModel`` step) is *not* applied — a reaction with
  isozymes is treated as "all listed genes must share its compartment". To get
  per-isozyme placement, call :func:`ravengem.manipulation.expand_model` first.
* Transports are routed through ``default_compartment`` (RAVEN's simplification).
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field

import cobra
import pandas as pd
from optlang.symbolics import Real, add, mul

from ravengem.localization.scores import LocalizationScores

_UNCERTAIN = "uncertain"


@dataclass
class LocalizationProposal:
    """What :func:`predict_localization` proposes, before applying it.

    All DataFrames have one row per item. Use this with ``apply=False`` to preview
    changes; pass it back to :func:`apply_localization` to commit, or diff against a
    curator's expectations.
    """

    moved: pd.DataFrame                       # rxn_id, from_compartment, to_compartment
    added_transports: pd.DataFrame            # met_id, compartment (other than default)
    gene_compartments: dict[str, list[str]]   # gene_id → list of compartments assigned
    unplaced_reactions: list[str] = field(default_factory=list)  # had no scored gene support
    objective: float = 0.0


@dataclass
class LocalizationResult:
    """Outcome of :func:`predict_localization` (when ``apply=True``)."""

    model: cobra.Model
    proposal: LocalizationProposal
    added_transports: list[cobra.Reaction] = field(default_factory=list)


# --------------------------------------------------------------------- helpers

def _reaction_compartment(rxn: cobra.Reaction) -> str | None:
    """Single compartment id if all metabolites share one, else ``None`` (transport)."""
    comps = {m.compartment for m in rxn.metabolites if m.compartment}
    return next(iter(comps)) if len(comps) == 1 else None


def _auto_relocate_set(model: cobra.Model) -> set[str]:
    """Default scope when the user passes no explicit list:

    a reaction is relocated iff ``notes['localization'] == 'uncertain'``. (Boundary
    reactions and reactions that already span multiple compartments are always pinned.)
    """
    out: set[str] = set()
    for r in model.reactions:
        if r.boundary or _reaction_compartment(r) is None:
            continue  # boundary / existing transport → pinned
        if (r.notes or {}).get("localization") == _UNCERTAIN:
            out.add(r.id)
    return out


def _reaction_genes(rxn: cobra.Reaction) -> list[str]:
    """Genes on the reaction's GPR (flat list; no AND/OR distinction in this v1)."""
    return [g.id for g in rxn.genes]


# --------------------------------------------------------------------- the MILP

def predict_localization(
    model: cobra.Model,
    scores: LocalizationScores,
    *,
    reactions_to_relocate: Iterable[str] | None = None,
    default_compartment: str = "c",
    transport_cost: float | Mapping[str, float] = 0.5,
    multi_compartment_genes: bool = False,
    multi_compartment_penalty: float = 0.5,
    apply: bool = True,
    mip_gap: float | None = None,
    time_limit: float | None = None,
) -> LocalizationResult | LocalizationProposal:
    """Place reactions in compartments via MILP. Returns a :class:`LocalizationProposal`
    (when ``apply=False``) or a :class:`LocalizationResult` (when ``apply=True``).

    ``reactions_to_relocate``: subset of reaction ids to (re-)place. When ``None``,
    reactions flagged ``notes['localization']='uncertain'`` are relocated and everything
    else is pinned. Boundary reactions and existing multi-compartment transports are
    always pinned.

    ``transport_cost``: either a scalar (same cost per added transport) or a mapping
    ``{metabolite_id_base: cost}`` (where the base id strips the compartment suffix,
    e.g. ``"glc__D"`` matches ``"glc__D_c"``/``"glc__D_e"``). Negative costs *favour*
    adding the transport.

    With ``multi_compartment_genes=True`` a gene can land in several compartments at
    once (e.g. mitochondrial/cytosolic isoforms); each extra compartment beyond the
    first costs ``multi_compartment_penalty`` in the objective.
    """
    # ---- 1. Scope: which reactions move, which are pinned. -----------------
    if reactions_to_relocate is None:
        to_relocate = _auto_relocate_set(model)
    else:
        to_relocate = set(reactions_to_relocate)
        # Filter out boundaries / transports that we always pin.
        to_relocate -= {r.id for r in model.reactions
                        if r.boundary or _reaction_compartment(r) is None}
    if not to_relocate:
        return _empty_result(model, apply)

    # ---- 2. Compartments universe (model + scores). ------------------------
    compartments = sorted(set(model.compartments) | set(scores.compartments))
    if default_compartment not in compartments:
        raise ValueError(f"default_compartment={default_compartment!r} not in known "
                         f"compartments {compartments}")

    # ---- 3. Gather genes for the relocate-set, build score lookup. ---------
    # Genes only mentioned by pinned reactions don't enter the MILP.
    moving = [model.reactions.get_by_id(rid) for rid in sorted(to_relocate)]
    genes_in_scope: set[str] = set()
    unplaced: list[str] = []
    for r in moving:
        gs = _reaction_genes(r)
        scored = [g for g in gs if g in scores.df.index]
        if not gs:
            # GPR-less reaction: place it freely (no gene coupling). Allowed.
            continue
        if not scored:
            # All genes absent from predictor → no signal; report and skip.
            unplaced.append(r.id)
            continue
        genes_in_scope.update(scored)
    # Remove reactions we can't score from the placement set.
    placeable = [r for r in moving if r.id not in set(unplaced)]
    if not placeable:
        # Everything in the relocate set lacks scored genes — return a proposal with
        # only the unplaced list.
        prop = LocalizationProposal(
            moved=pd.DataFrame(columns=["rxn_id", "from_compartment", "to_compartment"]),
            added_transports=pd.DataFrame(columns=["met_id", "compartment"]),
            gene_compartments={}, unplaced_reactions=unplaced, objective=0.0)
        return prop if not apply else LocalizationResult(model, prop)

    # ---- 4. Per-metabolite transport cost. ---------------------------------
    def _met_cost(m_id: str) -> float:
        if not isinstance(transport_cost, (int, float)):
            base = m_id.rsplit("_", 1)[0]
            return float(transport_cost.get(base, transport_cost.get(m_id, 0.5)))
        return float(transport_cost)

    # ---- 5. Build the MILP. ------------------------------------------------
    model.solver  # noqa: B018 — ensure the solver is initialised so model.problem works
    prob = model.problem
    opt = prob.Model()

    # x[r, c] = 1 iff reaction r placed in c (only for r ∈ placeable)
    x: dict[tuple[str, str], object] = {
        (r.id, c): prob.Variable(f"x_{r.id}_{c}", type="binary")
        for r in placeable for c in compartments
    }
    # y[g, c] = 1 iff gene g assigned to c
    y: dict[tuple[str, str], object] = {
        (g, c): prob.Variable(f"y_{g}_{c}", type="binary")
        for g in genes_in_scope for c in compartments
    }
    # t[m_id, c] = 1 iff metabolite m (with id including its current compartment suffix)
    # needs a transport to compartment c (c ≠ default). One per (base met, c).
    met_keys: set[tuple[str, str]] = set()
    for r in placeable:
        for m in r.metabolites:
            for c in compartments:
                if c != default_compartment:
                    met_keys.add((m.id, c))
    t: dict[tuple[str, str], object] = {
        k: prob.Variable(f"t_{k[0]}_{k[1]}", type="binary") for k in met_keys
    }

    cons: list = []
    # 5a. Each placeable reaction goes to exactly one compartment.
    for r in placeable:
        cons.append(prob.Constraint(add([mul([Real(1.0), x[r.id, c]]) for c in compartments]),
                                     lb=1.0, ub=1.0, name=f"place_{r.id}"))
    # 5b. Gene-reaction coupling: if r placed in c, every scored gene of r must be in c.
    for r in placeable:
        for g in _reaction_genes(r):
            if g not in genes_in_scope:
                continue
            for c in compartments:
                # x[r,c] − y[g,c] ≤ 0
                cons.append(prob.Constraint(x[r.id, c] - y[g, c], ub=0.0,
                                             name=f"gene_{r.id}_{g}_{c}"))
    # 5c. Gene assignment.
    for g in genes_in_scope:
        s = add([mul([Real(1.0), y[g, c]]) for c in compartments])
        if multi_compartment_genes:
            # ≥1 (gene must land somewhere if any of its reactions is placed); penalty handled
            # in the objective.
            cons.append(prob.Constraint(s, lb=1.0, name=f"gene_one_{g}"))
        else:
            cons.append(prob.Constraint(s, lb=1.0, ub=1.0, name=f"gene_one_{g}"))
    # 5d. Transport requirement: t[m,c] ≥ x[r,c] whenever r touches m and c ≠ default.
    for r in placeable:
        for m in r.metabolites:
            for c in compartments:
                if c == default_compartment:
                    continue
                # x[r,c] − t[m,c] ≤ 0
                cons.append(prob.Constraint(x[r.id, c] - t[m.id, c], ub=0.0,
                                             name=f"trans_{r.id}_{m.id}_{c}"))

    opt.add(list(x.values()) + list(y.values()) + list(t.values()) + cons)

    # 5e. Objective.
    obj_terms = []
    # + per-gene per-compartment localisation score (rows missing → 0)
    score_lookup = scores.df  # gene_id × compartment → float
    for g in genes_in_scope:
        for c in compartments:
            s = float(score_lookup.at[g, c]) if c in score_lookup.columns and not pd.isna(score_lookup.at[g, c]) else 0.0
            if s:
                obj_terms.append(mul([Real(s), y[g, c]]))
    # − transport cost per added transport
    for (m_id, _c), tvar in t.items():
        cost = _met_cost(m_id)
        if cost:
            obj_terms.append(mul([Real(-cost), tvar]))
    # − multi-compartment penalty per extra compartment (linearised: see design §2.3)
    if multi_compartment_genes and multi_compartment_penalty:
        for yvar in y.values():
            obj_terms.append(mul([Real(-multi_compartment_penalty), yvar]))
        # Net effect: penalty * (Σ_c y[g,c] − 0) — constant offset matters only if scores ≤ penalty;
        # design §2.3 explains the absorbed constant.

    opt.objective = prob.Objective(add(obj_terms) if obj_terms else Real(0.0), direction="max")
    if time_limit is not None:
        opt.configuration.timeout = int(time_limit)
    if mip_gap is not None:
        try:  # Gurobi-specific
            opt.problem.Params.MIPGap = mip_gap
        except Exception:  # noqa: BLE001
            pass

    opt.optimize()
    if opt.status not in ("optimal", "feasible", "suboptimal", "time_limit"):
        raise RuntimeError(f"localisation MILP did not solve (status: {opt.status}).")

    # ---- 6. Read the solution into a proposal. -----------------------------
    moved_rows: list[dict] = []
    for r in placeable:
        chosen = None
        for c in compartments:
            if (x[r.id, c].primal or 0.0) >= 0.5:
                chosen = c
                break
        from_c = _reaction_compartment(r)
        if chosen and chosen != from_c:
            moved_rows.append({"rxn_id": r.id, "from_compartment": from_c,
                                "to_compartment": chosen})
    moved = pd.DataFrame(moved_rows, columns=["rxn_id", "from_compartment", "to_compartment"])

    transp_rows: list[dict] = []
    for (m_id, c), tvar in t.items():
        if (tvar.primal or 0.0) >= 0.5:
            transp_rows.append({"met_id": m_id, "compartment": c})
    added_transports = pd.DataFrame(transp_rows, columns=["met_id", "compartment"])

    gene_comps: dict[str, list[str]] = {}
    for g in genes_in_scope:
        in_c = [c for c in compartments if (y[g, c].primal or 0.0) >= 0.5]
        gene_comps[g] = in_c

    proposal = LocalizationProposal(
        moved=moved, added_transports=added_transports, gene_compartments=gene_comps,
        unplaced_reactions=unplaced, objective=float(opt.objective.value or 0.0))

    if not apply:
        return proposal
    new_model, transports = apply_localization(model, proposal, default_compartment=default_compartment)
    return LocalizationResult(model=new_model, proposal=proposal, added_transports=transports)


# --------------------------------------------------------------------- apply

def apply_localization(
    model: cobra.Model,
    proposal: LocalizationProposal,
    *,
    default_compartment: str = "c",
) -> tuple[cobra.Model, list[cobra.Reaction]]:
    """Apply a :class:`LocalizationProposal` to ``model``: move reactions, add the
    inter-compartment transports the proposal listed, and return ``(model_copy, added)``.

    The returned model is a deep copy of the input (original left untouched). Moved
    reactions get their metabolites' compartment suffix swapped (e.g. ``A_c → A_m``);
    new compartment-specific metabolite copies are added on demand. Each added
    transport is a passive diffusion ``M[default] ⇌ M[c]`` (RAVEN convention),
    named ``tr_<met>_<c>``.
    """
    out = model.copy()
    added: list[cobra.Reaction] = []

    # 1. Move each reaction by remapping its metabolites to the target compartment.
    for _, row in proposal.moved.iterrows():
        rxn = out.reactions.get_by_id(row["rxn_id"])
        target = row["to_compartment"]
        new_stoich: dict[cobra.Metabolite, float] = {}
        old = list(rxn.metabolites.items())
        # Clear current stoichiometry first so cobra updates the constraints cleanly.
        rxn.subtract_metabolites(dict(old))
        for m, coeff in old:
            m_new = _met_in_compartment(out, m, target)
            new_stoich[m_new] = coeff
        rxn.add_metabolites(new_stoich)

    # 2. Add transports between default and each requested compartment.
    for _, row in proposal.added_transports.iterrows():
        m_id, c = row["met_id"], row["compartment"]
        if m_id not in out.metabolites:
            continue
        m_src = out.metabolites.get_by_id(m_id)
        if m_src.compartment == c:
            continue  # already there; no transport needed
        m_default = _met_in_compartment(out, m_src, default_compartment)
        m_dest = _met_in_compartment(out, m_src, c)
        if m_default.id == m_dest.id:
            continue
        tr_id = f"tr_{_base_met_id(m_src)}_{c}"
        if tr_id in out.reactions:
            continue
        tr = cobra.Reaction(tr_id, lower_bound=-1000, upper_bound=1000)
        tr.add_metabolites({m_default: -1.0, m_dest: 1.0})
        tr.notes["localization"] = "added by predict_localization"
        out.add_reactions([tr])
        added.append(out.reactions.get_by_id(tr_id))

    return out, added


def _base_met_id(m: cobra.Metabolite) -> str:
    """Strip the trailing ``_<compartment>`` suffix (or return id as-is)."""
    if m.compartment and m.id.endswith(f"_{m.compartment}"):
        return m.id[: -(len(m.compartment) + 1)]
    return m.id


def _met_in_compartment(model: cobra.Model, source: cobra.Metabolite,
                        compartment: str) -> cobra.Metabolite:
    """Return (creating if needed) the copy of ``source`` in ``compartment``."""
    if source.compartment == compartment:
        return source
    base = _base_met_id(source)
    new_id = f"{base}_{compartment}"
    if new_id in model.metabolites:
        return model.metabolites.get_by_id(new_id)
    new_met = cobra.Metabolite(new_id, name=source.name, compartment=compartment,
                               formula=source.formula, charge=source.charge)
    new_met.notes = dict(source.notes or {})
    model.add_metabolites([new_met])
    return new_met


def _empty_result(model: cobra.Model, apply_flag: bool):
    proposal = LocalizationProposal(
        moved=pd.DataFrame(columns=["rxn_id", "from_compartment", "to_compartment"]),
        added_transports=pd.DataFrame(columns=["met_id", "compartment"]),
        gene_compartments={}, unplaced_reactions=[], objective=0.0)
    return proposal if not apply_flag else LocalizationResult(model.copy(), proposal)
