"""Add reactions to a model from equation strings.

Most of the equivalent MATLAB code is struct-of-arrays bookkeeping (padding parallel
``rxnNames`` / ``lb`` / ``ub`` / ``grRules`` / ... fields) that does not exist in
cobra, where each ``Reaction`` carries its own attributes. cobra also already
covers a large part of the *behaviour*:

* ``Reaction.build_reaction_from_string`` parses equation strings, coefficients,
  and reversibility arrows (``<=>``, ``-->``, ``=>``) and creates unknown
  metabolites — but only matching metabolites **by ID**, and it leaves new
  metabolites with ``compartment=None``.
* assigning ``reaction.gene_reaction_rule`` auto-creates ``Gene`` objects.

So this port keeps only the parts cobra lacks:

* **name-based matching** — interpret equation tokens as metabolite *names*
  (RAVEN eqnType 2) or as ``name[comp]`` (eqnType 3), not just IDs;
* **correct compartment** assignment for newly created metabolites;
* **strict policies** — optionally *error* (rather than silently create) on
  unknown metabolites or genes, and always error on a duplicate reaction ID
  (cobra silently ignores those).

Instead of RAVEN's ``eqnType`` integer (1/2/3) the matching mode is a readable
keyword: ``mets_by="id"`` or ``mets_by="name"``, with ``name[comp]`` recognised
automatically. See IMPROVEMENTS.md (A-series) for the rationale.
"""
from __future__ import annotations

import re
import warnings
from collections import OrderedDict
from collections.abc import Mapping, Sequence

import cobra
from cobra import Metabolite, Reaction
from cobra.core.gene import GPR

from raven_python.utils.parse import parse_name_comp, subsystem_to_str

# Reversibility arrows. ``<=>`` must be tried before ``=>`` (it contains it).
_REVERSIBLE_ARROWS = ("<=>",)
_FORWARD_ARROWS = ("-->", "->", "=>")


def _split_equation(equation: str) -> tuple[str, str, bool]:
    """Split an equation into (lhs, rhs, reversible) on its arrow."""
    for arrow in _REVERSIBLE_ARROWS:
        if arrow in equation:
            lhs, rhs = equation.split(arrow, 1)
            return lhs, rhs, True
    for arrow in _FORWARD_ARROWS:
        if arrow in equation:
            lhs, rhs = equation.split(arrow, 1)
            return lhs, rhs, False
    raise ValueError(f"No reaction arrow (<=>, -->, =>) found in equation: {equation!r}")


def _parse_side(side: str) -> list[tuple[float, str, str | None]]:
    """Parse one side of an equation into ``[(coefficient, token, fallback), ...]``.

    The ``fallback`` slot is for the ambiguous ``"<number> <rest>"`` shape: when
    matching by name, ``"2 oxoglutarate"`` could be either ``coeff=2, name="oxoglutarate"``
    or ``coeff=1, name="2 oxoglutarate"`` (a real chemistry name). We return the
    coefficient-split form as the primary and the full term as the fallback; the
    resolver picks whichever matches an existing metabolite. Pure-number heads
    with no name (``"2"``) and pure-name terms (``"glucose"``) have no fallback.
    """
    terms: list[tuple[float, str, str | None]] = []
    for raw in side.split(" + "):
        term = raw.strip()
        if not term:
            continue
        head, _, tail = term.partition(" ")
        try:
            coeff = float(head)
            token = tail.strip()
        except ValueError:
            coeff, token = 1.0, term
            fallback = None
        else:
            # Coefficient-split succeeded. Keep the full term as a fallback when
            # the tail is non-empty so name-resolution can re-try it as one token.
            fallback = term if token else None
        if not token:
            raise ValueError(f"Missing metabolite after coefficient in term: {raw!r}")
        terms.append((coeff, token, fallback))
    return terms


def _new_met_id(model: cobra.Model, prefix: str) -> str:
    """Next free ``<prefix><int>`` metabolite ID (RAVEN m1, m2, ... scheme)."""
    pattern = re.compile(rf"^{re.escape(prefix)}(\d+)$")
    used = [int(m.group(1)) for met in model.metabolites if (m := pattern.match(met.id))]
    n = max(used) + 1 if used else 1
    while f"{prefix}{n}" in model.metabolites:
        n += 1
    return f"{prefix}{n}"


def _try_existing(
    model: cobra.Model, token: str, *, mets_by: str, compartment: str | None
) -> Metabolite | None:
    """Look up ``token`` as an existing metabolite (no creation, no side effects).

    Returns the matching metabolite or ``None``. Used by ``_stoichiometry`` to
    disambiguate the ``"<number> <rest>"`` shape: if a metabolite whose *name*
    (or id) literally contains a leading number exists, prefer it over splitting
    the number off as a coefficient.
    """
    name, comp = parse_name_comp(token)
    if mets_by == "id" and comp is None:
        return model.metabolites.get_by_id(token) if token in model.metabolites else None
    target_comp = comp if comp is not None else compartment
    if target_comp is None:
        return None
    for met in model.metabolites:
        if met.name == name and met.compartment == target_comp:
            return met
    return None


def _resolve_metabolite(
    model: cobra.Model,
    token: str,
    *,
    mets_by: str,
    compartment: str | None,
    allow_new_mets: bool,
    new_met_prefix: str,
) -> Metabolite:
    """Resolve an equation token to an existing or newly created Metabolite."""
    name, comp = parse_name_comp(token)

    if mets_by == "id" and comp is None:
        # token is a metabolite ID
        if token in model.metabolites:
            return model.metabolites.get_by_id(token)
        if not allow_new_mets:
            raise ValueError(
                f"Unknown metabolite ID {token!r}; pass allow_new_mets=True to create it."
            )
        if compartment is None:
            raise ValueError(
                f"Cannot create metabolite {token!r}: no compartment given."
            )
        _warn_unknown_compartment(model, compartment, token)
        met = Metabolite(token, compartment=compartment)
        model.add_metabolites([met])
        return met

    # name-based (mets_by="name") or explicit name[comp]
    target_comp = comp if comp is not None else compartment
    if target_comp is None:
        raise ValueError(
            f"Metabolite {token!r} matched by name needs a compartment; "
            "pass compartment=... or use the name[comp] syntax."
        )
    if comp is not None and target_comp not in model.compartments and not allow_new_mets:
        raise ValueError(f"Compartment {target_comp!r} is not in the model.")

    matches = [
        met
        for met in model.metabolites
        if met.name == name and met.compartment == target_comp
    ]
    if matches:
        return matches[0]
    if not allow_new_mets:
        raise ValueError(
            f"No metabolite named {name!r} in compartment {target_comp!r}; "
            "pass allow_new_mets=True to create it."
        )
    _warn_unknown_compartment(model, target_comp, name)
    met = Metabolite(_new_met_id(model, new_met_prefix), name=name, compartment=target_comp)
    model.add_metabolites([met])
    return met


def _warn_unknown_compartment(model: cobra.Model, compartment: str, identifier: str) -> None:
    """Warn when a new metabolite would be born into a not-yet-registered compartment.

    Both ``mets_by`` paths previously created the metabolite without validating
    the compartment, so a typo (``"cyto"`` for ``"c"``) silently produced a
    one-metabolite ghost compartment. cobra inherits the compartment from the
    first metabolite assigned to it, so the fix is a warning, not a hard error.
    """
    known = set(model.compartments) | set(model._compartments)
    if compartment not in known:
        warnings.warn(
            f"Creating metabolite {identifier!r} in unregistered compartment "
            f"{compartment!r} (existing: {sorted(known) or 'none'}); "
            "add the compartment first or check for a typo.",
            stacklevel=5,
        )


def _stoichiometry(
    model: cobra.Model,
    equation: str,
    *,
    mets_by: str,
    compartment: str | None,
    allow_new_mets: bool,
    new_met_prefix: str,
) -> tuple[dict[Metabolite, float], bool]:
    """Parse an equation into a {Metabolite: net coefficient} dict + reversibility."""
    lhs, rhs, reversible = _split_equation(equation)
    coeffs: OrderedDict[Metabolite, float] = OrderedDict()
    had_terms = False
    for sign, side in ((-1.0, lhs), (1.0, rhs)):
        for coeff, token, fallback in _parse_side(side):
            had_terms = True
            # "<number> <name>" is ambiguous when the name itself starts with a
            # number (e.g. "2 oxoglutarate"). Prefer the full-term interpretation
            # when it matches an existing metabolite — otherwise fall through to
            # the coefficient-split form.
            met = None
            if fallback is not None:
                met = _try_existing(
                    model, fallback, mets_by=mets_by, compartment=compartment
                )
                if met is not None:
                    coeff = 1.0
            if met is None:
                met = _resolve_metabolite(
                    model,
                    token,
                    mets_by=mets_by,
                    compartment=compartment,
                    allow_new_mets=allow_new_mets,
                    new_met_prefix=new_met_prefix,
                )
            coeffs[met] = coeffs.get(met, 0.0) + sign * coeff
    # Drop metabolites that net to zero (present as both substrate and product).
    coeffs = OrderedDict((met, c) for met, c in coeffs.items() if c != 0.0)
    if had_terms and not coeffs:
        warnings.warn(
            f"Equation {equation!r} has no net metabolites (all terms cancelled); "
            "the reaction will be added with empty stoichiometry.",
            stacklevel=4,
        )
    return dict(coeffs), reversible


def add_reactions_from_equations(
    model: cobra.Model,
    reactions: Sequence[Mapping],
    *,
    mets_by: str = "id",
    compartment: str | None = None,
    allow_new_mets: bool = True,
    allow_new_genes: bool = True,
    new_met_prefix: str = "m",
) -> list[Reaction]:
    """Add reactions defined by equation strings, matching mets by ID or name.

    Parameters
    ----------
    model
        Target ``cobra.Model``, mutated in place.
    reactions
        Sequence of mappings, one per reaction. Recognised keys:

        * ``id`` (**required**) — reaction ID; must not already exist.
        * ``equation`` (**required**) — e.g. ``"atp_c + h2o_c <=> adp_c + pi_c"``.
          Use ``<=>`` for reversible, ``-->``/``->``/``=>`` for irreversible.
        * ``name`` — reaction name.
        * ``bounds`` — ``(lower, upper)`` tuple; overrides the arrow.
        * ``gene_reaction_rule`` — GPR string.
        * ``subsystem`` — subsystem name.
    mets_by
        How bare equation tokens (without ``[comp]``) are matched:
        ``"id"`` (RAVEN eqnType 1) or ``"name"`` (eqnType 2). A ``name[comp]``
        token (eqnType 3) is always matched by name + compartment.
    compartment
        Default compartment for new metabolites and for name-matched tokens
        without an explicit ``[comp]``.
    allow_new_mets
        If True (default), create metabolites not found. New metabolites get
        ``compartment`` (id mode) or an auto ID ``m1``, ``m2``, ... (name mode).
        If False, an unknown metabolite raises.
    allow_new_genes
        If True (default), genes in a GPR are auto-created by cobra. If False,
        a GPR referencing a gene not already in the model raises.
    new_met_prefix
        Prefix for auto-generated metabolite IDs in name mode (default ``"m"``).

    Returns
    -------
    list of cobra.Reaction
        The reactions added, in input order.
    """
    if mets_by not in ("id", "name"):
        raise ValueError(f"mets_by must be 'id' or 'name', got {mets_by!r}")

    known_genes = {gene.id for gene in model.genes}
    added: list[Reaction] = []

    for spec in reactions:
        if "id" not in spec:
            raise ValueError(f"Reaction spec missing required 'id': {spec!r}")
        rxn_id = spec["id"]
        if rxn_id in model.reactions:
            raise ValueError(
                f"Reaction {rxn_id!r} already exists. To change its stoichiometry use "
                "change_reaction_equations; to replace it, remove it first with "
                "model.remove_reactions([...])."
            )
        if "equation" not in spec:
            raise ValueError(f"Reaction {rxn_id!r} spec missing required 'equation'.")

        coeffs, reversible = _stoichiometry(
            model,
            spec["equation"],
            mets_by=mets_by,
            compartment=compartment,
            allow_new_mets=allow_new_mets,
            new_met_prefix=new_met_prefix,
        )

        rxn = Reaction(rxn_id, name=spec.get("name", ""))
        if "bounds" in spec:
            rxn.bounds = tuple(spec["bounds"])
        else:
            config = cobra.Configuration()
            lower = config.lower_bound if reversible else 0.0
            rxn.bounds = (lower, config.upper_bound)
        if "subsystem" in spec:
            rxn.subsystem = subsystem_to_str(spec["subsystem"])

        model.add_reactions([rxn])
        rxn.add_metabolites(coeffs)

        rule = spec.get("gene_reaction_rule", "")
        if rule:
            if not allow_new_genes:
                missing = sorted(set(GPR.from_string(rule).genes) - known_genes)
                if missing:
                    raise ValueError(
                        f"Reaction {rxn_id!r} references genes not in the model: "
                        f"{missing}. Set allow_new_genes=True or add them first."
                    )
            rxn.gene_reaction_rule = rule
            known_genes.update(gene.id for gene in rxn.genes)

        added.append(rxn)

    return added
