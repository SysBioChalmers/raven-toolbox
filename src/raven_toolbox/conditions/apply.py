"""Generic "apply this YAML condition" loader.

Schema (all keys optional):

.. code-block:: yaml

    name: my_condition
    description: free-form text.

    # Set exchange reactions to (0, ub) before any per-rxn override.
    # "in" / "out" (RAVEN's getExchangeRxns direction: the boundary metabolite
    # is the reaction's substrate / product respectively) reset only that
    # direction; "all" / "both", or any other truthy value, reset every
    # exchange in either direction.
    prelude:
      reset_exchanges: out                # "in" | "out" | "all" | "both" | bool truthy

    # Remove metabolites from a pseudoreaction and rebalance H+
    # (or any other "balance" met) so total charge stays zero.
    cofactor_pseudoreaction:
      rxn_id: r_4598
      remove_mets:
        - { met: s_3714 }                 # comment field is informational
      charge_balance_met: s_0794          # optional

    # Add (combine=True) to a reaction's stoichiometry.
    biomass_stoichiometry_delta:
      rxn_id: r_4041
      add:
        - { met: s_0689, coef:  0.08 }
        - { met: s_0687, coef: -0.08 }
        - { met: s_0794, coef: -0.16 }

    # Per-reaction bounds. lb / ub omitted means "leave unchanged".
    bounds:
      - { rxn: r_1654, lb: -1000 }
      - { rxn: r_1992, lb: 0 }
      - { rxn: r_1663, lb: 0, ub: 0 }

    # Sanity-check counter. Emits a warning if the actual count of
    # ``lb == -1000`` bound entries doesn't match.
    expected_uptake_count: 15

The schema is intentionally narrow: deterministic, idempotent, easy to
diff in code review. Project-specific extensions (e.g. yeast-GEM's
``amino_acid_ratio``) are handled by the *caller* before / after this
function — keep the upstream generic.
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import cobra
from ruamel.yaml import YAML

# A safe loader keeps the parsed document as plain dict / list / scalars,
# which matches what callers expect from ``load_condition``. ruamel.yaml
# is already a transitive dependency via cobra, so we don't take on
# PyYAML on top.
_SAFE_YAML = YAML(typ="safe")

#: ``prelude.reset_exchanges`` puts every exchange reaction at this
#: upper bound (and lower bound = 0).
DEFAULT_RESET_EXCHANGES_UPPER_BOUND: float = 1000.0


def load_condition(path: str | Path) -> dict[str, Any]:
    """Parse a condition YAML file into a plain dict."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Condition file not found: {path}")
    with open(path) as f:
        return _SAFE_YAML.load(f)


def apply_condition(
    model: cobra.Model,
    condition: dict[str, Any] | str | Path,
) -> cobra.Model:
    """Apply a parsed condition (or a YAML file path) to ``model`` in place.

    Returns the same model object for chaining.
    """
    cfg = condition if isinstance(condition, dict) else load_condition(condition)

    _apply_prelude(model, cfg)
    _apply_cofactor_pseudoreaction(model, cfg)
    _apply_biomass_stoichiometry_delta(model, cfg)
    n_uptake = _apply_bounds(model, cfg)
    _check_uptake_count(cfg, n_uptake)

    return model


def set_reaction_bounds(rxn: cobra.Reaction, lb: float, ub: float) -> None:
    """Set both bounds, bypassing cobra's ``lb <= ub`` validator.

    Some conditions intentionally land on an infeasible ``lb > ub`` state
    (e.g. forcing flux via a sentinel bound). cobra's setter rejects
    that. This helper writes the underlying private attributes when
    needed so the model can match the caller's exact intent.
    """
    if lb > ub:
        rxn._lower_bound = lb  # noqa: SLF001 — intentional bypass
        rxn._upper_bound = ub  # noqa: SLF001 — intentional bypass
    else:
        rxn.bounds = (lb, ub)


# --- step implementations ---------------------------------------------

def _exchange_direction(rxn: cobra.Reaction) -> str | None:
    """RAVEN ``getExchangeRxns``' in/out rule, direction-of-boundary-metabolite based.

    "out": nothing is produced within the model at all (``sum(S > 0) == 0``) --- the
    boundary metabolite is implicitly the reaction's product. "in": nothing is
    consumed (``sum(S < 0) == 0``) --- the boundary metabolite is implicitly the
    substrate. Neither (an ordinary internal reaction with both) returns ``None``.
    Not restricted to single-metabolite reactions, matching RAVEN's own rule exactly
    rather than cobra's narrower ``Reaction.boundary``.
    """
    has_product = any(c > 0 for c in rxn.metabolites.values())
    has_reactant = any(c < 0 for c in rxn.metabolites.values())
    if not has_product:
        return "out"
    if not has_reactant:
        return "in"
    return None


def _apply_prelude(model: cobra.Model, cfg: dict[str, Any]) -> None:
    prelude = cfg.get("prelude") or {}
    value = prelude.get("reset_exchanges")
    if not value:
        return
    direction = value.lower() if isinstance(value, str) else None
    if direction in ("in", "out"):
        # A named direction resets only that direction, matching RAVEN's
        # applyCondition, which forwards this same value straight through to
        # getExchangeRxns as its direction filter.
        targets = [r for r in model.reactions if _exchange_direction(r) == direction]
    else:
        # "all" / "both", or any other truthy value that isn't a direction
        # keyword (e.g. a bare `true`): reset every reaction RAVEN's
        # getExchangeRxns would call an exchange in either direction.
        targets = [r for r in model.reactions if _exchange_direction(r) is not None]
    for rxn in targets:
        rxn.lower_bound = 0
        rxn.upper_bound = DEFAULT_RESET_EXCHANGES_UPPER_BOUND


def _apply_cofactor_pseudoreaction(model: cobra.Model, cfg: dict[str, Any]) -> None:
    cp = cfg.get("cofactor_pseudoreaction")
    if not cp:
        return
    rxn = model.reactions.get_by_id(cp["rxn_id"])
    for entry in cp.get("remove_mets", []):
        met = model.metabolites.get_by_id(entry["met"])
        _set_coefficient(rxn, met, 0.0)
    balance_met_id = cp.get("charge_balance_met")
    if balance_met_id:
        balance_met = model.metabolites.get_by_id(balance_met_id)
        _set_coefficient(rxn, balance_met, 0.0)
        total_charge = sum(
            (m.charge or 0) * coef
            for m, coef in rxn.metabolites.items()
        )
        _set_coefficient(rxn, balance_met, -total_charge)


def _apply_biomass_stoichiometry_delta(model: cobra.Model, cfg: dict[str, Any]) -> None:
    delta = cfg.get("biomass_stoichiometry_delta")
    if not delta:
        return
    rxn = model.reactions.get_by_id(delta["rxn_id"])
    for entry in delta.get("add", []):
        met = model.metabolites.get_by_id(entry["met"])
        rxn.add_metabolites({met: float(entry["coef"])}, combine=True)


def _apply_bounds(model: cobra.Model, cfg: dict[str, Any]) -> int:
    n_uptake = 0
    for entry in cfg.get("bounds", []):
        try:
            rxn = model.reactions.get_by_id(entry["rxn"])
        except KeyError:
            warnings.warn(
                f"Reaction {entry['rxn']!r} not found in model; skipping.",
                stacklevel=3,
            )
            continue
        new_lb = float(entry["lb"]) if "lb" in entry else rxn.lower_bound
        new_ub = float(entry["ub"]) if "ub" in entry else rxn.upper_bound
        set_reaction_bounds(rxn, new_lb, new_ub)
        if entry.get("lb") == -1000:
            n_uptake += 1
    return n_uptake


def _check_uptake_count(cfg: dict[str, Any], n_uptake: int) -> None:
    expected = cfg.get("expected_uptake_count")
    if expected is None:
        return
    if n_uptake != expected:
        warnings.warn(
            f"Expected {expected} uptake reactions, applied {n_uptake}. "
            "Some referenced reactions may be missing from the model.",
            stacklevel=3,
        )


# --- helpers ----------------------------------------------------------

def _set_coefficient(rxn: cobra.Reaction, met: cobra.Metabolite, value: float) -> None:
    """Set the stoichiometric coefficient of ``met`` in ``rxn`` to ``value``.

    cobra's only public mutation point is ``add_metabolites``; we go via
    ``combine=True`` with ``(value - current)`` to land on the desired
    value without depending on cobra's removal-on-zero behaviour.
    """
    current = rxn.metabolites.get(met, 0.0)
    delta = float(value) - current
    if delta != 0:
        rxn.add_metabolites({met: delta}, combine=True)
