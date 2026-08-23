"""Tests for convert_to_irreversible (RAVEN convertToIrrev.m).

Adopted from geckopy's tests/test_preprocess.py (the convert_to_irreversible subset).
Exchange reactions are excluded from the split, matching MATLAB behavior.
"""
import cobra

from raven_toolbox.manipulation import convert_to_irreversible


def _build_model_with_bounds(
    reactions: list[tuple[str, dict[str, float], float, float]],
) -> cobra.Model:
    """Build from (rxn_id, {met_id: coef}, lb, ub) tuples."""
    model = cobra.Model("test")
    mets: dict[str, cobra.Metabolite] = {}
    for _, stoich, _, _ in reactions:
        for met_id in stoich:
            if met_id not in mets:
                mets[met_id] = cobra.Metabolite(met_id, compartment="c")

    for rxn_id, stoich, lb, ub in reactions:
        rxn = cobra.Reaction(rxn_id)
        rxn.lower_bound = lb
        rxn.upper_bound = ub
        rxn.add_metabolites({mets[m]: c for m, c in stoich.items()})
        model.add_reactions([rxn])
    return model


def test_splits_single_reversible_non_exchange():
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),
    ])

    added = convert_to_irreversible(model)
    assert added == ["r1_REV"]

    fwd = model.reactions.get_by_id("r1")
    rev = model.reactions.get_by_id("r1_REV")

    assert fwd.bounds == (0.0, 1000.0)
    assert {m.id: c for m, c in fwd.metabolites.items()} == {"A": -1.0, "B": 1.0}

    assert rev.bounds == (0.0, 500.0)
    assert {m.id: c for m, c in rev.metabolites.items()} == {"A": 1.0, "B": -1.0}


def test_does_not_split_forward_only_reaction():
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0),
    ])
    added = convert_to_irreversible(model)
    assert added == []
    assert "r1_REV" not in {r.id for r in model.reactions}


def test_does_not_split_exchange_reaction_even_if_reversible():
    """Exchange reactions (one metabolite) are explicitly excluded from
    the irreversibility step in MATLAB, regardless of bounds."""
    model = _build_model_with_bounds([
        ("EX_A", {"A": -1.0}, -1000.0, 1000.0),
    ])
    added = convert_to_irreversible(model)
    assert added == []
    ex = model.reactions.get_by_id("EX_A")
    assert ex.bounds == (-1000.0, 1000.0)


def test_splits_multiple_mixed_reactions():
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),   # split
        ("r2", {"B": -2.0, "C": 3.0}, 0.0, 1000.0),      # forward only
        ("EX_A", {"A": -1.0}, -1000.0, 1000.0),          # exchange
        ("r3", {"C": -1.0, "D": 1.0}, -200.0, 200.0),    # split
    ])

    added = convert_to_irreversible(model)
    assert added == ["r1_REV", "r3_REV"]

    assert model.reactions.get_by_id("r1").bounds == (0.0, 1000.0)
    assert model.reactions.get_by_id("r1_REV").bounds == (0.0, 500.0)
    assert model.reactions.get_by_id("r2").bounds == (0.0, 1000.0)
    assert model.reactions.get_by_id("EX_A").bounds == (-1000.0, 1000.0)
    assert model.reactions.get_by_id("r3").bounds == (0.0, 200.0)
    assert model.reactions.get_by_id("r3_REV").bounds == (0.0, 200.0)


def test_reverse_reaction_inherits_gpr():
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),
    ])
    model.reactions.get_by_id("r1").gene_reaction_rule = "g1 and g2"

    convert_to_irreversible(model)

    rev = model.reactions.get_by_id("r1_REV")
    assert rev.gene_reaction_rule == "g1 and g2"
    assert {g.id for g in rev.genes} == {"g1", "g2"}


def test_forward_reaction_lb_is_clamped_to_zero():
    """After splitting, the original reaction should have lb = 0,
    which is what MATLAB's convertToIrrev does."""
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),
    ])
    convert_to_irreversible(model)
    assert model.reactions.get_by_id("r1").lower_bound == 0.0


def test_no_reverse_reaction_has_negative_bound():
    """After conversion, no non-exchange reaction may carry negative flux."""
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),
        ("r2", {"B": -1.0, "C": 1.0}, -1000.0, 0.0),      # blocked reverse
        ("EX_A", {"A": -1.0}, -1000.0, 1000.0),
    ])
    convert_to_irreversible(model)
    for rxn in model.reactions:
        if rxn.boundary:
            continue
        assert rxn.lower_bound >= 0, f"{rxn.id} still has lb < 0"


def test_returns_empty_list_when_nothing_to_split():
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0),
        ("EX_A", {"A": -1.0}, -1000.0, 1000.0),
    ])
    assert convert_to_irreversible(model) == []


def test_conversion_is_idempotent_after_first_pass():
    """Running convert_to_irreversible twice should not create
    `_REV_REV` reactions, because the first pass already clamped
    all non-exchange lb to 0."""
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),
    ])
    convert_to_irreversible(model)
    second = convert_to_irreversible(model)
    assert second == []
    assert "r1_REV_REV" not in {r.id for r in model.reactions}


def test_reverse_reaction_inherits_annotation_subsystem_and_notes():
    """MATLAB's convertToIrrev copies the per-reaction fields (eccodes,
    rxnMiriams, subSystems, rxnNotes, ...) onto the reverse reaction.
    Losing them here leaves e.g. GECKO's reverse reactions without an EC
    code, and therefore without a kcat."""
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),
    ])
    forward = model.reactions.get_by_id("r1")
    forward.annotation["ec-code"] = "1.1.1.1"
    forward.subsystem = "Glycolysis"
    forward.notes["origin"] = "curated"

    convert_to_irreversible(model)

    rev = model.reactions.get_by_id("r1_REV")
    assert rev.annotation["ec-code"] == "1.1.1.1"
    assert rev.subsystem == "Glycolysis"
    assert rev.notes["origin"] == "curated"


def test_reverse_reaction_annotation_is_a_copy_not_an_alias():
    """Editing the reverse reaction's annotations must not write through
    to the forward reaction."""
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),
    ])
    model.reactions.get_by_id("r1").annotation["ec-code"] = "1.1.1.1"

    convert_to_irreversible(model)

    model.reactions.get_by_id("r1_REV").annotation["ec-code"] = "2.2.2.2"
    assert model.reactions.get_by_id("r1").annotation["ec-code"] == "1.1.1.1"
