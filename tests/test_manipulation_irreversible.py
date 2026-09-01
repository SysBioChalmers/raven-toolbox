"""Tests for convert_to_irreversible (RAVEN convertToIrrev.m).

Every reversible reaction is split, including exchange reactions, matching
MATLAB behavior.
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


def test_splits_exchange_reaction_when_reversible():
    """Exchange reactions (one metabolite) are split like any other
    reversible reaction in MATLAB, which does not special-case them."""
    model = _build_model_with_bounds([
        ("EX_A", {"A": -1.0}, -1000.0, 1000.0),
    ])
    added = convert_to_irreversible(model)
    assert added == ["EX_A_REV"]
    ex = model.reactions.get_by_id("EX_A")
    ex_rev = model.reactions.get_by_id("EX_A_REV")
    assert ex.bounds == (0.0, 1000.0)
    assert ex_rev.bounds == (0.0, 1000.0)
    assert {m.id: c for m, c in ex_rev.metabolites.items()} == {"A": 1.0}


def test_splits_multiple_mixed_reactions():
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),   # split
        ("r2", {"B": -2.0, "C": 3.0}, 0.0, 1000.0),      # forward only
        ("EX_A", {"A": -1.0}, -1000.0, 1000.0),          # exchange, split too
        ("r3", {"C": -1.0, "D": 1.0}, -200.0, 200.0),    # split
    ])

    added = convert_to_irreversible(model)
    assert added == ["EX_A_REV", "r1_REV", "r3_REV"]

    assert model.reactions.get_by_id("r1").bounds == (0.0, 1000.0)
    assert model.reactions.get_by_id("r1_REV").bounds == (0.0, 500.0)
    assert model.reactions.get_by_id("r2").bounds == (0.0, 1000.0)
    assert model.reactions.get_by_id("EX_A").bounds == (0.0, 1000.0)
    assert model.reactions.get_by_id("EX_A_REV").bounds == (0.0, 1000.0)
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


def test_no_reaction_has_negative_bound():
    """After conversion, no reaction -- exchange or not -- may carry
    negative flux."""
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),
        ("r2", {"B": -1.0, "C": 1.0}, -1000.0, 0.0),      # blocked reverse
        ("EX_A", {"A": -1.0}, -1000.0, 1000.0),
    ])
    convert_to_irreversible(model)
    for rxn in model.reactions:
        assert rxn.lower_bound >= 0, f"{rxn.id} still has lb < 0"


def test_returns_empty_list_when_nothing_to_split():
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0),
        ("EX_A", {"A": -1.0}, 0.0, 1000.0),
    ])
    assert convert_to_irreversible(model) == []


def test_conversion_is_idempotent_after_first_pass():
    """Running convert_to_irreversible twice should not create
    `_REV_REV` reactions, because the first pass already clamped
    all lb to 0."""
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


def test_negative_objective_coefficient_moves_to_reverse_reaction():
    """A negative objective coefficient on a split reaction is moved onto
    the new reverse reaction (sign-flipped positive) and zeroed on the
    forward reaction, matching MATLAB's convertToIrrev."""
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),
    ])
    model.objective = {model.reactions.get_by_id("r1"): -1.0}

    convert_to_irreversible(model)

    fwd = model.reactions.get_by_id("r1")
    rev = model.reactions.get_by_id("r1_REV")
    assert fwd.objective_coefficient == 0.0
    assert rev.objective_coefficient == 1.0


def test_positive_objective_coefficient_stays_on_forward_reaction():
    """A non-negative objective coefficient is left untouched on the
    forward reaction; the new reverse reaction gets no credit for it."""
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),
    ])
    model.objective = {model.reactions.get_by_id("r1"): 1.0}

    convert_to_irreversible(model)

    fwd = model.reactions.get_by_id("r1")
    rev = model.reactions.get_by_id("r1_REV")
    assert fwd.objective_coefficient == 1.0
    assert rev.objective_coefficient == 0.0


# --------------------------------------------------------------------------- #
# rxns parameter (MATLAB: convertToIrrev(model, 'rxns', rxns))
# --------------------------------------------------------------------------- #

def test_rxns_param_restricts_conversion_to_named_reactions():
    """Reversible reactions outside the given rxns list are left alone,
    even though they would otherwise be split -- this is how MATLAB's
    makeEcModel.m keeps exchange reactions unsplit (it passes
    nonExchRxns, the model's reactions minus its exchanges)."""
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),   # split
        ("EX_A", {"A": -1.0}, -1000.0, 1000.0),          # excluded
    ])

    added = convert_to_irreversible(model, rxns=["r1"])

    assert added == ["r1_REV"]
    assert "EX_A_REV" not in {r.id for r in model.reactions}
    assert model.reactions.get_by_id("EX_A").bounds == (-1000.0, 1000.0)


def test_rxns_param_none_converts_all_reversible_reactions():
    """The default (no rxns given) still splits every reversible
    reaction, including exchanges -- matching convertToIrrev's own
    model.rxns default."""
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, -500.0, 1000.0),
        ("EX_A", {"A": -1.0}, -1000.0, 1000.0),
    ])

    added = convert_to_irreversible(model)

    assert added == ["EX_A_REV", "r1_REV"]


def test_rxns_param_forward_only_reaction_in_list_is_not_split():
    """Being named in rxns is necessary but not sufficient -- the
    reaction still needs lb < 0 to be split."""
    model = _build_model_with_bounds([
        ("r1", {"A": -1.0, "B": 1.0}, 0.0, 1000.0),
    ])
    added = convert_to_irreversible(model, rxns=["r1"])
    assert added == []
