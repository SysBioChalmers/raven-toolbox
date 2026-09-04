"""Tests for flux-network navigation (analysis/walk.py, walkFluxes port)."""
import cobra
import pandas as pd
import pytest

from raven_toolbox.analysis import FluxWalker, walk_fluxes


@pytest.fixture
def model_and_fluxes():
    """R0 (A -> B, + a below-cutoff touch on H) is the reaction to walk from.

    A's other carriers: R_supply_A (produces A, flux 10) and R4 (A -> B,
    consumes A *and* produces B -- present in both groups, so it must keep
    the same display number rather than being renumbered under B).
    B's other carriers: R1 (consumes, flux 6), R2 (consumes, flux 4), R4
    (produces, flux 3, see above), and R3 (consumes, flux 1e-10 -- below the
    per-neighbour cutoff, must not appear at all). With max_per_met=2, B's
    list keeps only the two largest by |flux| (R1, R2); R4 is truncated out
    of B's own list even though it still holds the number it got under A.
    """
    m = cobra.Model("t")
    A, B, C, D, E, H = (cobra.Metabolite(x, compartment="c") for x in "ABCDEH")
    m.add_metabolites([A, B, C, D, E, H])

    r0 = cobra.Reaction("R0", lower_bound=-1000, upper_bound=1000)
    r0.add_metabolites({A: -1, B: 1, H: 1e-10})
    supply_a = cobra.Reaction("R_supply_A", lower_bound=-1000, upper_bound=1000)
    supply_a.add_metabolites({A: 1})
    r4 = cobra.Reaction("R4", lower_bound=-1000, upper_bound=1000)
    r4.add_metabolites({A: -1, B: 1})
    r1 = cobra.Reaction("R1", lower_bound=-1000, upper_bound=1000)
    r1.add_metabolites({B: -1, C: 1})
    r2 = cobra.Reaction("R2", lower_bound=-1000, upper_bound=1000)
    r2.add_metabolites({B: -1, D: 1})
    r3 = cobra.Reaction("R3", lower_bound=-1000, upper_bound=1000)
    r3.add_metabolites({B: -1, E: 1})
    m.add_reactions([r0, supply_a, r4, r1, r2, r3])

    fluxes = pd.Series(
        {"R0": 10.0, "R_supply_A": 10.0, "R4": 3.0, "R1": 6.0, "R2": 4.0, "R3": 1e-10},
    )
    return m, fluxes


def test_groups_and_roles(model_and_fluxes):
    model, fluxes = model_and_fluxes
    walker = FluxWalker(model, fluxes, "R0", max_per_met=2)
    groups = {g.metabolite: g for g in walker.groups}

    assert set(groups) == {"A", "B"}  # H filtered by the outer cutoff

    a = groups["A"]
    assert a.role == "consumed"
    assert [n.reaction for n in a.neighbors] == ["R_supply_A", "R4"]
    assert [n.role for n in a.neighbors] == ["produces", "consumes"]

    b = groups["B"]
    assert b.role == "produced"
    # R4 is a real candidate (flux 3) but truncated out by max_per_met=2,
    # even though it already holds a display number from group A.
    assert [n.reaction for n in b.neighbors] == ["R1", "R2"]


def test_neighbor_numbering_reused_across_groups(model_and_fluxes):
    model, fluxes = model_and_fluxes
    walker = FluxWalker(model, fluxes, "R0", max_per_met=2)
    assert walker.neighbor_ids == ["R_supply_A", "R4", "R1", "R2"]
    numbers = {n.reaction: n.number for g in walker.groups for n in g.neighbors}
    assert numbers == {"R_supply_A": 1, "R4": 2, "R1": 3, "R2": 4}


def test_below_cutoff_neighbor_excluded(model_and_fluxes):
    model, fluxes = model_and_fluxes
    walker = FluxWalker(model, fluxes, "R0")  # default max_per_met=8, no truncation
    b = next(g for g in walker.groups if g.metabolite == "B")
    assert "R3" not in [n.reaction for n in b.neighbors]  # flux 1e-10 <= default cutoff


def test_step_and_back(model_and_fluxes):
    model, fluxes = model_and_fluxes
    walker = FluxWalker(model, fluxes, "R0", max_per_met=2)
    walker.step(1)  # -> R_supply_A
    assert walker.current == "R_supply_A"
    assert walker.history == ["R0"]
    assert walker.back() is True
    assert walker.current == "R0"
    assert walker.history == []
    assert walker.back() is False  # nothing left to go back to


def test_step_out_of_range_raises(model_and_fluxes):
    model, fluxes = model_and_fluxes
    walker = FluxWalker(model, fluxes, "R0", max_per_met=2)
    with pytest.raises(ValueError, match="between 1 and 4"):
        walker.step(5)


def test_start_by_index(model_and_fluxes):
    model, fluxes = model_and_fluxes
    walker = FluxWalker(model, fluxes, 0)  # R0 is reactions[0]
    assert walker.current == "R0"


# --- walk_fluxes (REPL) ------------------------------------------------

def test_repl_step_back_quit(model_and_fluxes):
    model, fluxes = model_and_fluxes
    responses = iter(["1", "b", "q"])
    printed: list[str] = []
    walk_fluxes(
        model, fluxes, "R0", max_per_met=2,
        input_fn=lambda _: next(responses),
        print_fn=printed.append,
    )
    text = "\n".join(printed)
    assert "[R0]" in text
    assert "Navigator closed." in text


def test_repl_invalid_choice_reprompts(model_and_fluxes):
    model, fluxes = model_and_fluxes
    responses = iter(["99", "q"])
    printed: list[str] = []
    walk_fluxes(
        model, fluxes, "R0", max_per_met=2,
        input_fn=lambda _: next(responses),
        print_fn=printed.append,
    )
    assert any("Enter a number 1-4" in line for line in printed)
