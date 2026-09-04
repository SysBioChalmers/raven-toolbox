"""Tests for generate_new_ids."""
import cobra
import pytest

from raven_toolbox.utils import generate_new_ids


def _model_with_reactions(ids):
    m = cobra.Model("t")
    for rid in ids:
        m.add_reactions([cobra.Reaction(rid)])
    return m


def test_generate_new_ids_starts_at_one_when_none_exist():
    m = _model_with_reactions([])
    assert generate_new_ids(m, "reactions", "r_") == ["r_0001"]


def test_generate_new_ids_continues_after_existing_max():
    m = _model_with_reactions(["r_0001", "r_0003", "unrelated"])
    assert generate_new_ids(m, "reactions", "r_", quantity=2) == ["r_0004", "r_0005"]


def test_generate_new_ids_adopts_existing_width_over_requested_num_length():
    # Existing ids are zero-padded to 6 digits; the default/requested
    # num_length=4 is overridden by that, matching generateNewIds.m.
    m = _model_with_reactions(["r_000042"])
    assert generate_new_ids(m, "reactions", "r_", num_length=4) == ["r_000043"]


def test_generate_new_ids_metabolites():
    m = cobra.Model("t")
    m.add_metabolites([cobra.Metabolite("s_0001")])
    assert generate_new_ids(m, "metabolites", "s_") == ["s_0002"]


def test_generate_new_ids_rejects_unknown_entity_type():
    m = _model_with_reactions([])
    with pytest.raises(ValueError, match="entity_type"):
        generate_new_ids(m, "genes", "g_")


def test_generate_new_ids_non_numeric_suffix_falls_back_to_zero():
    # A stripped id that isn't a plain integer (e.g. hand-edited) can't be
    # parsed as a starting point -- matching str2double producing NaN,
    # coerced to 0 in generateNewIds.m -- so numbering restarts at 1.
    m = _model_with_reactions(["r_abcd"])
    assert generate_new_ids(m, "reactions", "r_") == ["r_0001"]


def test_generate_new_ids_sorts_stripped_ids_as_strings_not_integers():
    # Matches generateNewIds.m exactly, quirk included: the existing max is
    # found by a string sort of the stripped ids, not a numeric one. Mixed
    # widths can make that disagree with a numeric sort -- "9" sorts after
    # "10" as an integer but before it as a string, so with both present the
    # last (and therefore adopted) width comes from whichever string is
    # alphabetically greatest, here "9" (width 1) rather than "10" (width 2).
    m = _model_with_reactions(["r_9", "r_10"])
    assert generate_new_ids(m, "reactions", "r_") == ["r_10"]
