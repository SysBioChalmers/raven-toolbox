"""Phase 4d.2: linear reaction merging (merge_linear + group_rxn_scores).

Oracles: RAVEN tinitTests T0004. testModel merges {R1,R2},{R3,R5},{R4,R6},{R7,R8},
{R9,R10}; testModel4 merges {R5,R6},{R7,R8},{R9,R10} with two reactions flipped.
"""
import pytest
from tinit_oracles import (
    TEST_MODEL4_GROUP_IDS,
    TEST_MODEL4_MERGED_REV,
    TEST_MODEL4_REVERSED_RXNS,
    TEST_MODEL_GROUP_IDS,
    TEST_MODEL_GROUPED_SCORES,
    TEST_MODEL_MERGED_LB,
    TEST_MODEL_MERGED_REV,
    TEST_MODEL_SCORES,
    make_test_model,
    make_test_model4,
)

from ravengem.init import group_rxn_scores, merge_linear


def test_test_model_group_ids():
    _, orig_ids, group_ids, _ = merge_linear(make_test_model())
    assert orig_ids == [f"R{i}" for i in range(1, 11)]
    assert group_ids == TEST_MODEL_GROUP_IDS  # [1,1,2,3,2,3,4,4,5,5]


def test_test_model_reduced_shape():
    reduced, _, _, _ = merge_linear(make_test_model())
    # Five merged reactions, survivors keep the producer's id, original order.
    assert [r.id for r in reduced.reactions] == ["R1", "R3", "R4", "R7", "R9"]
    assert [int(r.lower_bound < 0) for r in reduced.reactions] == TEST_MODEL_MERGED_REV
    assert [r.lower_bound for r in reduced.reactions] == TEST_MODEL_MERGED_LB


def test_test_model_grouped_scores():
    reduced, orig_ids, group_ids, _ = merge_linear(make_test_model())
    scores = dict(zip(orig_ids, TEST_MODEL_SCORES, strict=True))
    grouped = group_rxn_scores(reduced, scores, orig_ids, group_ids,
                               to_zero={"R1", "R2", "R8"})
    got = [grouped[r.id] for r in reduced.reactions]
    assert got == pytest.approx(TEST_MODEL_GROUPED_SCORES)  # [0,-0.5,7.5,-1,0.5]


def test_test_model4_group_ids_and_flips():
    reduced, orig_ids, group_ids, reversed_rxns = merge_linear(make_test_model4())
    assert group_ids == TEST_MODEL4_GROUP_IDS  # [0,0,0,0,1,1,2,2,3,3,0]
    assert [int(r.lower_bound < 0) for r in reduced.reactions] == TEST_MODEL4_MERGED_REV
    flipped = {oid for oid, rev in zip(orig_ids, reversed_rxns, strict=True) if rev}
    assert flipped == set(TEST_MODEL4_REVERSED_RXNS)  # {R6, R9}


def test_merge_preserves_feasible_space():
    """The reduced model admits flux through the merged export path, like the original.

    The reduced model carries no objective (merging drops genes and objective; ftINIT
    sets its own from scores), so we set one on the surviving export reaction. R8
    (e[s]=>) was merged into R7 (grp4), so R7 is the reduced export.
    """
    original = make_test_model()
    assert original.slim_optimize() > 1e-9  # exports e via R8
    reduced, _, _, _ = merge_linear(original)
    reduced.objective = "R7"
    assert reduced.slim_optimize() > 1e-9


def test_no_merge_blocks_merging():
    """A reaction in no_merge keeps its own group (id 0) and is not contracted."""
    _, orig_ids, group_ids, _ = merge_linear(make_test_model(), no_merge=["R2"])
    g = dict(zip(orig_ids, group_ids, strict=True))
    assert g["R2"] == 0  # R2 never merged
    # R1 was only mergeable with R2, so it stays unmerged too.
    assert g["R1"] == 0


def test_group_scores_zero_handling():
    """Genuine-zero score → 0.01; a group cancelling to zero with nonzero members → 0.01."""
    reduced, orig_ids, group_ids, _ = merge_linear(make_test_model())
    # Give group {R3,R5} scores that cancel: R3=+1, R5=-1 -> sum 0 but members nonzero.
    scores = dict.fromkeys(orig_ids, 0.0)
    scores["R3"], scores["R5"] = 1.0, -1.0
    grouped = group_rxn_scores(reduced, scores, orig_ids, group_ids)
    assert grouped["R3"] == pytest.approx(0.01)        # cancelled group rescued
    assert grouped["R4"] == pytest.approx(0.02)         # {R4,R6} both genuine-0 → 0.01+0.01
