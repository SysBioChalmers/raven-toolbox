"""Phase 4d.3b: the staged ftINIT pipeline (prep_init_model + get_init_steps + ftinit).

Oracles: RAVEN tinitTests T0001/T0002 on testModel with the default '1+1' schedule.
"""

from tinit_oracles import (
    TEST_MODEL_FTINIT_NO_TASKS,
    TEST_MODEL_FTINIT_SPONT_R7_R10,
    TEST_MODEL_FTINIT_WITH_TASK,
    TEST_MODEL_SCORES,
    TEST_MODEL_TASK_ESSENTIAL_MERGED,
    expr_for_rxn_score,
    make_test_model,
    make_test_task,
)

from ravengem.init import (
    classify_reactions,
    ftinit,
    get_init_steps,
    prep_init_model,
    score_reactions_from_genes,
)
from ravengem.init.score import gene_scores_from_expression


def _scores(model):
    return score_reactions_from_genes(
        model, gene_scores_from_expression(expr_for_rxn_score(TEST_MODEL_SCORES), 1.0)
    )


# --------------------------------------------------------------------------- #
# classify_reactions (the toIgnore masks) — tinitTests T0001 mask oracle.
# --------------------------------------------------------------------------- #
def test_classify_exchange_and_transport():
    masks = classify_reactions(make_test_model(), ext_comp="s")
    assert masks.exchange == {"R1", "R8"}        # boundary reactions
    assert masks.import_rxns == {"R2"}           # a[s] <=> a[c], no GPR, into ext comp
    assert masks.no_gpr == {"R1", "R2", "R8"}
    assert "R7" not in masks.import_rxns         # R7 has a GPR -> not a transport category


def test_classify_spontaneous():
    masks = classify_reactions(make_test_model(), ext_comp="s", spontaneous=["R7", "R10"])
    assert masks.exchange | masks.spontaneous == {"R1", "R7", "R8", "R10"}


def test_get_init_steps_default():
    steps = get_init_steps("1+1")
    assert len(steps) == 2
    assert steps[0].how_to_use_prev == "ignore"
    assert steps[0].ignore_mask == (1, 1, 1, 1, 1, 1, 1, 0)
    assert steps[1].how_to_use_prev == "essential"
    assert steps[1].ignore_mask == (1, 0, 0, 0, 1, 0, 0, 0)
    assert len(get_init_steps("full")) == 1


# --------------------------------------------------------------------------- #
# Full '1+1' pipeline — T0001 (no tasks) and T0002 (with task).
# --------------------------------------------------------------------------- #
def test_ftinit_no_tasks_matches_oracle():
    """T0001: testModel, no tasks, '1+1' → {R1,R4,R6,R8,R9,R10}."""
    model = make_test_model()
    prep = prep_init_model(model, ext_comp="s")
    out = ftinit(prep, _scores(model))
    assert {r.id for r in out.reactions} == set(TEST_MODEL_FTINIT_NO_TASKS)


def test_ftinit_with_spontaneous_matches_oracle():
    """T0001 variant: R7,R10 spontaneous → the path through R2/R7, {R1,R2,R4,R6,R7,R8}."""
    model = make_test_model()
    prep = prep_init_model(model, ext_comp="s", spontaneous=["R7", "R10"])
    out = ftinit(prep, _scores(model))
    assert {r.id for r in out.reactions} == set(TEST_MODEL_FTINIT_SPONT_R7_R10)


def test_ftinit_with_task_matches_oracle():
    """T0002: task 'make e[s] from a[s]' makes R2,R7 essential → {R1,R2,R4,R6,R7,R8,R9,R10}."""
    model = make_test_model()
    prep = prep_init_model(model, [make_test_task()], ext_comp="s")
    # Essentials map to merged ids {R1, R7} (RAVEN T0002).
    assert prep.essential_rxns == set(TEST_MODEL_TASK_ESSENTIAL_MERGED)
    out = ftinit(prep, _scores(model))
    assert {r.id for r in out.reactions} == set(TEST_MODEL_FTINIT_WITH_TASK)


def test_full_series_runs():
    """The single-step 'full' series also produces a feasible subnetwork."""
    model = make_test_model()
    prep = prep_init_model(model, ext_comp="s")
    out = ftinit(prep, _scores(model), series="full")
    assert len(out.reactions) >= 1
