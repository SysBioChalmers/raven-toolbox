"""The staged ftINIT pipeline (prep_init_model + get_init_steps + ftinit).

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

from raven_toolbox.init import (
    classify_reactions,
    ftinit,
    get_init_steps,
    prep_init_model,
    score_reactions_from_genes,
)
from raven_toolbox.init.score import gene_scores_from_expression


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


def test_pipeline_with_gene_scores_and_tasks_wires_up():
    """ftinit accepts gene_scores (gene pruning) + tasks (gap-fill) without breaking T0002.

    The toy's GPRs are single-gene (nothing to prune) and the task is feasible in the
    extracted model (nothing to gap-fill), so the reaction set is unchanged — this
    confirms the integration wiring (the pruning/gap-fill logic is unit-tested
    separately in test_init_genes / test_init_taskfill).
    """
    model = make_test_model()
    gene_scores = gene_scores_from_expression(expr_for_rxn_score(TEST_MODEL_SCORES), 1.0)
    prep = prep_init_model(model, [make_test_task()], ext_comp="s")
    out = ftinit(prep, _scores(model), gene_scores=gene_scores)
    assert {r.id for r in out.reactions} == set(TEST_MODEL_FTINIT_WITH_TASK)


def test_orient_forward_reverses_a_reversible_reaction():
    """_orient_forward(rxn, -1) flips stoichiometry and makes it irreversible forward."""
    import cobra

    from raven_toolbox.init.prep import _orient_forward

    m = cobra.Model("o")
    a, b = (cobra.Metabolite(x, compartment="s") for x in "ab")
    m.add_metabolites([a, b])
    r = cobra.Reaction("R", lower_bound=-800, upper_bound=1000)
    r.add_metabolites({a: -1, b: 2})  # a <=> 2 b
    m.add_reactions([r])

    _orient_forward(r, -1)  # forced reverse → becomes forward
    assert r.bounds == (0, 800)  # [-800,1000] → flip [-1000,800] → lb→0
    assert {mt.id: c for mt, c in r.metabolites.items()} == {"a": 1, "b": -2}  # 2 b => a

    fwd = cobra.Reaction("F", lower_bound=-500, upper_bound=900)
    fwd.add_metabolites({a: -1})
    m.add_reactions([fwd])
    _orient_forward(fwd, 1)  # forced forward → just made irreversible
    assert fwd.bounds == (0, 900)


def test_essential_merged_away_is_skipped():
    """An essential reaction whose merge group collapses away imposes no constraint.

    REV sits between two exchanges, so it merges with them into a trivial source→sink
    that is removed; its group has no survivor. prep_init_model must skip it, not crash.
    """
    import cobra

    from raven_toolbox.tasks import Task

    m = cobra.Model("collapse")
    a, b = (cobra.Metabolite(x, name=x, compartment="s") for x in "ab")
    m.add_metabolites([a, b])
    r = cobra.Reaction("REV", lower_bound=-1000, upper_bound=1000)
    r.add_metabolites({a: -1, b: 1})
    r.gene_reaction_rule = "g1"
    exchanges = []
    for met in (a, b):
        ex = cobra.Reaction(f"EX_{met.id}", lower_bound=-1000, upper_bound=1000)
        ex.add_metabolites({met: -1})
        exchanges.append(ex)
    m.add_reactions([r, *exchanges])
    m.objective = "REV"
    task = Task(id="mk_a", inputs=[("b[s]", 0.0, 1000.0)], outputs=[("a[s]", 1.0, 1.0)])

    prep = prep_init_model(m, [task], ext_comp="s")  # must not raise
    assert "REV" not in prep.essential_rxns  # merged into a collapsed group


# --------------------------------------------------------------------------- #
# prove_abs_gap / resolve_ties (deterministic extraction).
# These are opt-in; the default path stays exact-RAVEN. On the toy oracle they must
# reproduce T0001 (no regression) and be run-to-run identical.
# --------------------------------------------------------------------------- #
def test_ftinit_prove_abs_gap_matches_oracle():
    """prove_abs_gap: one near-proven-optimal solve per step, still gives T0001."""
    model = make_test_model()
    prep = prep_init_model(model, ext_comp="s")
    out = ftinit(prep, _scores(model), prove_abs_gap=0.05)
    assert {r.id for r in out.reactions} == set(TEST_MODEL_FTINIT_NO_TASKS)


def test_ftinit_resolve_ties_matches_oracle_and_is_stable():
    """resolve_ties (+ prove_abs_gap): preserves T0001, identical across repeated runs."""
    model = make_test_model()
    prep = prep_init_model(model, ext_comp="s")
    out1 = ftinit(prep, _scores(model), prove_abs_gap=0.05, resolve_ties=True)
    out2 = ftinit(prep, _scores(model), prove_abs_gap=0.05, resolve_ties=True)
    assert {r.id for r in out1.reactions} == set(TEST_MODEL_FTINIT_NO_TASKS)
    assert {r.id for r in out1.reactions} == {r.id for r in out2.reactions}


# --------------------------------------------------------------------------- #
# reference_reactions (stability under similar input) — axis-B follow-up.
# resolve_ties/prove_abs_gap give determinism (same input -> same output); this
# targets stability (similar input -> similar output) by preferring the tied
# solution closest to a reference build instead of the sparsest/lowest-id one.
# --------------------------------------------------------------------------- #
def _degenerate_prep():
    """R1/R2 are mutually redundant producers of an essential intermediate.

    Forced essential directly on the built PrepData (RAVEN's own task-essential-reaction
    detection correctly finds NEITHER individually essential here -- removing just one
    still leaves the task feasible through the other, so this isolates the tie-break
    question under test rather than fighting that, correct, detection).
    """
    import cobra

    m = cobra.Model("degen")
    a, mm, p = (cobra.Metabolite(x, name=x, compartment="s") for x in ("a", "m", "p"))
    m.add_metabolites([a, mm, p])
    R1 = cobra.Reaction("R1", lower_bound=0, upper_bound=1000)
    R1.add_metabolites({a: -1, mm: 1})
    R2 = cobra.Reaction("R2", lower_bound=0, upper_bound=1000)
    R2.add_metabolites({a: -1, mm: 1})
    E = cobra.Reaction("E", lower_bound=0, upper_bound=1000)
    E.add_metabolites({mm: -1, p: 1})
    EXa = cobra.Reaction("EX_a", lower_bound=-1000, upper_bound=1000)
    EXa.add_metabolites({a: -1})
    EXp = cobra.Reaction("EX_p", lower_bound=-1000, upper_bound=1000)
    EXp.add_metabolites({p: -1})
    m.add_reactions([R1, R2, E, EXa, EXp])
    m.objective = "E"
    prep = prep_init_model(m, ext_comp="s")
    prep.essential_rxns = {"E"}
    return prep


def test_reference_reactions_redirects_the_tie_break():
    """reference_reactions flips a degenerate choice, ahead of the id-rank default.

    Exercises the full staged ftinit() pipeline (prep_init_model + merge id bookkeeping),
    not just the single-step run_ftinit MILP: this is the layer that has to translate a
    reference given in original ids into each step's merged-reaction id space.
    """
    prep = _degenerate_prep()
    scores = {"R1": -1.0, "R2": -1.0}

    baseline = ftinit(prep, scores, resolve_ties=True, fill_gaps=False)
    baseline_ids = {r.id for r in baseline.reactions}
    assert "R1" in baseline_ids and "R2" not in baseline_ids  # default: lower id wins

    anchored = ftinit(prep, scores, resolve_ties=True, fill_gaps=False,
                      reference_reactions={"R2", "E", "EX_a", "EX_p"})
    anchored_ids = {r.id for r in anchored.reactions}
    assert "R2" in anchored_ids and "R1" not in anchored_ids  # reference wins instead


def test_reference_reactions_self_anchoring_is_idempotent():
    """Anchoring a build to its own output must reproduce it exactly."""
    model = make_test_model()
    prep = prep_init_model(model, ext_comp="s")
    baseline = ftinit(prep, _scores(model), resolve_ties=True, fill_gaps=False)
    baseline_ids = {r.id for r in baseline.reactions}

    anchored = ftinit(prep, _scores(model), resolve_ties=True, fill_gaps=False,
                      reference_reactions=baseline_ids)
    assert {r.id for r in anchored.reactions} == baseline_ids


def test_reference_reactions_requires_resolve_ties():
    """reference_reactions without resolve_ties is a usage error, not a silent no-op."""
    import pytest

    model = make_test_model()
    prep = prep_init_model(model, ext_comp="s")
    with pytest.raises(ValueError, match="resolve_ties=True"):
        ftinit(prep, _scores(model), reference_reactions={"R1"}, fill_gaps=False)


def test_reference_reactions_translates_through_merge_groups():
    """The reference->merged-id translation matches via ANY member of a merged group.

    testModel merges R1+R2 (survivor R1) and R3+R5 (survivor R3). Naming only the
    non-survivor member of each group must still count as a match.
    """
    from raven_toolbox.init.ftinit import _translate_reference

    model = make_test_model()
    prep = prep_init_model(model, ext_comp="s")
    assert prep.group_of["R1"] == prep.group_of["R2"] != 0
    assert prep.group_of["R3"] == prep.group_of["R5"] != 0

    matched = _translate_reference(prep, {"R2", "R5", "unknown_reaction_id"})
    assert matched == {"R1", "R3"}
