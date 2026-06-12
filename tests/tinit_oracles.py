"""Shared (ft)INIT test oracles, ported from RAVEN's ``tinitTests.m``.

These toy models have **defined reaction scores** and **known ftINIT outputs**, so
they serve as exact correctness oracles for the Phase 4d port (see
docs/ftinit_review_and_plan.md). Building them here once lets every sub-phase
(essential-reaction discovery, the MILP, linear merge, staging) check against the
same RAVEN-verified answers.

Reaction scores are injected through gene expression using :func:`expr_for_rxn_score`
(RAVEN's ``getExprForRxnScore``): each toy reaction ``Ri`` has at most one gene
``Gi``, so an expression of ``exp(score_i/5)`` reproduces the desired score exactly
(no-gene reactions get ``no_gene_score = -2`` regardless).
"""
from __future__ import annotations

import math

import cobra


def expr_for_rxn_score(scores, threshold: float = 1.0) -> dict:
    """RAVEN ``getExprForRxnScore``: gene expression giving a target single-gene score.

    Inverts ``score = 5·ln(level/threshold)`` → ``level = threshold·exp(score/5)``.
    Returns ``{Gi: level}`` for i = 1..len(scores) (gene name ``"G{i}"``), mirroring the
    1-reaction-1-gene layout of the toy models.
    """
    return {f"G{i + 1}": threshold * math.exp(s / 5) for i, s in enumerate(scores)}


def _build(model_id, mets, reactions, objective):
    """mets: {id: (name, compartment)}; reactions: {id: (stoich, lb, ub, gpr)}."""
    m = cobra.Model(model_id)
    met_objs = {
        mid: cobra.Metabolite(mid, name=name, compartment=comp)
        for mid, (name, comp) in mets.items()
    }
    m.add_metabolites(list(met_objs.values()))
    for rid, (stoich, lb, ub, gpr) in reactions.items():
        r = cobra.Reaction(rid, lower_bound=lb, upper_bound=ub)
        r.add_metabolites({met_objs[mid]: coeff for mid, coeff in stoich.items()})
        m.add_reactions([r])
        if gpr:
            r.gene_reaction_rule = gpr
    m.objective = objective
    return m


# --------------------------------------------------------------------------- #
# testModel — RAVEN getTstModel(): 8 mets, 10 rxns. a[s] -> ... -> e[s] export.
# --------------------------------------------------------------------------- #
_TEST_METS = {
    "as": ("a", "s"), "ac": ("a", "c"), "bc": ("b", "c"), "cc": ("c", "c"),
    "dc": ("d", "c"), "ec": ("e", "c"), "es": ("e", "s"), "fc": ("f", "c"),
}
_TEST_RXNS = {
    "R1": ({"as": 1}, 0, 1000, ""),                       # -> a[s]   (exchange, no GPR)
    "R2": ({"as": -1, "ac": 1}, -1000, 1000, ""),         # a[s] <=> a[c]  (transport, no GPR)
    "R3": ({"ac": -1, "bc": 1, "cc": 1}, -1000, 1000, "G3"),
    "R4": ({"ac": -1, "dc": 2}, -1000, 1000, "G4"),
    "R5": ({"bc": -1, "cc": -1, "ec": 1}, 0, 1000, "G5"),
    "R6": ({"dc": -2, "ec": 1}, 0, 1000, "G6"),
    "R7": ({"ec": -1, "es": 1}, 0, 1000, "G7"),           # transport, with GPR
    "R8": ({"es": -1}, 0, 1000, ""),                      # e[s] ->  (exchange, no GPR)
    "R9": ({"ac": -1, "fc": 1}, -1000, 1000, "G9"),
    "R10": ({"fc": -1, "ec": 1}, -1000, 1000, "G10"),
}
# RAVEN getTstModelRxnScores(), R1..R10.
TEST_MODEL_SCORES = [-2, -2, -1, 7, 0.5, 0.5, -1, -2, -3, 3.5]


def make_test_model() -> cobra.Model:
    return _build("testModel", _TEST_METS, _TEST_RXNS, "R8")


# Oracles (RAVEN tinitTests):
# T0001 ftINIT, no tasks, default '1+1':
TEST_MODEL_FTINIT_NO_TASKS = ["R1", "R4", "R6", "R8", "R9", "R10"]
# T0001 with R7,R10 spontaneous:
TEST_MODEL_FTINIT_SPONT_R7_R10 = ["R1", "R2", "R4", "R6", "R7", "R8"]
# T0002 with task "gen e[s] from a[s]": essential rxns (pre-merge ids) and output:
TEST_MODEL_TASK_ESSENTIAL_PREMERGE = ["R2", "R7"]
TEST_MODEL_TASK_ESSENTIAL_MERGED = ["R1", "R7"]
TEST_MODEL_FTINIT_WITH_TASK = ["R1", "R2", "R4", "R6", "R7", "R8", "R9", "R10"]
# T0004 mergeLinear(testModel): merges {R1,R2},{R3,R5},{R4,R6},{R7,R8},{R9,R10}
TEST_MODEL_GROUP_IDS = [1, 1, 2, 3, 2, 3, 4, 4, 5, 5]
TEST_MODEL_MERGED_REV = [0, 0, 0, 0, 1]
TEST_MODEL_MERGED_LB = [0, 0, 0, 0, -1000]
# groupRxnScores with R1,R2,R8 zeroed (toIgnore): -> per merged group
TEST_MODEL_GROUPED_SCORES = [0, -0.5, 7.5, -1, 0.5]


# The task: generate e[s] from a[s] (RAVEN getTstModelTasks()).
def make_test_task():
    """RAVEN getTstModelTasks(): make e[s] from a[s]."""
    from raven_toolbox.tasks import Task

    return Task(
        id="Gen e[s] from a[s]",
        description="Gen e[s] from a[s]",
        inputs=[("a[s]", 0.0, math.inf)],   # (token, LBin, UBin)
        outputs=[("e[s]", 1.0, 1.0)],       # (token, LBout, UBout)
    )


# --------------------------------------------------------------------------- #
# testModel4 — RAVEN getTstModel4(): partial linear merges + flips.
# --------------------------------------------------------------------------- #
_TEST4_METS = {
    "a": ("a", "s"), "b": ("b", "s"), "d": ("d", "s"), "e": ("e", "s"),
    "f": ("f", "s"), "g": ("g", "s"), "h": ("h", "s"),
}
_TEST4_RXNS = {
    "R1": ({"a": -1}, -1000, 1000, "G1"),                 # a[s] <=>
    "R2": ({"a": -1, "b": 1}, 0, 1000, "G2"),             # a[s] -> b[s]
    "R3": ({"a": -1, "b": 1}, -1000, 1000, "G3"),         # a[s] <=> b[s]
    "R4": ({"b": -1}, 0, 1000, "G4"),                     # b[s] ->
    "R5": ({"a": -5, "d": 5}, -1000, 1000, "G5"),         # 5 a[s] <=> 5 d[s]
    "R6": ({"e": -1, "d": 1}, -1000, 1000, "G6"),         # e[s] <=> d[s]
    "R7": ({"f": -1, "g": -1, "e": 1}, -1000, 1000, "G7"),  # f[s]+g[s] <=> e[s]
    "R8": ({"b": -1, "f": 1}, -1000, 1000, "G8"),         # b[s] <=> f[s]
    "R9": ({"h": -1, "g": 1}, -1000, 1000, "G9"),         # h[s] <=> g[s]
    "R10": ({"h": -1}, 0, 1000, "G10"),                   # h[s] ->
    "R11": ({"e": -1, "g": 1}, 0, 1000, "G11"),           # e[s] -> g[s]
}
TEST_MODEL4_SCORES = [-1, -1, 2, -1, 0.5, -2, 1, 1.3, -0.5, -0.4, 8]
# T0004 mergeLinear(testModel4): merges {R5,R6},{R7,R8},{R9,R10}; rest unmerged.
TEST_MODEL4_GROUP_IDS = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 0]
TEST_MODEL4_MERGED_REV = [1, 0, 1, 0, 1, 1, 0, 0]
TEST_MODEL4_REVERSED_RXNS = ["R6", "R9"]  # flipped direction when made irreversible


def make_test_model4() -> cobra.Model:
    return _build("testModel4", _TEST4_METS, _TEST4_RXNS, "R4")


# --------------------------------------------------------------------------- #
# testModel5 — RAVEN getTstModel5(): testModel + an unmerged parallel path R11-R14.
# --------------------------------------------------------------------------- #
def make_test_model5() -> cobra.Model:
    m = make_test_model()
    m.id = "testModel5"
    m.add_metabolites([cobra.Metabolite("gc", name="g", compartment="c")])
    gc = m.metabolites.get_by_id("gc")
    ac = m.metabolites.get_by_id("ac")
    ec = m.metabolites.get_by_id("ec")
    extra = {
        "R11": ({ac: -1, gc: 1}, -1000, 1000, "G11"),
        "R12": ({ac: -1, gc: 1}, -1000, 1000, "G12"),
        "R13": ({gc: -1, ec: 1}, -1000, 1000, "G13"),
        "R14": ({gc: -1, ec: 1}, -1000, 1000, "G14"),
    }
    for rid, (stoich, lb, ub, gpr) in extra.items():
        r = cobra.Reaction(rid, lower_bound=lb, upper_bound=ub)
        r.add_metabolites(stoich)
        m.add_reactions([r])
        r.gene_reaction_rule = gpr
    return m


TEST_MODEL5_SCORES = [*TEST_MODEL_SCORES, -1, -1.5, -1, -1.5]
