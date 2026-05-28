"""Tests for ravengem.localization — predictor loaders + the MILP + apply (Phase 7)."""
from __future__ import annotations

from textwrap import dedent

import cobra
import pandas as pd
import pytest

from ravengem.localization import (
    LocalizationProposal,
    LocalizationResult,
    LocalizationScores,
    apply_localization,
    load_deeploc,
    load_wolfpsort,
    predict_localization,
)

# --------------------------------------------------------------------- loaders

def test_load_wolfpsort_basic(tmp_path):
    p = tmp_path / "wolf.txt"
    p.write_text(dedent("""\
        # header comment
        Gene1 cyto 13, nucl 7, mito 4
        Gene2: treating 9 X's as Glycines
        Gene3 mito 20, cyto 2
    """))
    s = load_wolfpsort(p)
    assert "Gene1" in s.genes
    assert "Gene2" not in s.genes      # the 'treating' line is skipped
    assert "Gene3" in s.genes
    # row-normalised to max=1:
    assert s.df.loc["Gene1", "cyto"] == pytest.approx(1.0)   # 13/13
    assert s.df.loc["Gene1", "nucl"] == pytest.approx(7 / 13)
    assert s.df.loc["Gene3", "mito"] == pytest.approx(1.0)
    assert s.df.loc["Gene3", "cyto"] == pytest.approx(0.1)


def test_load_deeploc_csv(tmp_path):
    p = tmp_path / "deeploc.csv"
    p.write_text(dedent("""\
        Protein_ID,Localizations,Signals,Cytoplasm,Nucleus,Mitochondrion
        G1,Cytoplasm,,0.8,0.1,0.05
        G2,Mitochondrion,SP,0.05,0.15,0.9
    """))
    s = load_deeploc(p)
    assert set(s.compartments) == {"Cytoplasm", "Nucleus", "Mitochondrion"}
    # row-max → 1.0
    assert s.df.loc["G1", "Cytoplasm"] == pytest.approx(1.0)
    assert s.df.loc["G2", "Mitochondrion"] == pytest.approx(1.0)


def test_localization_scores_with_compartments_rename():
    df = pd.DataFrame({"cyto": [1.0], "mito": [0.2]}, index=pd.Index(["g1"], name="gene_id"))
    s = LocalizationScores(df).with_compartments({"cyto": "c", "mito": "m"})
    assert list(s.compartments) == ["c", "m"]


# ----------------------------------------------------------------- predict (toy)

def _toy_two_compartment_model() -> cobra.Model:
    """Single-compartment draft (everything in 'c'):

    A_c -(r1)-> B_c -(r2)-> C_c        (r2 should move to 'm' per scores below)
    Boundary EX_A imports A; EX_C drains C.
    """
    m = cobra.Model("toy")
    A, B, C = (cobra.Metabolite(x + "_c", name=x, compartment="c") for x in "ABC")
    m.add_metabolites([A, B, C])

    def rxn(rid, lb, ub, mets, gpr=None):
        r = cobra.Reaction(rid, lower_bound=lb, upper_bound=ub)
        r.add_metabolites(mets)
        if gpr:
            r.gene_reaction_rule = gpr
        return r
    m.add_reactions([rxn("EX_A", -1000, 0, {A: -1}),
                     rxn("EX_C", 0, 1000, {C: -1}),
                     rxn("r1", 0, 1000, {A: -1, B: 1}, "g1"),
                     rxn("r2", 0, 1000, {B: -1, C: 1}, "g2")])
    return m


def test_predict_no_uncertain_reactions_returns_empty(tmp_path):
    """With no reaction flagged 'uncertain' and no explicit list, nothing moves."""
    m = _toy_two_compartment_model()
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0, 0.0], "m": [0.0, 1.0]}, index=pd.Index(["g1", "g2"], name="gene_id")))
    res = predict_localization(m, scores, apply=False)
    assert isinstance(res, LocalizationProposal)
    assert res.moved.empty


def test_predict_relocates_only_flagged_reaction():
    """Flag r2 uncertain; only r2 is placed (it goes to 'm' per scores)."""
    m = _toy_two_compartment_model()
    m.reactions.r2.notes["localization"] = "uncertain"
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0, 0.0], "m": [0.0, 1.0]}, index=pd.Index(["g1", "g2"], name="gene_id")))
    res = predict_localization(m, scores, default_compartment="c", apply=False,
                                transport_cost=0.1)
    assert isinstance(res, LocalizationProposal)
    assert set(res.moved["rxn_id"]) == {"r2"}
    assert res.moved.iloc[0]["to_compartment"] == "m"


def test_predict_apply_creates_compartment_metabolites_and_transports():
    """apply=True should mutate the (copy) model: r2 in m, and B/C transports added."""
    m = _toy_two_compartment_model()
    m.reactions.r2.notes["localization"] = "uncertain"
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0, 0.0], "m": [0.0, 1.0]}, index=pd.Index(["g1", "g2"], name="gene_id")))
    res = predict_localization(m, scores, default_compartment="c", apply=True,
                                transport_cost=0.05)
    assert isinstance(res, LocalizationResult)
    out = res.model
    r2 = out.reactions.r2
    assert {mt.compartment for mt in r2.metabolites} == {"m"}   # both substrates now in m
    # B_m and C_m metabolite copies must exist:
    assert "B_m" in out.metabolites and "C_m" in out.metabolites
    # Transports tr_B_m and tr_C_m must be added (default c ↔ m):
    transport_ids = {t.id for t in res.added_transports}
    assert "tr_B_m" in transport_ids
    assert "tr_C_m" in transport_ids
    # Original model untouched (we copied).
    assert m.reactions.r2.metabolites != r2.metabolites
    assert "B_m" not in m.metabolites


def test_predict_explicit_relocate_set_overrides_auto():
    """Passing reactions_to_relocate explicitly ignores the 'uncertain' flag mechanism."""
    m = _toy_two_compartment_model()
    # No flag set; explicit list still kicks in.
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0, 0.0], "m": [0.0, 1.0]}, index=pd.Index(["g1", "g2"], name="gene_id")))
    # transport_cost=0.1 so placing r2 in 'm' (gene score 1.0, 2 transports = -0.2) beats
    # placing in 'c' (gene score 0.0, no transports). Default 0.5 would tie at 0.0.
    res = predict_localization(m, scores, reactions_to_relocate=["r2"], apply=False,
                                transport_cost=0.1)
    assert isinstance(res, LocalizationProposal)
    assert list(res.moved["rxn_id"]) == ["r2"]


def test_predict_unplaced_reaction_when_no_scored_gene():
    """A relocate-set reaction whose genes are all absent from scores is reported, not crashed."""
    m = _toy_two_compartment_model()
    m.reactions.r2.notes["localization"] = "uncertain"
    # Only g1 has scores; g2 (r2's gene) is absent → r2 is unplaceable.
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0], "m": [0.0]}, index=pd.Index(["g1"], name="gene_id")))
    res = predict_localization(m, scores, apply=False)
    assert isinstance(res, LocalizationProposal)
    assert "r2" in res.unplaced_reactions
    assert res.moved.empty   # nothing actually placed


def test_predict_boundary_reactions_always_pinned():
    """Even if a boundary reaction is in reactions_to_relocate, it's filtered out."""
    m = _toy_two_compartment_model()
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0, 0.0], "m": [0.0, 1.0]}, index=pd.Index(["g1", "g2"], name="gene_id")))
    res = predict_localization(m, scores,
                                reactions_to_relocate=["EX_A", "EX_C", "r2"], apply=False)
    # Only r2 should appear in the proposal — boundaries dropped.
    assert set(res.moved["rxn_id"]) <= {"r2"}


def test_predict_default_compartment_validated():
    m = _toy_two_compartment_model()
    m.reactions.r2.notes["localization"] = "uncertain"
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0], "m": [0.0]}, index=pd.Index(["g2"], name="gene_id")))
    with pytest.raises(ValueError, match="default_compartment"):
        predict_localization(m, scores, default_compartment="x", apply=False)


def test_apply_localization_idempotent_on_empty_proposal():
    """An empty proposal (no moves, no transports) shouldn't change the model."""
    m = _toy_two_compartment_model()
    empty = LocalizationProposal(
        moved=pd.DataFrame(columns=["rxn_id", "from_compartment", "to_compartment"]),
        added_transports=pd.DataFrame(columns=["met_id", "compartment"]),
        gene_compartments={})
    out, added = apply_localization(m, empty)
    assert len(out.reactions) == len(m.reactions)
    assert added == []
