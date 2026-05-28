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


def test_predict_empty_relocate_set_is_no_op():
    """An empty relocate set short-circuits to an empty proposal."""
    m = _toy_two_compartment_model()
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0, 0.0], "m": [0.0, 1.0]}, index=pd.Index(["g1", "g2"], name="gene_id")))
    res = predict_localization(m, scores, reactions_to_relocate=[], apply=False)
    assert isinstance(res, LocalizationProposal)
    assert res.moved.empty


def test_predict_places_single_reaction():
    """Pass r2 in the relocate set; it goes to 'm' per scores."""
    m = _toy_two_compartment_model()
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0, 0.0], "m": [0.0, 1.0]}, index=pd.Index(["g1", "g2"], name="gene_id")))
    res = predict_localization(m, scores, ["r2"], default_compartment="c", apply=False,
                                transport_cost=0.1)
    assert isinstance(res, LocalizationProposal)
    assert set(res.moved["rxn_id"]) == {"r2"}
    assert res.moved.iloc[0]["to_compartment"] == "m"


def test_predict_apply_creates_compartment_metabolites_and_transports():
    """apply=True should mutate the (copy) model: r2 in m, and B/C transports added."""
    m = _toy_two_compartment_model()
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0, 0.0], "m": [0.0, 1.0]}, index=pd.Index(["g1", "g2"], name="gene_id")))
    res = predict_localization(m, scores, ["r2"], default_compartment="c", apply=True,
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


def test_predict_unplaced_reaction_when_no_scored_gene():
    """A relocate-set reaction whose genes are all absent from scores is reported, not crashed."""
    m = _toy_two_compartment_model()
    # Only g1 has scores; g2 (r2's gene) is absent → r2 is unplaceable.
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0], "m": [0.0]}, index=pd.Index(["g1"], name="gene_id")))
    res = predict_localization(m, scores, ["r2"], apply=False)
    assert isinstance(res, LocalizationProposal)
    assert "r2" in res.unplaced_reactions
    assert res.moved.empty   # nothing actually placed


def test_predict_boundary_reactions_always_pinned():
    """Boundary reactions in the relocate set are silently filtered out."""
    m = _toy_two_compartment_model()
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0, 0.0], "m": [0.0, 1.0]}, index=pd.Index(["g1", "g2"], name="gene_id")))
    res = predict_localization(m, scores, ["EX_A", "EX_C", "r2"], apply=False,
                                transport_cost=0.1)
    # Only r2 should appear in the proposal — boundaries dropped.
    assert set(res.moved["rxn_id"]) == {"r2"}


def test_predict_default_compartment_validated():
    m = _toy_two_compartment_model()
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0], "m": [0.0]}, index=pd.Index(["g2"], name="gene_id")))
    with pytest.raises(ValueError, match="default_compartment"):
        predict_localization(m, scores, ["r2"], default_compartment="x", apply=False)


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


# ----------------------------------------- multi-compartment scoring (NEW)

def test_predict_multi_compartment_when_secondary_score_beats_penalty():
    """Dual-localised gene: secondary compartment score 0.8 > penalty 0.3 → gene lands in
    both compartments. Two reactions sharing one gene each placed in their best
    compartment without contradicting the gene assignment."""
    m = cobra.Model("dual")
    A_c = cobra.Metabolite("A_c", compartment="c")
    B_c = cobra.Metabolite("B_c", compartment="c")
    m.add_metabolites([A_c, B_c])
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({A_c: -1, B_c: 1})
    r1.gene_reaction_rule = "g_dual"
    r2 = cobra.Reaction("r2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({A_c: -1, B_c: 1})
    r2.gene_reaction_rule = "g_dual"
    m.add_reactions([r1, r2])
    # g_dual scores: c=1.0 (primary), m=0.8 (secondary). Penalty 0.3 — secondary worth it.
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0], "m": [0.8]}, index=pd.Index(["g_dual"], name="gene_id")))
    res = predict_localization(m, scores, ["r2"], default_compartment="c", apply=False,
                                transport_cost=0.0, multi_compartment_penalty=0.3)
    # The gene should land in BOTH c and m (primary free + 0.8 - 0.3 > 0 for secondary).
    assert set(res.gene_compartments["g_dual"]) == {"c", "m"}


def test_predict_mono_when_secondary_score_below_penalty():
    """Same as above but penalty 0.9 > secondary score 0.8 → gene mono-localises (c only)."""
    m = cobra.Model("mono")
    A_c = cobra.Metabolite("A_c", compartment="c")
    B_c = cobra.Metabolite("B_c", compartment="c")
    m.add_metabolites([A_c, B_c])
    r2 = cobra.Reaction("r2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({A_c: -1, B_c: 1})
    r2.gene_reaction_rule = "g_dual"
    m.add_reactions([r2])
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0], "m": [0.8]}, index=pd.Index(["g_dual"], name="gene_id")))
    res = predict_localization(m, scores, ["r2"], default_compartment="c", apply=False,
                                transport_cost=0.0, multi_compartment_penalty=0.9)
    # Penalty exceeds secondary score → only the primary compartment.
    assert res.gene_compartments["g_dual"] == ["c"]


def test_predict_high_penalty_forces_mono_localisation():
    """Very high penalty effectively bans extra compartments."""
    m = _toy_two_compartment_model()
    scores = LocalizationScores(pd.DataFrame(
        {"c": [0.4, 0.3], "m": [0.5, 0.6]},
        index=pd.Index(["g1", "g2"], name="gene_id")))
    res = predict_localization(m, scores, ["r1", "r2"], default_compartment="c", apply=False,
                                transport_cost=0.0, multi_compartment_penalty=1000.0)
    # With penalty so high, every gene gets exactly one compartment.
    for g, comps in res.gene_compartments.items():
        assert len(comps) == 1, f"{g} landed in {comps} with prohibitive penalty"
