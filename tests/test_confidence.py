"""Tests for per-reaction confidence: the data model, notes round-trip (YAML + SBML), the localization
scorer, and the invariant that a confidence-annotated model still loads and solves in plain cobra."""
import cobra
import pandas as pd
import pytest

from raven_toolbox.confidence import (
    ConfidenceEntry,
    ReactionConfidence,
    confidence_report,
    equation_exempt,
    facet_summary,
    gene_association_exempt,
    get_confidence,
    mark_curated,
    read_confidence,
    score_equation_confidence,
    score_gene_association_confidence,
    score_localization_confidence,
    set_confidence,
)
from raven_toolbox.localization import AssignmentProposal
from raven_toolbox.localization.scores import LocalizationScores


def _model():
    m = cobra.Model("t")
    a = cobra.Metabolite("A_c", name="A", compartment="c", formula="C6H12O6", charge=0)
    b = cobra.Metabolite("B_c", name="B", compartment="c", formula="C6H12O6", charge=0)
    m.add_metabolites([a, b])
    ex = cobra.Reaction("EX_A", lower_bound=-10, upper_bound=1000)
    ex.add_metabolites({a: -1})
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({a: -1, b: 1})
    r1.gene_reaction_rule = "g1"
    bio = cobra.Reaction("bio", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({b: -1})
    m.add_reactions([ex, r1, bio])
    m.objective = "bio"
    return m


# --------------------------------------------------------------------------- data model

def test_entry_clamps_score_and_serialises():
    e = ConfidenceEntry(score=1.5, basis="deeploc", level="strong", note="hi")
    assert e.score == 1.0                                   # clamped to [0, 1]
    d = e.to_dict()
    assert d == {"score": 1.0, "level": "strong", "basis": "deeploc", "note": "hi"}  # None fields dropped
    assert ConfidenceEntry.from_dict(d).score == 1.0
    assert ConfidenceEntry(score=-0.2).score == 0.0


def test_reaction_confidence_overall_is_weakest_facet():
    rc = ReactionConfidence(facets={"localization": ConfidenceEntry(0.9),
                                    "equation": ConfidenceEntry(0.3)})
    assert rc.overall == 0.3
    assert ReactionConfidence().overall is None             # empty -> None


# --------------------------------------------------------------------------- storage

def test_set_get_clear_confidence():
    m = _model()
    r = m.reactions.r1
    assert get_confidence(r).facets == {}                   # nothing yet
    set_confidence(r, "localization", ConfidenceEntry(0.8, basis="deeploc", level="strong"))
    set_confidence(r, "equation", ConfidenceEntry(1.0, basis="mass-balanced"))
    rc = get_confidence(r)
    assert set(rc.facets) == {"localization", "equation"}
    assert rc.facets["localization"].score == 0.8
    assert rc.overall == 0.8
    from raven_toolbox.confidence import clear_confidence
    clear_confidence(r, "equation")
    assert set(get_confidence(r).facets) == {"localization"}
    clear_confidence(r)
    assert get_confidence(r).facets == {} and "raven_confidence" not in r.notes


def test_yaml_round_trip_and_cobra_still_solves(tmp_path):
    m = _model()
    set_confidence(m.reactions.r1, "localization",
                   ConfidenceEntry(0.82, basis="deeploc", level="strong", note="x"))
    p = tmp_path / "m.yml"
    cobra.io.save_yaml_model(m, str(p))
    m2 = cobra.io.load_yaml_model(str(p))
    rc = get_confidence(m2.reactions.r1)
    assert rc.facets["localization"].score == 0.82
    assert rc.facets["localization"].basis == "deeploc"
    assert m2.slim_optimize() > 0                            # confidence annotation does not break FBA


def test_sbml_round_trip_survives_html_escaping(tmp_path):
    m = _model()
    set_confidence(m.reactions.r1, "localization",
                   ConfidenceEntry(0.7, basis="deeploc", level="moderate", note='has "quotes"'))
    p = tmp_path / "m.xml"
    cobra.io.write_sbml_model(m, str(p))
    m2 = cobra.io.read_sbml_model(str(p))
    # SBML stores notes as an HTML-escaped string; get_confidence unescapes + json.loads it
    rc = get_confidence(m2.reactions.r1)
    assert rc.facets["localization"].score == 0.7
    assert m2.slim_optimize() > 0


def test_malformed_notes_value_is_ignored():
    m = _model()
    m.reactions.r1.notes["raven_confidence"] = "not json {["
    assert get_confidence(m.reactions.r1).facets == {}      # unreadable -> empty, no crash


# --------------------------------------------------------------------------- localization scorer

def _scored_proposal():
    """r1 placed in m with strong gene evidence; r2 placed in c by connectivity (no scored gene)."""
    m = _model()
    r2 = cobra.Reaction("r2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({m.metabolites.B_c: -1})
    r2.gene_reaction_rule = "g2"                            # g2 is NOT in the score table
    m.add_reactions([r2])
    scores = LocalizationScores(pd.DataFrame({"c": [0.1], "m": [0.9]}, index=["g1"]))
    proposal = AssignmentProposal(placements={"r1": ["m"], "r2": ["c"]},
                                  unplaced_reactions=["r2"], certified=True, status="certified")
    return m, proposal, scores


def test_score_localization_confidence_reflects_evidence():
    m, proposal, scores = _scored_proposal()
    n = score_localization_confidence(m, proposal, scores)
    loc1 = get_confidence(m.reactions.r1).facets["localization"]
    assert loc1.score == 0.9 and loc1.level == "strong" and "deeploc" in loc1.basis
    assert "fba-certified" in loc1.basis                    # proposal was certified
    # r2's gene is absent from the score table: no measurement was possible, so the scorer ABSTAINS
    # rather than writing 0.0 -- which, under overall=min, would veto r2 on a missing input.
    assert n == 1
    assert "localization" not in get_confidence(m.reactions.r2).facets


def test_localization_abstention_clears_a_stale_score():
    m, proposal, scores = _scored_proposal()
    set_confidence(m.reactions.r2, "localization", ConfidenceEntry(0.4, basis="deeploc"))
    score_localization_confidence(m, proposal, scores)
    assert "localization" not in get_confidence(m.reactions.r2).facets


def test_localization_zero_is_a_measurement_not_ignorance():
    """A gene that IS scored, but scores 0.0 at the assigned compartment, is evidence *against* the
    placement -- that earns a real 0.0, unlike the unmeasurable case above."""
    m = _model()
    scores = LocalizationScores(pd.DataFrame({"c": [0.0], "m": [1.0]}, index=["g1"]))
    proposal = AssignmentProposal(placements={"r1": ["c"]}, certified=False)
    assert score_localization_confidence(m, proposal, scores) == 1
    loc = get_confidence(m.reactions.r1).facets["localization"]
    assert loc.score == 0.0 and loc.level == "none" and loc.basis == "deeploc"


def test_mark_curated_and_scorer_respects_it():
    m, proposal, scores = _scored_proposal()
    mark_curated(m.reactions.r1, source="curator:eduardk")
    curated = get_confidence(m.reactions.r1).facets["localization"]
    assert curated.score == 1.0 and curated.level == "curated"
    # a later automated pass leaves the curated call alone...
    score_localization_confidence(m, proposal, scores)
    assert get_confidence(m.reactions.r1).facets["localization"].level == "curated"
    # ...unless explicitly overwritten
    score_localization_confidence(m, proposal, scores, overwrite_curated=True)
    assert get_confidence(m.reactions.r1).facets["localization"].level == "strong"


def test_mark_curated_can_pin_any_facet():
    # mark_curated used to hardcode "localization", which made every other scorer's `curated` guard
    # unreachable through the public API.
    m = _model()
    mark_curated(m.reactions.r1, facet="equation", source="curator:eduardk")
    assert get_confidence(m.reactions.r1).facets["equation"].level == "curated"
    score_equation_confidence(m)                            # r1 balances, so it would score 1.0/"strong"
    assert get_confidence(m.reactions.r1).facets["equation"].level == "curated"
    score_equation_confidence(m, overwrite_curated=True)
    assert get_confidence(m.reactions.r1).facets["equation"].score == 1.0


def test_confidence_report_is_lowest_first():
    m, proposal, scores = _scored_proposal()
    score_localization_confidence(m, proposal, scores)
    set_confidence(m.reactions.r1, "equation", ConfidenceEntry(0.4, basis="imbalanced"))
    set_confidence(m.reactions.r2, "equation", ConfidenceEntry(0.0, basis="mass-imbalanced"))
    rep = confidence_report(m)
    assert list(rep.columns) == ["reaction", "equation", "localization", "overall"]
    assert rep["overall"].is_monotonic_increasing          # weakest first
    assert rep.iloc[0]["reaction"] == "r2"                 # r2 (0.0) is the least confident
    assert set(read_confidence(m)) == {"r1", "r2"}


# --------------------------------------------------------------- exemptions (which facets apply)

def _sbo_model():
    """One reaction of each kind the exemptions must tell apart."""
    m = cobra.Model("sbo")
    atp = cobra.Metabolite("atp_c", formula="C10H12N5O13P3", charge=-4, compartment="c")
    adp = cobra.Metabolite("adp_c", formula="C10H12N5O10P2", charge=-3, compartment="c")
    pi = cobra.Metabolite("pi_c", formula="HO4P", charge=-2, compartment="c")
    h2o = cobra.Metabolite("h2o_c", formula="H2O", charge=0, compartment="c")
    h = cobra.Metabolite("h_c", formula="H", charge=1, compartment="c")
    m.add_metabolites([atp, adp, pi, h2o, h])

    ngam = cobra.Reaction("NGAM")                     # ATP + H2O -> ADP + Pi + H : real chemistry
    ngam.add_metabolites({atp: -1, h2o: -1, adp: 1, pi: 1, h: 1})
    ngam.annotation["sbo"] = "SBO:0000630"            # ATP maintenance
    ngam.name = "non-growth associated maintenance reaction"

    ex = cobra.Reaction("EX_h2o", lower_bound=-10)    # boundary: one metabolite
    ex.add_metabolites({h2o: -1})

    bio = cobra.Reaction("BIOMASS")
    bio.add_metabolites({atp: -1, adp: 1})
    bio.annotation["sbo"] = "SBO:0000629"

    slime = cobra.Reaction("SLIME")
    slime.add_metabolites({atp: -1, adp: 1})
    slime.annotation["sbo"] = "SBO:0000395"

    spont = cobra.Reaction("SPONT")
    spont.add_metabolites({h2o: -1, h: 1})            # mass-imbalanced on purpose
    spont.annotation["sbo"] = "SBO:0000672"

    h2o_m = cobra.Metabolite("h2o_m", formula="H2O", charge=0, compartment="m")
    m.add_metabolites([h2o_m])
    transport = cobra.Reaction("T_h2o")               # transport, no gene: a real curation gap
    transport.add_metabolites({h2o: -1, h2o_m: 1})
    transport.annotation["sbo"] = "SBO:0000655"

    m.add_reactions([ngam, ex, bio, slime, spont, transport])
    return m


def test_equation_exempts_boundary_biomass_pseudo_but_never_atp_maintenance():
    m = _sbo_model()
    assert equation_exempt(m.reactions.EX_h2o) == "boundary"
    assert equation_exempt(m.reactions.BIOMASS) == "biomass"
    assert equation_exempt(m.reactions.SLIME) == "pseudoreaction"
    # ATP maintenance is real chemistry that must balance -- a name regex on "non-growth associated
    # maintenance reaction" would wrongly silence the check.
    assert equation_exempt(m.reactions.NGAM) is None
    assert equation_exempt(m.reactions.T_h2o) is None


def test_gene_association_exempts_maintenance_and_spontaneous_but_not_transport():
    m = _sbo_model()
    assert gene_association_exempt(m.reactions.NGAM) == "maintenance"   # real chemistry, no catalyst
    assert gene_association_exempt(m.reactions.SPONT) == "spontaneous"
    assert gene_association_exempt(m.reactions.BIOMASS) == "biomass"
    # a transport reaction with no transporter gene is a genuine gap, not a convention
    assert gene_association_exempt(m.reactions.T_h2o) is None


def test_exempt_reactions_are_not_scored_and_stale_scores_are_cleared():
    m = _sbo_model()
    set_confidence(m.reactions.EX_h2o, "equation", ConfidenceEntry(0.0, basis="mass-imbalanced"))
    n = score_equation_confidence(m)
    # scored: NGAM, SPONT, T_h2o. exempt: EX_h2o (boundary), BIOMASS, SLIME.
    assert n == 3
    assert "equation" not in get_confidence(m.reactions.EX_h2o).facets   # stale score removed
    assert "equation" not in get_confidence(m.reactions.BIOMASS).facets
    assert get_confidence(m.reactions.NGAM).facets["equation"].score == 1.0   # balances -> scored 1.0


def test_abstention_never_deletes_a_curators_call():
    """`overwrite_curated=True` licenses *recomputing* over a curated entry, not deleting one. A scorer
    that cannot measure a reaction has learned nothing that invalidates a human's assertion about it --
    and silently dropping the annotation from the model file would be unrecoverable."""
    m = _sbo_model()
    mark_curated(m.reactions.EX_h2o, facet="equation", note="checked by hand")   # exempt: boundary
    score_equation_confidence(m)
    assert get_confidence(m.reactions.EX_h2o).facets["equation"].level == "curated"
    score_equation_confidence(m, overwrite_curated=True)
    assert get_confidence(m.reactions.EX_h2o).facets["equation"].level == "curated"


def test_scorers_are_idempotent():
    m = _sbo_model()
    first = (score_equation_confidence(m), score_gene_association_confidence(m))
    snapshot = {r.id: dict(r.notes) for r in m.reactions}
    second = (score_equation_confidence(m), score_gene_association_confidence(m))
    assert first == second
    assert {r.id: dict(r.notes) for r in m.reactions} == snapshot   # no drift, no accumulation


def test_no_sbo_model_warns_that_pseudo_reactions_cannot_be_told_from_defects():
    m = _model()  # no reaction carries an SBO term
    with pytest.warns(UserWarning, match="SBO"):
        score_equation_confidence(m)


# --------------------------------------------------------------------------- equation facet

def _eq_model():
    m = cobra.Model("eq")
    a = cobra.Metabolite("a_c", formula="C6H12O6", charge=0, compartment="c")
    b = cobra.Metabolite("b_c", formula="C6H12O6", charge=0, compartment="c")
    c = cobra.Metabolite("c_c", formula="C3H6O3", charge=0, compartment="c")
    d = cobra.Metabolite("d_c", formula="C6H12O6", charge=-1, compartment="c")
    m.add_metabolites([a, b, c, d])

    def rx(rid, stoich):
        r = cobra.Reaction(rid, lower_bound=0, upper_bound=1000)
        r.add_metabolites(stoich)
        m.add_reactions([r])
        return r

    rx("balanced", {a: -1, b: 1})
    rx("mass_imbal", {a: -1, c: 1})
    rx("charge_imbal", {a: -1, d: 1})
    return m


def test_equation_bands_and_defect_ordering():
    m = _eq_model()
    assert score_equation_confidence(m) == 3
    f = {r.id: get_confidence(r).facets["equation"] for r in m.reactions}
    assert (f["balanced"].score, f["balanced"].basis) == (1.0, "balanced")
    assert (f["mass_imbal"].score, f["mass_imbal"].basis) == (0.0, "mass-imbalanced")
    assert (f["charge_imbal"].score, f["charge_imbal"].basis) == (0.1, "charge-imbalanced")
    # a proven mass error is worse than a proven charge error (often just a protonation convention)
    assert f["mass_imbal"].score < f["charge_imbal"].score
    assert "C-3" in f["mass_imbal"].note and "charge" in f["charge_imbal"].note


def test_missing_formula_is_ignorance_and_outranks_every_proven_defect():
    """The invariant: max(defect) < min(ignorance). A reaction we cannot check must never sort above
    one we checked and found broken."""
    m = _eq_model()
    m.metabolites.b_c.formula = None
    score_equation_confidence(m)
    f = {r.id: get_confidence(r).facets["equation"] for r in m.reactions}
    assert (f["balanced"].score, f["balanced"].basis) == (0.3, "formula-missing")
    assert "b_c" in f["balanced"].note
    assert f["mass_imbal"].score < f["balanced"].score      # 0.0 defect  <  0.3 ignorance
    assert f["charge_imbal"].score < f["balanced"].score    # 0.1 defect  <  0.3 ignorance


def test_unparseable_polymer_formula_is_unknown_not_a_crash():
    m = _eq_model()
    m.metabolites.b_c.formula = "(C5H8)n"                   # glycogen-style; cobra raises on this
    score_equation_confidence(m)
    entry = get_confidence(m.reactions.balanced).facets["equation"]
    assert (entry.score, entry.basis) == (0.3, "formula-unparseable")


def test_generic_r_group_residual_is_not_a_proven_imbalance():
    """A residual landing in an R group is uninterpretable, not a defect -- and cobra's element table
    contains U and F, so a hand-rolled element set would read FULLR2's letters as uranium/fluorine."""
    m = _eq_model()
    m.metabolites.b_c.formula = "C6H12O6R"                  # R group only on the product side
    score_equation_confidence(m)
    entry = get_confidence(m.reactions.balanced).facets["equation"]
    assert (entry.score, entry.basis) == (0.3, "formula-generic")
    assert "R" in entry.note


def test_r_group_that_cancels_is_noted_but_still_scores_balanced():
    m = _eq_model()
    m.metabolites.a_c.formula = "C6H12O6R"
    m.metabolites.b_c.formula = "C6H12O6R"                  # R cancels: the reaction does balance
    score_equation_confidence(m)
    entry = get_confidence(m.reactions.balanced).facets["equation"]
    assert entry.score == 1.0 and entry.basis == "balanced"
    assert "R" in entry.note                                # recorded, but not scored against


def test_charge_residual_is_never_fabricated_from_a_partial_sum():
    """cobra accumulates check_mass_balance()['charge'] over only the non-None metabolites, inventing a
    residual from half the reaction. With any charge unset the verdict must be `charge-unknown`."""
    m = _eq_model()
    m.metabolites.a_c.charge = None                         # b_c keeps charge 0
    score_equation_confidence(m)
    entry = get_confidence(m.reactions.balanced).facets["equation"]
    assert (entry.score, entry.basis) == (0.6, "charge-unknown")
    assert entry.level == "moderate"                        # mass proven; only charge is unverifiable


# --------------------------------------------------------------------- gene_association facet

def test_gene_association_bands_track_thiele_palsson():
    m = _eq_model()
    m.reactions.balanced.gene_reaction_rule = "g1 or g2"
    m.reactions.balanced.annotation["pubmed"] = "12345"      # genetic + literature  ~ TP 3
    m.reactions.mass_imbal.gene_reaction_rule = "g3 and g4"  # genetic only          ~ TP 2
    #  charge_imbal keeps no GPR                             # no evidence           ~ TP 0/1
    assert score_gene_association_confidence(m) == 3
    f = {r.id: get_confidence(r).facets["gene_association"] for r in m.reactions}
    assert (f["balanced"].score, f["balanced"].basis) == (0.9, "gpr+literature")
    assert (f["mass_imbal"].score, f["mass_imbal"].basis) == (0.6, "gpr")
    assert (f["charge_imbal"].score, f["charge_imbal"].basis) == (0.2, "no-gpr")
    # 1.0 is reserved for a curator's call, so an inferred score never ties one
    assert max(e.score for e in f.values()) < 1.0
    # GPR shape is recorded, never scored from
    assert f["balanced"].note == "2 isozymes"
    assert f["mass_imbal"].note == "complex of 2"


@pytest.mark.parametrize("recorded", ["3", 3, 3.0, "3.0", " 3 "])
def test_gene_association_flags_a_model_contradicting_its_own_recorded_confidence(recorded):
    # A YAML round-trip readily turns the note into 3.0 / "3.0"; parsing must not silently drop it.
    m = _eq_model()
    m.reactions.balanced.notes["Confidence Level"] = recorded  # claims genetic+biochemical evidence...
    score_gene_association_confidence(m)                       # ...but carries no GPR
    entry = get_confidence(m.reactions.balanced).facets["gene_association"]
    assert entry.score == 0.2 and "recorded Confidence Level 3 but no GPR" in entry.note


@pytest.mark.parametrize("recorded", ["high", "", None, ["3"]])
def test_unparseable_recorded_confidence_is_ignored(recorded):
    m = _eq_model()
    m.reactions.balanced.notes["Confidence Level"] = recorded
    score_gene_association_confidence(m)
    assert get_confidence(m.reactions.balanced).facets["gene_association"].note is None


def test_defects_outrank_ignorance_across_both_facets():
    """The single invariant that makes `overall` a usable review queue: `overall == 0.0` means evidence
    contradicts the model, never that evidence is missing."""
    from raven_toolbox.confidence import (
        _EQ_CHARGE_IMBALANCED,
        _EQ_CHARGE_UNKNOWN,
        _EQ_FORMULA_UNKNOWN,
        _EQ_MASS_IMBALANCED,
        _GA_GPR,
        _GA_GPR_LITERATURE,
        _GA_NO_GPR,
    )
    defects = {_EQ_MASS_IMBALANCED, _EQ_CHARGE_IMBALANCED}
    ignorance = {_EQ_FORMULA_UNKNOWN, _EQ_CHARGE_UNKNOWN, _GA_NO_GPR, _GA_GPR, _GA_GPR_LITERATURE}
    assert max(defects) < min(ignorance)
    assert max(ignorance) < 1.0                              # 1.0 stays reserved


def test_facet_summary_separates_scored_from_exempt():
    m = _sbo_model()
    score_equation_confidence(m)
    score_gene_association_confidence(m)
    s = facet_summary(m)
    assert list(s.columns) == ["facet", "basis", "score", "level", "n"]
    # Only NGAM, SPONT and T_h2o carry an equation facet; the exempt 3 carry none. So "verified
    # balanced" reads 2 of 3 *checked* (NGAM, T_h2o) -- not 5 of 6, which is what writing 1.0 for the
    # exempt reactions would have made it look like.
    eq = s[s["facet"] == "equation"]
    assert eq["n"].sum() == 3
    assert int(eq[eq["basis"] == "balanced"]["n"].iloc[0]) == 2
    assert int(eq[eq["basis"] == "mass-imbalanced"]["n"].iloc[0]) == 1   # SPONT: H2O -> H
