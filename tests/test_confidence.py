"""Tests for per-reaction confidence: the data model, notes round-trip (YAML + SBML), the localization
scorer, and the invariant that a confidence-annotated model still loads and solves in plain cobra."""
import cobra
import pandas as pd

from raven_toolbox.confidence import (
    ConfidenceEntry,
    ReactionConfidence,
    confidence_report,
    get_confidence,
    mark_curated,
    read_confidence,
    score_localization_confidence,
    set_confidence,
)
from raven_toolbox.localization import AssignmentProposal
from raven_toolbox.localization.scores import LocalizationScores


def _model():
    m = cobra.Model("t")
    a = cobra.Metabolite("A_c", name="A", compartment="c", formula="C6H12O6")
    b = cobra.Metabolite("B_c", name="B", compartment="c", formula="C6H12O6")
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
    assert n == 2
    loc1 = get_confidence(m.reactions.r1).facets["localization"]
    assert loc1.score == 0.9 and loc1.level == "strong" and "deeploc" in loc1.basis
    assert "fba-certified" in loc1.basis                    # proposal was certified
    loc2 = get_confidence(m.reactions.r2).facets["localization"]
    assert loc2.score == 0.0 and loc2.level == "none" and loc2.basis == "connectivity"


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


def test_confidence_report_is_lowest_first():
    m, proposal, scores = _scored_proposal()
    score_localization_confidence(m, proposal, scores)
    set_confidence(m.reactions.r1, "equation", ConfidenceEntry(0.4, basis="imbalanced"))
    rep = confidence_report(m)
    assert list(rep.columns) == ["reaction", "equation", "localization", "overall"]
    assert rep["overall"].is_monotonic_increasing          # weakest first
    assert rep.iloc[0]["reaction"] == "r2"                 # r2 (0.0) is the least confident
    assert set(read_confidence(m)) == {"r1", "r2"}
