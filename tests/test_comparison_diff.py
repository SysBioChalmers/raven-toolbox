"""Tests for raven_python.comparison.diff — strict 2-model equality."""
from __future__ import annotations

import cobra
import pytest

from raven_python.comparison import DiffReport, diff_models
from raven_python.comparison.diff import _normalise_gpr


def _mini_model(model_id: str = "m") -> cobra.Model:
    """Tiny but realistic model: one reaction, one extracellular met."""
    m = cobra.Model(model_id)
    a = cobra.Metabolite("a_c", name="A", compartment="c", charge=0, formula="C1")
    b = cobra.Metabolite("b_e", name="B", compartment="e", charge=0, formula="C2")
    m.add_metabolites([a, b])
    r = cobra.Reaction("r1", lower_bound=-1000, upper_bound=1000)
    r.add_metabolites({a: -1, b: 1})
    r.gene_reaction_rule = "g1 AND g2"
    r.annotation["sbo"] = "SBO:0000176"
    m.add_reactions([r])
    return m


def test_model_equal_to_itself():
    m = _mini_model()
    report = diff_models(m, m)
    assert isinstance(report, DiffReport)
    assert report.equal


def test_independent_copies_are_equal():
    a = _mini_model()
    b = _mini_model()
    assert diff_models(a, b).equal


def test_dropped_reaction_is_detected():
    a = _mini_model()
    b = _mini_model()
    b.remove_reactions([b.reactions[0]])
    report = diff_models(a, b)
    assert not report.equal
    assert any("reactions only in A" in d for d in report.differences)


def test_bound_diff_is_detected():
    a = _mini_model()
    b = _mini_model()
    b.reactions[0].lower_bound = -42
    report = diff_models(a, b)
    assert not report.equal
    assert any("bounds" in d for d in report.differences)


def test_stoich_within_tolerance_is_equal():
    a = _mini_model()
    b = _mini_model()
    met = next(iter(b.reactions[0].metabolites))
    b.reactions[0].add_metabolites({met: 1e-12}, combine=True)
    assert diff_models(a, b, stoichiometry_tol=1e-9).equal


def test_stoich_outside_tolerance_is_not_equal():
    a = _mini_model()
    b = _mini_model()
    met = next(iter(b.reactions[0].metabolites))
    b.reactions[0].add_metabolites({met: 5.0}, combine=False)
    report = diff_models(a, b, stoichiometry_tol=1e-9)
    assert not report.equal
    assert any("coef[" in d for d in report.differences)


def test_annotation_diff_detected():
    a = _mini_model()
    b = _mini_model()
    b.reactions[0].annotation["sbo"] = "SBO:0009999"
    report = diff_models(a, b)
    assert not report.equal
    assert any("annotation['sbo']" in d for d in report.differences)


def test_ignore_annotations_suppresses_diff():
    a = _mini_model()
    b = _mini_model()
    b.reactions[0].annotation["sbo"] = "SBO:0009999"
    assert diff_models(a, b, ignore_annotations={"sbo"}).equal


def test_extra_annotations_picked_up():
    a = _mini_model()
    b = _mini_model()
    a.reactions[0].annotation["custom-key"] = "v1"
    b.reactions[0].annotation["custom-key"] = "v2"
    # Default key set ignores custom-key → equal.
    assert diff_models(a, b).equal
    # Pull it in → diff.
    assert not diff_models(a, b, extra_annotations={"custom-key"}).equal


@pytest.mark.parametrize(
    "ga, gb",
    [
        ("A and B", "a AND b"),
        ("A  and   B", "A and B"),
        ("(A or B) and C", "(a OR b) AND c"),
    ],
)
def test_gpr_normalisation(ga, gb):
    assert _normalise_gpr(ga) == _normalise_gpr(gb)


def test_max_per_category_truncates():
    """The per-category cap kicks in on per-reaction (and per-met) diffs,
    not on the id-set summary line. Build a scenario with many bound
    diffs across shared reactions."""
    a = cobra.Model("a")
    b = cobra.Model("b")
    mets_a = [cobra.Metabolite(f"m{i}_c", compartment="c") for i in range(60)]
    mets_b = [cobra.Metabolite(f"m{i}_c", compartment="c") for i in range(60)]
    a.add_metabolites(mets_a)
    b.add_metabolites(mets_b)
    for i in range(60):
        ra = cobra.Reaction(f"r{i}", lower_bound=0, upper_bound=1000)
        ra.add_metabolites({mets_a[i]: -1})
        a.add_reactions([ra])
        rb = cobra.Reaction(f"r{i}", lower_bound=-7, upper_bound=1000)  # differs
        rb.add_metabolites({mets_b[i]: -1})
        b.add_reactions([rb])
    report = diff_models(a, b, max_per_category=5)
    assert not report.equal
    assert any("truncated at 5" in d for d in report.differences)


def test_report_str_lists_differences():
    a = _mini_model()
    b = _mini_model()
    b.remove_reactions([b.reactions[0]])
    text = str(diff_models(a, b))
    assert "Models differ" in text
    assert "reactions only in A" in text


def test_report_str_when_equal():
    m = _mini_model()
    assert str(diff_models(m, m)) == "Models are semantically equal."
