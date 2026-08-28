"""``review_identity`` reports the near misses without letting them in.

The filters are set strict on purpose: a wrongly transferred reaction is hard to
find later and harder to remove, while a missing one can be gap-filled. That
argument justifies rejecting borderline hits -- it does not justify throwing away
the evidence for them, which is what these tests pin down.
"""
from __future__ import annotations

import cobra
import pandas as pd
import pytest

from raven_toolbox.reconstruction.homology import get_model_from_homology
from raven_toolbox.reconstruction.homology.hits import HIT_COLUMNS


def _template() -> cobra.Model:
    model = cobra.Model("tmpl")
    mets = {m: cobra.Metabolite(m, compartment="c") for m in ("a", "b", "c")}
    model.add_metabolites(list(mets.values()))
    for rid, stoich, gene in (
        ("STRONG", {"a": -1, "b": 1}, "gStrong"),
        ("WEAK", {"b": -1, "c": 1}, "gWeak"),
    ):
        rxn = cobra.Reaction(rid, lower_bound=0.0, upper_bound=1000.0)
        rxn.add_metabolites({mets[m]: v for m, v in stoich.items()})
        model.add_reactions([rxn])
        rxn.gene_reaction_rule = gene
    return model


def _hits(weak_identity: float) -> pd.DataFrame:
    """One comfortable ortholog and one whose identity is the variable."""
    rows = []
    for template_gene, target_gene, identity in (
        ("gStrong", "tStrong", 90.0),
        ("gWeak", "tWeak", weak_identity),
    ):
        rows.append(("tgt", "tmpl", target_gene, template_gene, 0.0, identity, 500, 400.0, 100.0))
        rows.append(("tmpl", "tgt", template_gene, target_gene, 0.0, identity, 500, 400.0, 100.0))
    return pd.DataFrame(rows, columns=HIT_COLUMNS)


def test_no_candidates_reported_unless_asked():
    result = get_model_from_homology([_template()], _hits(30.0), "tgt")

    assert result.candidates is None
    assert {r.id for r in result.model.reactions} == {"STRONG"}


def test_borderline_reaction_is_reported_but_not_transferred():
    result = get_model_from_homology(
        [_template()], _hits(30.0), "tgt", review_identity=25
    )

    assert {r.id for r in result.model.reactions} == {"STRONG"}, (
        "a candidate must never reach the model -- that is the whole point"
    )
    assert list(result.candidates.reaction) == ["WEAK"]

    row = result.candidates.iloc[0]
    assert row.template_gene == "gWeak"
    assert row.target_gene == "tWeak"
    assert row.identity == pytest.approx(30.0)
    assert row.n_support == 1


def test_reactions_that_pass_are_not_also_listed_as_candidates():
    result = get_model_from_homology(
        [_template()], _hits(80.0), "tgt", review_identity=25
    )

    assert {r.id for r in result.model.reactions} == {"STRONG", "WEAK"}
    assert result.candidates.empty, "nothing was rejected, so nothing is a near miss"


def test_candidates_are_ordered_strongest_first():
    """A curator reads from the top, so the best-supported must be there."""
    template = _template()
    extra = cobra.Reaction("WEAKER", lower_bound=0.0, upper_bound=1000.0)
    extra.add_metabolites({template.metabolites.get_by_id("c"): -1})
    template.add_reactions([extra])
    extra.gene_reaction_rule = "gWeaker"

    hits = pd.concat([
        _hits(30.0),
        pd.DataFrame(
            [
                ("tgt", "tmpl", "tWeaker", "gWeaker", 0.0, 27.0, 500, 300.0, 100.0),
                ("tmpl", "tgt", "gWeaker", "tWeaker", 0.0, 27.0, 500, 300.0, 100.0),
            ],
            columns=HIT_COLUMNS,
        ),
    ])

    result = get_model_from_homology([template], hits, "tgt", review_identity=25)

    assert list(result.candidates.reaction) == ["WEAK", "WEAKER"]
    assert list(result.candidates.identity) == pytest.approx([30.0, 27.0])


def test_review_identity_above_min_identity_is_rejected():
    """It exists to catch what min_identity turns away; above it, it catches nothing."""
    with pytest.raises(ValueError, match="must be below min_identity"):
        get_model_from_homology(
            [_template()], _hits(30.0), "tgt", min_identity=40, review_identity=45
        )


def test_reported_identity_is_the_direction_that_did_the_rejecting():
    """A one-sided strong hit must not be reported as though it were the reason.

    Matching is bidirectional: a pair has to clear the threshold both ways. If
    the forward hit is 44 % and the reverse 30 %, the pair is rejected on the 30.
    Reporting the 44 would show a curator a comfortable match with no explanation
    of why it was turned away.
    """
    template = _template()
    hits = pd.DataFrame(
        [
            # STRONG transfers normally.
            ("tgt", "tmpl", "tStrong", "gStrong", 0.0, 90.0, 500, 400.0, 100.0),
            ("tmpl", "tgt", "gStrong", "tStrong", 0.0, 90.0, 500, 400.0, 100.0),
            # WEAK looks comfortable one way and marginal the other.
            ("tgt", "tmpl", "tWeak", "gWeak", 0.0, 44.0, 500, 400.0, 100.0),
            ("tmpl", "tgt", "gWeak", "tWeak", 0.0, 30.0, 500, 300.0, 100.0),
        ],
        columns=HIT_COLUMNS,
    )

    result = get_model_from_homology([template], hits, "tgt", review_identity=25)

    assert {r.id for r in result.model.reactions} == {"STRONG"}
    assert list(result.candidates.reaction) == ["WEAK"]
    assert result.candidates.iloc[0].identity == pytest.approx(30.0), (
        "the limiting direction, not the flattering one"
    )
