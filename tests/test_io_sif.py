"""Tests for ravengem.io.sif (exportModelToSIF port)."""
import cobra
import pytest

from ravengem.io import export_model_to_sif
from ravengem.manipulation import add_reactions_from_equations


@pytest.fixture
def model():
    m = cobra.Model("t")
    m.add_metabolites([cobra.Metabolite(x, compartment="c") for x in ("a", "b", "c")])
    add_reactions_from_equations(
        m,
        [
            {"id": "R1", "equation": "a --> b"},
            {"id": "R2", "equation": "b --> c"},
        ],
    )
    return m


def _lines(path):
    return [ln.split("\t") for ln in path.read_text().splitlines()]


def test_reaction_compound(model, tmp_path):
    out = tmp_path / "g.sif"
    export_model_to_sif(model, out, "rc")
    rows = {r[0]: (r[1], set(r[2:])) for r in _lines(out)}
    assert rows["R1"] == ("rc", {"a", "b"})
    assert rows["R2"] == ("rc", {"b", "c"})


def test_reaction_reaction(model, tmp_path):
    out = tmp_path / "g.sif"
    export_model_to_sif(model, out, "rr")
    rows = {r[0]: set(r[2:]) for r in _lines(out)}
    # R1 and R2 share metabolite b
    assert rows["R1"] == {"R2"}
    assert rows["R2"] == {"R1"}


def test_compound_compound(model, tmp_path):
    out = tmp_path / "g.sif"
    export_model_to_sif(model, out, "cc")
    rows = {r[0]: set(r[2:]) for r in _lines(out)}
    # a is a substrate of R1 (a->b): a links to product b
    assert "b" in rows.get("a", set())
    # b is substrate of R2 (b->c): b links to c
    assert "c" in rows.get("b", set())


def test_custom_labels(model, tmp_path):
    out = tmp_path / "g.sif"
    export_model_to_sif(model, out, "rc", reaction_labels={"R1": "Reaction1"})
    sources = {r[0] for r in _lines(out)}
    assert "Reaction1" in sources
    assert "R1" not in sources


def test_bad_graph_type(model, tmp_path):
    with pytest.raises(ValueError, match="graph_type"):
        export_model_to_sif(model, tmp_path / "g.sif", "xx")


def test_cc_does_not_mutate_input(model, tmp_path):
    n_before = len(model.reactions)
    export_model_to_sif(model, tmp_path / "g.sif", "cc")
    assert len(model.reactions) == n_before  # convert_to_irreversible ran on a copy


# --- regression: label-map collision (known_issues.md B4) ------------------

def test_collapsing_label_map_warns(model, tmp_path):
    """A label map that sends two distinct ids to the same label silently merges
    nodes during the target-side dedup. Now warns so the user sees it."""
    with pytest.warns(UserWarning, match="multiple ids to the same label"):
        export_model_to_sif(
            model, tmp_path / "g.sif", "rc",
            reaction_labels={"R1": "shared", "R2": "shared"},
        )
