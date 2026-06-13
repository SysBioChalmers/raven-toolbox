"""Tests for sort_identifiers and write_yaml_model(sort_ids=True)."""
import cobra

from raven_toolbox.io import read_yaml_model, write_yaml_model
from raven_toolbox.manipulation import add_reactions_from_equations
from raven_toolbox.utils import sort_identifiers


def _model():
    m = cobra.Model("t")
    m.add_metabolites([cobra.Metabolite(x, compartment="c") for x in ("b_c", "a_c")])
    add_reactions_from_equations(
        m,
        [
            {"id": "R2", "equation": "a_c --> b_c", "gene_reaction_rule": "GB"},
            {"id": "R1", "equation": "b_c --> a_c", "gene_reaction_rule": "GA"},
        ],
    )
    return m


def test_sort_identifiers_orders_everything():
    m = _model()
    sort_identifiers(m)
    assert [r.id for r in m.reactions] == ["R1", "R2"]
    assert [x.id for x in m.metabolites] == ["a_c", "b_c"]
    assert [g.id for g in m.genes] == ["GA", "GB"]
    # lookup index still intact after sorting
    assert m.reactions.get_by_id("R2").id == "R2"


def test_write_yaml_sort_ids_does_not_mutate(tmp_path):
    m = _model()
    order_before = [r.id for r in m.reactions]
    out = tmp_path / "m.yml"
    write_yaml_model(m, out, sort_ids=True)
    assert [r.id for r in m.reactions] == order_before  # model untouched
    # but the file is sorted
    text = out.read_text()
    assert text.index("R1") < text.index("R2")
    reloaded = read_yaml_model(out)
    assert [r.id for r in reloaded.reactions] == ["R1", "R2"]
