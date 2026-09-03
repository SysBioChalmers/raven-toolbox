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


def _model_with_compartments():
    m = _model()
    # A compartment only reaches the written file if some metabolite is
    # actually in it (cobra's own model_to_dict derives the section from
    # metabolite membership, not from model.compartments verbatim) --- so
    # an unused mitochondrial metabolite is added purely to make "m" appear
    # alongside "c". Insertion order is deliberately not alphabetical, so a
    # sorted-output assertion actually exercises the sort rather than
    # passing by luck.
    m.compartments = {"m": "mitochondria", "c": "cytosol"}
    m.add_metabolites([cobra.Metabolite("x_m", compartment="m")])
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


def test_write_yaml_sort_ids_keeps_compartments_omap_tag(tmp_path):
    """A regression test for a sort_ids=True-only bug: sorting compartments
    into a plain dict (rather than an OrderedDict) drops the !!omap tag
    ruamel would otherwise give it, so RAVEN's line-based readYAMLmodel.m
    --- keyed on that tag --- silently fails to recognise the section. The
    unsorted path (test_output_carries_omap_tags in test_io_yaml_parity.py)
    never exercised this, since it never took the sort_ids=True branch."""
    m = _model_with_compartments()
    out = tmp_path / "m.yml"
    write_yaml_model(m, out, sort_ids=True)
    text = out.read_text()

    assert "- compartments: !!omap" in text
    assert text.index("c: cytosol") < text.index("m: mitochondria")

    reloaded = read_yaml_model(out)
    assert reloaded.compartments == {"c": "cytosol", "m": "mitochondria"}
    assert {x.id for x in reloaded.metabolites} == {"a_c", "b_c", "x_m"}
