"""Tests for raven_toolbox.io.git (exportForGit port)."""
import cobra
import pytest

from raven_toolbox.io import export_for_git
from raven_toolbox.manipulation import add_reactions_from_equations


@pytest.fixture
def model():
    m = cobra.Model("yeastGEM")
    m.compartments = {"c": "cytoplasm"}
    m.add_metabolites(
        [cobra.Metabolite("atp_c", name="ATP", compartment="c"),
         cobra.Metabolite("adp_c", name="ADP", compartment="c")]
    )
    add_reactions_from_equations(m, [{"id": "R1", "equation": "atp_c <=> adp_c"}])
    return m


def test_standard_gem_layout(model, tmp_path):
    root = export_for_git(model, tmp_path, prefix="yeast", formats=("yml", "xml", "mat", "xlsx", "txt"))
    assert root == tmp_path / "model"
    assert (root / "yml" / "yeast.yml").exists()
    assert (root / "xml" / "yeast.xml").exists()
    assert (root / "mat" / "yeast.mat").exists()
    assert (root / "xlsx" / "yeast.xlsx").exists()
    assert (root / "txt" / "yeast.txt").exists()
    assert (root / "dependencies.txt").exists()


def test_dependencies_file(model, tmp_path):
    root = export_for_git(model, tmp_path, formats=("yml",))
    deps = (root / "dependencies.txt").read_text()
    assert "python\t" in deps
    assert "cobra\t" in deps
    assert "raven_toolbox\t" in deps


def test_flat_layout(model, tmp_path):
    root = export_for_git(model, tmp_path, formats=("yml",), sub_dirs=False)
    assert root == tmp_path
    assert (tmp_path / "model.yml").exists()


def test_subset_of_formats(model, tmp_path):
    root = export_for_git(model, tmp_path, formats=("yml", "xml"))
    assert (root / "yml" / "model.yml").exists()
    assert not (root / "mat").exists()
    assert not (root / "xlsx").exists()


def test_does_not_mutate_model(model, tmp_path):
    order_before = [r.id for r in model.reactions]
    export_for_git(model, tmp_path, formats=("yml",))
    assert [r.id for r in model.reactions] == order_before


def test_txt_table_content(model, tmp_path):
    root = export_for_git(model, tmp_path, formats=("txt",))
    txt = (root / "txt" / "model.txt").read_text()
    assert txt.splitlines()[0].startswith("Rxn name\t")
    assert "R1" in txt
    assert "ATP[c]" in txt


def test_bad_format(model, tmp_path):
    with pytest.raises(ValueError, match="Unknown format"):
        export_for_git(model, tmp_path, formats=("yml", "json"))


def test_mat_export_with_subsystem_lists(tmp_path):
    """A model read from YAML keeps subsystems as lists; the MATLAB export still
    works and writes them as strings, without changing the model passed in."""
    m = cobra.Model("m")
    a = cobra.Metabolite("a_c", compartment="c")
    b = cobra.Metabolite("b_c", compartment="c")
    r1, r2 = cobra.Reaction("r1"), cobra.Reaction("r2")
    r1.add_metabolites({a: -1, b: 1})
    r2.add_metabolites({b: -1})
    m.add_reactions([r1, r2])
    r1.subsystem = ["Glycolysis"]
    r2.subsystem = ["Transport", "Exchange"]
    export_for_git(m, tmp_path, prefix="m", formats=("mat",), sub_dirs=False)
    back = cobra.io.load_matlab_model(str(tmp_path / "m.mat"))
    assert back.reactions.get_by_id("r1").subsystem == "Glycolysis"
    assert back.reactions.get_by_id("r2").subsystem == "Transport;Exchange"
    assert r1.subsystem == ["Glycolysis"]
