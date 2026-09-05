"""Tests for replace_metabolite (manipulation/replace.py, replaceMets port)."""
import cobra
import pytest

from raven_toolbox.manipulation import replace_metabolite


def _model():
    m = cobra.Model("t")
    a = cobra.Metabolite("a", name="oxygen", compartment="c", formula="O2", charge=0)
    b = cobra.Metabolite("b", name="o2", compartment="c")
    c = cobra.Metabolite("c", compartment="c")
    m.add_metabolites([a, b, c])
    r1 = cobra.Reaction("r1", lower_bound=-1000, upper_bound=1000)
    r1.add_metabolites({a: -1, c: 1})
    m.add_reactions([r1])
    return m, a, b, c


# --- identifiers=True (id-based, always 1:1) ---------------------------

def test_identifiers_merges_stoichiometry_and_removes_old():
    m, a, b, c = _model()
    replace_metabolite(m, "a", "b", identifiers=True)
    assert "a" not in m.metabolites
    r1 = m.reactions.get_by_id("r1")
    assert r1.get_coefficient("b") == -1
    assert "c" not in {mm.id for mm in r1.metabolites if mm.id == "a"}


def test_identifiers_survivor_keeps_its_own_identity():
    # replaceMets.m copies the replacement's identity onto the old
    # metabolite's row, but that row is deleted immediately after (a dead
    # step) -- the survivor's own pre-existing identity is what's actually
    # observable afterward, unaffected by the old metabolite's identity.
    m, a, b, c = _model()
    a.formula = "not-o2-formula"
    a.annotation["kegg.compound"] = "WRONG"
    b.formula = "O2"
    b.annotation["kegg.compound"] = "C00007"
    replace_metabolite(m, "a", "b", identifiers=True)
    survivor = m.metabolites.get_by_id("b")
    assert survivor.formula == "O2"
    assert survivor.annotation.get("kegg.compound") == "C00007"


def test_identifiers_missing_metabolite_raises():
    m, a, b, c = _model()
    with pytest.raises(ValueError, match="cannot be found"):
        replace_metabolite(m, "nonexistent", "b", identifiers=True)


def test_identifiers_missing_replacement_raises():
    m, a, b, c = _model()
    with pytest.raises(ValueError, match="cannot be found"):
        replace_metabolite(m, "a", "nonexistent", identifiers=True)


def test_identifiers_removes_reactions_that_became_duplicates():
    m, a, b, c = _model()
    r2 = cobra.Reaction("r2", lower_bound=-1000, upper_bound=1000)
    r2.add_metabolites({b: -1, c: 1})  # identical to r1 once a -> b
    m.add_reactions([r2])
    removed = replace_metabolite(m, "a", "b", identifiers=True)
    assert set(removed) == {"r2"}
    assert "r1" in m.reactions
    assert "r2" not in m.reactions


def test_identifiers_leaves_unrelated_preexisting_duplicates_alone():
    m, a, b, c = _model()
    d = cobra.Metabolite("d", compartment="c")
    m.add_metabolites([d])
    unrelated1 = cobra.Reaction("unrelated1", lower_bound=-1000, upper_bound=1000)
    unrelated1.add_metabolites({d: -1})
    unrelated2 = cobra.Reaction("unrelated2", lower_bound=-1000, upper_bound=1000)
    unrelated2.add_metabolites({d: -1})  # duplicate of unrelated1, untouched by a->b
    m.add_reactions([unrelated1, unrelated2])
    removed = replace_metabolite(m, "a", "b", identifiers=True)
    assert removed == []
    assert {"unrelated1", "unrelated2"} <= {r.id for r in m.reactions}


# --- identifiers=False (name-based) -------------------------------------

def test_name_based_merges_same_compartment_duplicates():
    m, a, b, c = _model()
    removed_names = {mm.name for mm in m.metabolites}
    assert removed_names == {"oxygen", "o2", ""}
    replace_metabolite(m, "oxygen", "o2")
    # a (renamed to "o2") and b share (name="o2", compartment="c") -> merged
    assert len(m.metabolites) == 2
    survivor = next(mm for mm in m.metabolites if mm.name == "o2")
    r1 = m.reactions.get_by_id("r1")
    assert r1.get_coefficient(survivor) == -1


def test_name_based_removes_reactions_that_became_duplicates():
    """Regression test: the replacement metabolite itself can be the one
    merged away (the first-encountered metabolite in the (name, compartment)
    group survives, which need not be `replacement`'s own object) -- the
    post-merge duplicate-reaction cleanup must key off the actual survivor,
    not a stale reference to a metabolite that's just been deleted."""
    m = cobra.Model("t")
    a = cobra.Metabolite("a", name="oxygen", compartment="c")
    b = cobra.Metabolite("b", name="o2", compartment="c")
    x = cobra.Metabolite("x", compartment="c")
    m.add_metabolites([a, b, x])
    r1 = cobra.Reaction("r1", lower_bound=-1000, upper_bound=1000)
    r1.add_metabolites({a: -1, x: 1})
    r2 = cobra.Reaction("r2", lower_bound=-1000, upper_bound=1000)
    r2.add_metabolites({b: -1, x: 1})  # identical to r1 once a/b are merged
    m.add_reactions([r1, r2])

    removed = replace_metabolite(m, "oxygen", "o2")

    assert set(removed) == {"r2"}
    assert "r1" in m.reactions
    assert "r2" not in m.reactions


def test_name_based_rename_without_collision_just_renames():
    m = cobra.Model("t")
    a = cobra.Metabolite("a", name="oxygen", compartment="c")
    b = cobra.Metabolite("b", name="o2", compartment="m")  # different compartment
    m.add_metabolites([a, b])
    replace_metabolite(m, "oxygen", "o2")
    assert len(m.metabolites) == 2  # no collision -> no merge
    assert m.metabolites.get_by_id("a").name == "o2"
    assert m.metabolites.get_by_id("a").compartment == "c"


def test_name_based_missing_metabolite_raises():
    m, a, b, c = _model()
    with pytest.raises(ValueError, match="cannot be found"):
        replace_metabolite(m, "nonexistent", "o2")


def test_name_based_missing_replacement_raises():
    m, a, b, c = _model()
    with pytest.raises(ValueError, match="cannot be found"):
        replace_metabolite(m, "oxygen", "nonexistent")


def test_verbose_does_not_crash(capsys):
    m, a, b, c = _model()
    replace_metabolite(m, "a", "b", identifiers=True, verbose=True)
    captured = capsys.readouterr()
    assert "r1" in captured.out
