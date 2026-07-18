"""Tests for raven_toolbox.annotation (SBO terms + ΔG CSV persistence)."""
from __future__ import annotations

import math

import cobra
import pandas as pd

from raven_toolbox.annotation import (
    DEFAULT_BIOMASS_MET_NAMES,
    add_sbo_terms,
    load_delta_g_csv,
    save_delta_g_csv,
)
from raven_toolbox.annotation.sbo import _default_transport_detector

# --- shared tiny model -------------------------------------------------

def _toy_model() -> cobra.Model:
    """Toy yeast-flavoured model with: extracellular exchange, transport,
    metabolic, biomass, NGAM, generic pseudo, lipid-backbone biomass met."""
    m = cobra.Model("toy")
    m.compartments = {"c": "cytoplasm", "e": "extracellular", "x": "peroxisome"}

    mets = {
        "atp_c":      cobra.Metabolite("atp_c",      name="ATP",      compartment="c", charge=-4, formula="C10H12N5O13P3"),
        "atp_x":      cobra.Metabolite("atp_x",      name="ATP",      compartment="x", charge=-4, formula="C10H12N5O13P3"),
        "glc_e":      cobra.Metabolite("glc_e",      name="glucose",  compartment="e", charge=0,  formula="C6H12O6"),
        "biomass":    cobra.Metabolite("biomass_c",  name="biomass",  compartment="c"),
        "lipid_bb":   cobra.Metabolite("lbb_c",      name="lipid backbone", compartment="c"),
    }
    m.add_metabolites(list(mets.values()))

    # Exchange (single met in extracellular)
    ex = cobra.Reaction("EX_glc", lower_bound=-10, upper_bound=1000)
    ex.add_metabolites({mets["glc_e"]: -1})
    # Transport (ATP_c ↔ ATP_x — same met name, two compartments)
    tr = cobra.Reaction("ATP_tx", lower_bound=-1000, upper_bound=1000)
    tr.add_metabolites({mets["atp_c"]: -1, mets["atp_x"]: 1})
    # Normal metabolic reaction
    met = cobra.Reaction("metab1", lower_bound=0, upper_bound=1000)
    met.add_metabolites({mets["glc_e"]: -1, mets["atp_c"]: 1})
    # Biomass pseudoreaction (must be last for the bug-compat test)
    bio = cobra.Reaction("biomass_rxn", lower_bound=0, upper_bound=1000)
    bio.name = "biomass pseudoreaction"
    bio.add_metabolites({mets["atp_c"]: -1, mets["biomass"]: 1})
    m.add_reactions([ex, tr, met, bio])
    return m


# --- add_sbo_terms -----------------------------------------------------

def test_default_biomass_names_set_is_reasonable():
    assert {"biomass", "DNA", "RNA", "protein", "lipid"} <= DEFAULT_BIOMASS_MET_NAMES


def test_simple_chemical_gets_simple_sbo():
    m = _toy_model()
    add_sbo_terms(m)
    assert m.metabolites.get_by_id("glc_e").annotation["sbo"] == "SBO:0000247"


def test_biomass_pseudo_metabolite_gets_biomass_sbo():
    m = _toy_model()
    add_sbo_terms(m)
    assert m.metabolites.get_by_id("biomass_c").annotation["sbo"] == "SBO:0000649"


def test_lipid_backbone_suffix_match():
    m = _toy_model()
    add_sbo_terms(m)
    assert m.metabolites.get_by_id("lbb_c").annotation["sbo"] == "SBO:0000649"


def test_exchange_gets_exchange_sbo():
    m = _toy_model()
    add_sbo_terms(m)
    assert m.reactions.get_by_id("EX_glc").annotation["sbo"] == "SBO:0000627"


def test_transport_gets_transport_sbo():
    m = _toy_model()
    add_sbo_terms(m)
    assert m.reactions.get_by_id("ATP_tx").annotation["sbo"] == "SBO:0000655"


def test_biomass_pseudoreaction_gets_pseudo_sbo_by_default():
    m = _toy_model()
    add_sbo_terms(m)
    assert m.reactions.get_by_id("biomass_rxn").annotation["sbo"] == "SBO:0000629"


def test_legacy_bug_flag_scopes_pseudo_to_last_reaction():
    """The yeast-GEM compat flag must restrict the pseudo override to
    the last reaction; if biomass_rxn is last, it still wins."""
    m = _toy_model()
    add_sbo_terms(m, only_last_reaction_for_pseudo=True)
    # biomass_rxn was added last → still gets the pseudo SBO.
    assert m.reactions.get_by_id("biomass_rxn").annotation["sbo"] == "SBO:0000629"


def test_legacy_bug_flag_skips_pseudo_when_not_last():
    """Reorder: a non-pseudo reaction last; the bug flag should leave
    the pseudoreaction with whatever non-pseudo SBO it already had."""
    m = _toy_model()
    # Swap last two reactions so the pseudo is NOT last.
    rxns = list(m.reactions)
    new_order = rxns[:-2] + rxns[-1:] + rxns[-2:-1]
    m.reactions._generate_index()  # noqa: SLF001 (rebuild before reorder)
    # Easier: copy bounds, set order via remove+re-add
    m.remove_reactions([rxns[-2], rxns[-1]])
    m.add_reactions(new_order[-2:])
    add_sbo_terms(m, only_last_reaction_for_pseudo=True)
    bio = m.reactions.get_by_id("biomass_rxn")
    # The biomass reaction has a single bounded reactant set (atp_c -1,
    # biomass +1) → 2 mets → falls to default SBO:0000176.
    assert bio.annotation["sbo"] == "SBO:0000176"


def test_fill_semantic_preserves_existing():
    m = _toy_model()
    rxn = m.reactions.get_by_id("metab1")
    rxn.annotation["sbo"] = "SBO:0009999"
    add_sbo_terms(m)
    assert rxn.annotation["sbo"] == "SBO:0009999"


def test_custom_transport_detector_overrides_default():
    m = _toy_model()
    add_sbo_terms(m, transport_detector=lambda _m: set())
    # ATP_tx is no longer flagged transport → falls through to default
    # (it has 2 metabolites so it is not exchange/sink/demand).
    assert m.reactions.get_by_id("ATP_tx").annotation["sbo"] == "SBO:0000176"


def test_default_transport_detector_finds_transport():
    m = _toy_model()
    transports = _default_transport_detector(m)
    assert "ATP_tx" in transports
    assert "EX_glc" not in transports


# --- ΔG CSV round-trip ------------------------------------------------

def test_save_then_load_round_trip(tmp_path):
    m = _toy_model()
    m.metabolites.get_by_id("atp_c").notes["deltaG"] = "-12.34"
    m.metabolites.get_by_id("glc_e").notes["deltaG"] = "0.5"

    csv = tmp_path / "met_dg.csv"
    written = save_delta_g_csv(m.metabolites, csv)
    assert written == len(m.metabolites)
    df = pd.read_csv(csv)
    assert list(df.columns) == ["Var1", "Var2"]

    # Reload on a fresh model
    fresh = _toy_model()
    stamped = load_delta_g_csv(fresh.metabolites, csv)
    assert stamped == 2
    assert fresh.metabolites.get_by_id("atp_c").notes["deltaG"] == "-12.34"
    assert fresh.metabolites.get_by_id("glc_e").notes["deltaG"] == "0.5"


def test_save_writes_nan_for_missing_notes(tmp_path):
    m = _toy_model()  # no notes set
    csv = tmp_path / "met_dg.csv"
    save_delta_g_csv(m.metabolites, csv)
    df = pd.read_csv(csv)
    assert len(df) == len(m.metabolites)
    assert df["Var2"].apply(lambda v: math.isnan(v)).all()


def test_load_skips_nan_rows(tmp_path):
    """A NaN in the CSV must NOT clobber the entity's existing note."""
    m = _toy_model()
    m.metabolites.get_by_id("atp_c").notes["deltaG"] = "preserved"

    csv = tmp_path / "met_dg.csv"
    pd.DataFrame(
        {"Var1": ["atp_c", "glc_e"], "Var2": [math.nan, 1.0]}
    ).to_csv(csv, index=False)

    load_delta_g_csv(m.metabolites, csv)
    assert m.metabolites.get_by_id("atp_c").notes["deltaG"] == "preserved"
    assert m.metabolites.get_by_id("glc_e").notes["deltaG"] == "1.0"


def test_load_skips_the_missing_sentinel(tmp_path):
    """The ΔG side-car tables write 10000000 for "no valid ΔG", and yeast-GEM's side-car carries it on
    777 of its 4102 reaction rows. Stamping it would record a physically impossible 10^7 kJ/mol as a
    measurement. yeast-GEM's own checkrxnDirection.m gates on the same value."""
    m = _toy_model()
    m.metabolites.get_by_id("atp_c").notes["deltaG"] = "preserved"

    csv = tmp_path / "met_dg.csv"
    pd.DataFrame({"Var1": ["atp_c", "glc_e"], "Var2": [10000000.0, 1.0]}).to_csv(csv, index=False)

    stamped = load_delta_g_csv(m.metabolites, csv)
    assert stamped == 1                                              # only the real value
    assert m.metabolites.get_by_id("atp_c").notes["deltaG"] == "preserved"
    assert m.metabolites.get_by_id("glc_e").notes["deltaG"] == "1.0"


def test_sentinel_skipping_can_be_opted_out_of(tmp_path):
    m = _toy_model()
    csv = tmp_path / "met_dg.csv"
    pd.DataFrame({"Var1": ["atp_c"], "Var2": [10000000.0]}).to_csv(csv, index=False)

    assert load_delta_g_csv(m.metabolites, csv, missing_value=None) == 1
    assert m.metabolites.get_by_id("atp_c").notes["deltaG"] == "10000000.0"


def test_sentinel_recognised_whatever_dtype_the_csv_round_trip_produces(tmp_path):
    """The same sentinel arrives as 10000000, 10000000.0 or "10000000.0" depending on the writer and
    the column's inferred dtype -- a string column appears as soon as one row holds text."""
    m = _toy_model()
    csv = tmp_path / "met_dg.csv"
    pd.DataFrame({"Var1": ["atp_c", "glc_e"], "Var2": ["10000000.0", "-2.5"]}).to_csv(csv, index=False)

    assert load_delta_g_csv(m.metabolites, csv) == 1
    assert "deltaG" not in m.metabolites.get_by_id("atp_c").notes
    assert m.metabolites.get_by_id("glc_e").notes["deltaG"] == "-2.5"


def test_custom_columns_and_note_key(tmp_path):
    m = _toy_model()
    m.metabolites.get_by_id("atp_c").notes["dG_kJ"] = "-30.5"
    csv = tmp_path / "out.csv"
    save_delta_g_csv(
        m.metabolites, csv,
        id_column="id", value_column="dG", note_key="dG_kJ",
    )
    df = pd.read_csv(csv)
    assert list(df.columns) == ["id", "dG"]

    fresh = _toy_model()
    load_delta_g_csv(
        fresh.metabolites, csv,
        id_column="id", value_column="dG", note_key="dG_kJ",
    )
    assert fresh.metabolites.get_by_id("atp_c").notes["dG_kJ"] == "-30.5"
