"""Tests for raven_toolbox.localization — predictor loaders + the MILP + apply (Phase 7)."""
from __future__ import annotations

from textwrap import dedent

import cobra
import pandas as pd
import pytest

from raven_toolbox.localization import (
    DEFAULT_COMPARTMENT_MAP,
    LocalizationProposal,
    LocalizationResult,
    LocalizationScores,
    apply_localization,
    load_compartments,
    load_deeploc,
    load_mulocdeep,
    load_uniprot,
    predict_localization,
)

# --------------------------------------------------------------------- loaders

def test_load_deeploc_csv(tmp_path):
    p = tmp_path / "deeploc.csv"
    p.write_text(dedent("""\
        Protein_ID,Localizations,Signals,Cytoplasm,Nucleus,Mitochondrion
        G1,Cytoplasm,,0.8,0.1,0.05
        G2,Mitochondrion,SP,0.05,0.15,0.9
    """))
    s = load_deeploc(p)
    # non-numeric metadata columns (Localizations, Signals) are dropped automatically
    assert set(s.compartments) == {"Cytoplasm", "Nucleus", "Mitochondrion"}
    # row-max → 1.0
    assert s.df.loc["G1", "Cytoplasm"] == pytest.approx(1.0)
    assert s.df.loc["G2", "Mitochondrion"] == pytest.approx(1.0)


def test_load_deeploc_with_compartment_map(tmp_path):
    p = tmp_path / "deeploc.csv"
    p.write_text(dedent("""\
        Protein_ID,Localizations,Signals,Cytoplasm,Nucleus,Mitochondrion,Plastid
        G1,Cytoplasm,,0.8,0.1,0.05,0.4
    """))
    s = load_deeploc(p, compartment_map=DEFAULT_COMPARTMENT_MAP)
    # predictor labels mapped to model ids; Plastid (unmapped, no fungal equivalent) dropped
    assert set(s.compartments) == {"c", "n", "m"}
    assert s.df.loc["G1", "c"] == pytest.approx(1.0)


def test_load_mulocdeep_wide(tmp_path):
    # MULocDeep-style wide table; tab-separated, id not necessarily first-named
    p = tmp_path / "muloc.tsv"
    p.write_text(dedent("""\
        protein\tCytoplasm\tMitochondrion\tPeroxisome
        G1\t0.2\t0.7\t0.1
        G2\t0.9\t0.05\t0.05
    """))
    s = load_mulocdeep(p, compartment_map=DEFAULT_COMPARTMENT_MAP)
    assert set(s.compartments) == {"c", "m", "p"}
    assert s.df.loc["G1", "m"] == pytest.approx(1.0)
    assert s.df.loc["G2", "c"] == pytest.approx(1.0)


def test_load_compartments_tsv(tmp_path):
    # COMPARTMENTS integrated_full layout: id, name, GO id, GO term name, confidence (0-5)
    p = tmp_path / "yeast_compartments_integrated_full.tsv"
    p.write_text(
        "YGL001C\tERG26\tGO:0005739\tMitochondrion\t4.5\n"
        "YGL001C\tERG26\tGO:0005829\tCytosol\t2.0\n"
        "YGL002W\tABC1\tGO:0005777\tPeroxisome\t3.0\n"
        "YGL002W\tABC1\tGO:0099999\tPlastid\t5.0\n"        # unmapped → dropped
    )
    s = load_compartments(p, compartment_map=DEFAULT_COMPARTMENT_MAP)
    assert set(s.compartments) == {"m", "c", "p"}
    assert s.df.loc["YGL001C", "m"] == pytest.approx(1.0)   # 4.5 normalised to row max
    assert s.df.loc["YGL001C", "c"] == pytest.approx(2.0 / 4.5)
    # min_confidence filters weak annotations: the 2.0 cytosol row is dropped, so no 'c' column
    s2 = load_compartments(p, compartment_map=DEFAULT_COMPARTMENT_MAP, min_confidence=3.0)
    assert "c" not in s2.compartments
    assert s2.df.loc["YGL001C", "m"] == pytest.approx(1.0)


def test_load_uniprot(tmp_path):
    # UniProtKB TSV export (accession / primary / ordered-locus / Subcellular location [CC])
    p = tmp_path / "uniprot.tsv"
    p.write_text(
        "Entry\tGene Names (primary)\tGene Names (ordered locus)\tSubcellular location [CC]\n"
        "P00890\tCIT1\tYNR001C\tSUBCELLULAR LOCATION: Mitochondrion matrix {ECO:0000269|PubMed:1}.\n"
        "P12345\tGENEX\tYAL001C\tSUBCELLULAR LOCATION: Cytoplasm. Nucleus {ECO:0000255}. "
        "Note=Shuttles to the mitochondrion under stress.\n"
        "P99999\tNOLOC\tYBR002C\t\n"
    )
    # use the ordered-locus column so ids match yeast-GEM ORF gene ids
    s = load_uniprot(p, id_column="Gene Names (ordered locus)")
    assert s.df.loc["YNR001C", "m"] == pytest.approx(1.0)
    assert s.df.loc["YAL001C", "c"] == pytest.approx(1.0)
    assert s.df.loc["YAL001C", "n"] == pytest.approx(1.0)
    # the "mitochondrion" mention inside Note=… is stripped, so YAL001C is not placed in m
    assert s.df.loc["YAL001C", "m"] == 0.0
    # a protein with no annotation is absent (treated as "no signal" downstream)
    assert "YBR002C" not in s.genes


def test_load_uniprot_autodetects_location_column(tmp_path):
    p = tmp_path / "u.tsv"
    p.write_text(
        "Entry\tSubcellular location [CC]\n"
        "P1\tSUBCELLULAR LOCATION: Peroxisome.\n"
    )
    s = load_uniprot(p)               # id_column defaults to first col, location col auto-detected
    assert s.df.loc["P1", "p"] == pytest.approx(1.0)


def test_fetch_uniprot_localization(monkeypatch):
    import urllib.request

    from raven_toolbox.localization import fetch_uniprot_localization

    tsv = (
        "Entry\tGene Names (primary)\tGene Names (ordered locus)\tSubcellular location [CC]\n"
        "P00890\tCIT1\tYNR001C\tSUBCELLULAR LOCATION: Mitochondrion matrix {ECO:1}.\n"
    )
    captured = {}

    class _Resp:
        def read(self):
            return tsv.encode("utf-8")

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    def _fake_urlopen(req, timeout=None):
        captured["url"] = req.full_url
        return _Resp()

    monkeypatch.setattr(urllib.request, "urlopen", _fake_urlopen)
    s = fetch_uniprot_localization(559292)            # default id_field=gene_oln → ORF ids
    assert s.df.loc["YNR001C", "m"] == pytest.approx(1.0)
    # query was built correctly (organism filter + reviewed-only by default)
    assert "organism_id" in captured["url"] and "559292" in captured["url"]
    assert "reviewed" in captured["url"]


def test_load_compartments_collapses_synonyms(tmp_path):
    # Lysosome and Vacuole both map to 'v' → collapsed by max
    p = tmp_path / "c.tsv"
    p.write_text(
        "g1\tn1\tGO:1\tVacuole\t2.0\n"
        "g1\tn1\tGO:2\tLysosome\t4.0\n"
    )
    s = load_compartments(p, compartment_map=DEFAULT_COMPARTMENT_MAP)
    assert s.compartments == ["v"]
    assert s.df.loc["g1", "v"] == pytest.approx(1.0)        # max(2,4)=4 → normalised 1.0


def test_localization_scores_with_compartments_rename():
    df = pd.DataFrame({"cyto": [1.0], "mito": [0.2]}, index=pd.Index(["g1"], name="gene_id"))
    s = LocalizationScores(df).with_compartments({"cyto": "c", "mito": "m"})
    assert list(s.compartments) == ["c", "m"]


# ----------------------------------------------------------------- predict (toy)

def _toy_two_compartment_model() -> cobra.Model:
    """Single-compartment draft (everything in 'c'):

    A_c -(r1)-> B_c -(r2)-> C_c        (r2 should move to 'm' per scores below)
    Boundary EX_A imports A; EX_C drains C.
    """
    m = cobra.Model("toy")
    A, B, C = (cobra.Metabolite(x + "_c", name=x, compartment="c") for x in "ABC")
    m.add_metabolites([A, B, C])

    def rxn(rid, lb, ub, mets, gpr=None):
        r = cobra.Reaction(rid, lower_bound=lb, upper_bound=ub)
        r.add_metabolites(mets)
        if gpr:
            r.gene_reaction_rule = gpr
        return r
    m.add_reactions([rxn("EX_A", -1000, 0, {A: -1}),
                     rxn("EX_C", 0, 1000, {C: -1}),
                     rxn("r1", 0, 1000, {A: -1, B: 1}, "g1"),
                     rxn("r2", 0, 1000, {B: -1, C: 1}, "g2")])
    return m


def test_predict_empty_relocate_set_is_no_op():
    """An empty relocate set short-circuits to an empty proposal."""
    m = _toy_two_compartment_model()
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0, 0.0], "m": [0.0, 1.0]}, index=pd.Index(["g1", "g2"], name="gene_id")))
    res = predict_localization(m, scores, reactions_to_relocate=[], apply=False)
    assert isinstance(res, LocalizationProposal)
    assert res.moved.empty


def test_predict_places_single_reaction():
    """Pass r2 in the relocate set; it goes to 'm' per scores."""
    m = _toy_two_compartment_model()
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0, 0.0], "m": [0.0, 1.0]}, index=pd.Index(["g1", "g2"], name="gene_id")))
    res = predict_localization(m, scores, ["r2"], default_compartment="c", apply=False,
                                transport_cost=0.1)
    assert isinstance(res, LocalizationProposal)
    assert set(res.moved["rxn_id"]) == {"r2"}
    assert res.moved.iloc[0]["to_compartment"] == "m"


def test_predict_apply_creates_compartment_metabolites_and_transports():
    """apply=True should mutate the (copy) model: r2 in m, and B/C transports added."""
    m = _toy_two_compartment_model()
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0, 0.0], "m": [0.0, 1.0]}, index=pd.Index(["g1", "g2"], name="gene_id")))
    res = predict_localization(m, scores, ["r2"], default_compartment="c", apply=True,
                                transport_cost=0.05)
    assert isinstance(res, LocalizationResult)
    out = res.model
    r2 = out.reactions.r2
    assert {mt.compartment for mt in r2.metabolites} == {"m"}   # both substrates now in m
    # B_m and C_m metabolite copies must exist:
    assert "B_m" in out.metabolites and "C_m" in out.metabolites
    # Transports tr_B_m and tr_C_m must be added (default c ↔ m):
    transport_ids = {t.id for t in res.added_transports}
    assert "tr_B_m" in transport_ids
    assert "tr_C_m" in transport_ids
    # Original model untouched (we copied).
    assert m.reactions.r2.metabolites != r2.metabolites
    assert "B_m" not in m.metabolites


def test_predict_unplaced_reaction_when_no_scored_gene():
    """A relocate-set reaction whose genes are all absent from scores is reported, not crashed."""
    m = _toy_two_compartment_model()
    # Only g1 has scores; g2 (r2's gene) is absent → r2 is unplaceable.
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0], "m": [0.0]}, index=pd.Index(["g1"], name="gene_id")))
    res = predict_localization(m, scores, ["r2"], apply=False)
    assert isinstance(res, LocalizationProposal)
    assert "r2" in res.unplaced_reactions
    assert res.moved.empty   # nothing actually placed


def test_predict_boundary_reactions_always_pinned():
    """Boundary reactions in the relocate set are silently filtered out."""
    m = _toy_two_compartment_model()
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0, 0.0], "m": [0.0, 1.0]}, index=pd.Index(["g1", "g2"], name="gene_id")))
    res = predict_localization(m, scores, ["EX_A", "EX_C", "r2"], apply=False,
                                transport_cost=0.1)
    # Only r2 should appear in the proposal — boundaries dropped.
    assert set(res.moved["rxn_id"]) == {"r2"}


def test_predict_default_compartment_validated():
    m = _toy_two_compartment_model()
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0], "m": [0.0]}, index=pd.Index(["g2"], name="gene_id")))
    with pytest.raises(ValueError, match="default_compartment"):
        predict_localization(m, scores, ["r2"], default_compartment="x", apply=False)


def test_apply_localization_idempotent_on_empty_proposal():
    """An empty proposal (no moves, no transports) shouldn't change the model."""
    m = _toy_two_compartment_model()
    empty = LocalizationProposal(
        moved=pd.DataFrame(columns=["rxn_id", "from_compartment", "to_compartment"]),
        added_transports=pd.DataFrame(columns=["met_id", "compartment"]),
        gene_compartments={})
    out, added = apply_localization(m, empty)
    assert len(out.reactions) == len(m.reactions)
    assert added == []


# ----------------------------------------- multi-compartment scoring (NEW)

def test_predict_multi_compartment_when_secondary_score_beats_penalty():
    """Dual-localised gene: secondary compartment score 0.8 > penalty 0.3 → gene lands in
    both compartments. Two reactions sharing one gene each placed in their best
    compartment without contradicting the gene assignment."""
    m = cobra.Model("dual")
    A_c = cobra.Metabolite("A_c", compartment="c")
    B_c = cobra.Metabolite("B_c", compartment="c")
    m.add_metabolites([A_c, B_c])
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({A_c: -1, B_c: 1})
    r1.gene_reaction_rule = "g_dual"
    r2 = cobra.Reaction("r2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({A_c: -1, B_c: 1})
    r2.gene_reaction_rule = "g_dual"
    m.add_reactions([r1, r2])
    # g_dual scores: c=1.0 (primary), m=0.8 (secondary). Penalty 0.3 — secondary worth it.
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0], "m": [0.8]}, index=pd.Index(["g_dual"], name="gene_id")))
    res = predict_localization(m, scores, ["r2"], default_compartment="c", apply=False,
                                transport_cost=0.0, multi_compartment_penalty=0.3)
    # The gene should land in BOTH c and m (primary free + 0.8 - 0.3 > 0 for secondary).
    assert set(res.gene_compartments["g_dual"]) == {"c", "m"}


def test_predict_mono_when_secondary_score_below_penalty():
    """Same as above but penalty 0.9 > secondary score 0.8 → gene mono-localises (c only)."""
    m = cobra.Model("mono")
    A_c = cobra.Metabolite("A_c", compartment="c")
    B_c = cobra.Metabolite("B_c", compartment="c")
    m.add_metabolites([A_c, B_c])
    r2 = cobra.Reaction("r2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({A_c: -1, B_c: 1})
    r2.gene_reaction_rule = "g_dual"
    m.add_reactions([r2])
    scores = LocalizationScores(pd.DataFrame(
        {"c": [1.0], "m": [0.8]}, index=pd.Index(["g_dual"], name="gene_id")))
    res = predict_localization(m, scores, ["r2"], default_compartment="c", apply=False,
                                transport_cost=0.0, multi_compartment_penalty=0.9)
    # Penalty exceeds secondary score → only the primary compartment.
    assert res.gene_compartments["g_dual"] == ["c"]


def test_predict_high_penalty_forces_mono_localisation():
    """Very high penalty effectively bans extra compartments."""
    m = _toy_two_compartment_model()
    scores = LocalizationScores(pd.DataFrame(
        {"c": [0.4, 0.3], "m": [0.5, 0.6]},
        index=pd.Index(["g1", "g2"], name="gene_id")))
    res = predict_localization(m, scores, ["r1", "r2"], default_compartment="c", apply=False,
                                transport_cost=0.0, multi_compartment_penalty=1000.0)
    # With penalty so high, every gene gets exactly one compartment.
    for g, comps in res.gene_compartments.items():
        assert len(comps) == 1, f"{g} landed in {comps} with prohibitive penalty"


# --------------------------------------------------- sequence input prep (DeepLoc 2.1 FASTA)
def _fake_urlopen_returning(text):
    class _Resp:
        def read(self):
            return text.encode("utf-8")

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    def _open(req, timeout=None):
        _open.url = req.full_url
        return _Resp()

    return _open


def test_write_fasta_basic(tmp_path):
    from raven_toolbox.localization import write_fasta
    p = tmp_path / "in.fasta"
    paths = write_fasta({"YAL001C": "MKVLAA", "YBR002W": "MQT"}, p, wrap=4)
    assert paths == [p]
    assert p.read_text() == ">YAL001C\nMKVL\nAA\n>YBR002W\nMQT\n"   # wrapped at 4 residues


def test_write_fasta_chunks_for_web_limit(tmp_path):
    from raven_toolbox.localization import write_fasta
    seqs = {f"G{i}": "AAAAAAAAAA" for i in range(5)}
    paths = write_fasta(seqs, tmp_path / "in.fasta", max_records_per_file=2)
    assert [p.name for p in paths] == ["in_001.fasta", "in_002.fasta", "in_003.fasta"]
    assert paths[0].read_text().count(">") == 2 and paths[2].read_text().count(">") == 1


def test_write_fasta_rejects_whitespace_id(tmp_path):
    from raven_toolbox.localization import write_fasta
    with pytest.raises(ValueError, match="whitespace"):
        write_fasta({"bad id": "AAAAAAAAAA"}, tmp_path / "x.fasta")


def test_fetch_protein_sequences(monkeypatch):
    import urllib.request

    from raven_toolbox.localization import fetch_protein_sequences
    tsv = (
        "Entry\tGene Names (primary)\tGene Names (ordered locus)\tSequence\n"
        "P1\tCIT1\tYNR001C\tMSEQUENCEAA\n"               # 11 aa -> kept
        "P2\tGENEB\tYBR002W YBR002W-A\tMQTPROTEINX\n"    # two ordered-locus tokens, both map
        "P3\tTINY\tYTI001C\tSHORT\n"                     # 5 aa -> dropped by min_length=10
    )
    monkeypatch.setattr(urllib.request, "urlopen", _fake_urlopen_returning(tsv))
    seqs = fetch_protein_sequences(559292, min_length=10)
    assert seqs["YNR001C"] == "MSEQUENCEAA"
    assert seqs["YBR002W"] == "MQTPROTEINX" and seqs["YBR002W-A"] == "MQTPROTEINX"
    assert "YTI001C" not in seqs                          # below DeepLoc's 10-aa minimum

    monkeypatch.setattr(urllib.request, "urlopen", _fake_urlopen_returning(tsv))
    only = fetch_protein_sequences(559292, genes=["YNR001C"], min_length=10)
    assert set(only) == {"YNR001C"}                       # restricted to requested genes


def test_prepare_deeploc_input_from_model(monkeypatch, tmp_path):
    import urllib.request

    from raven_toolbox.localization import prepare_deeploc_input
    tsv = (
        "Entry\tGene Names (primary)\tGene Names (ordered locus)\tSequence\n"
        "P1\tA\tYAL001C\tMAAAAAAAAAA\n"
        "P2\tB\tYBR002W\tMBBBBBBBBBB\n"
    )
    monkeypatch.setattr(urllib.request, "urlopen", _fake_urlopen_returning(tsv))
    # a cobra model exposes .genes; YDR999W has no UniProt entry -> reported missing, not dropped
    model = cobra.Model("toy")
    r = cobra.Reaction("r1")
    r.gene_reaction_rule = "YAL001C or YBR002W or YDR999W"
    model.add_reactions([r])

    res = prepare_deeploc_input(model, 559292, tmp_path / "yeast.fasta")
    assert res.n_requested == 3 and res.n_written == 2
    assert res.missing == ["YDR999W"]
    assert len(res.paths) == 1                            # 2 seqs < 500 web limit -> single file
    fasta = res.paths[0].read_text()
    # headers are the model's gene ids -> they become DeepLoc's Protein_ID and load_deeploc lines up
    assert ">YAL001C" in fasta and ">YBR002W" in fasta and "YDR999W" not in fasta
