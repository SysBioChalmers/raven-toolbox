"""Tests for rename_model_genes (curation/rename_genes.py, renameModelGenes port)."""
import cobra
import pandas as pd
import pytest

from raven_toolbox.curation import RenameGenesResult, rename_model_genes


@pytest.fixture
def model():
    m = cobra.Model("t")
    m.add_metabolites([cobra.Metabolite("a", compartment="c"), cobra.Metabolite("b", compartment="c")])
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({m.metabolites.a: -1, m.metabolites.b: 1})
    r1.gene_reaction_rule = "gene1 or gene10"
    m.add_reactions([r1])
    return m


def test_renames_genes_and_gpr(model):
    table = pd.DataFrame({"locus_tag": ["gene1", "gene10"], "gene_name": ["ALPHA1", "ALPHA10"]})
    result = rename_model_genes(model, table, "locus_tag", "gene_name")
    assert isinstance(result, RenameGenesResult)
    assert result.renamed == 2
    assert result.unmapped == []
    assert {g.id for g in model.genes} == {"ALPHA1", "ALPHA10"}
    # "gene1" must not have been matched inside "gene10" -- cobra's AST-based
    # renamer parses gene tokens rather than substring-searching, so it
    # can't suffer that failure mode in the first place.
    r1 = model.reactions.get_by_id("r1")
    assert set(g.id for g in r1.genes) == {"ALPHA1", "ALPHA10"}


def test_unmapped_genes_left_unchanged_and_reported(model):
    table = pd.DataFrame({"locus_tag": ["gene1"], "gene_name": ["ALPHA1"]})
    result = rename_model_genes(model, table, "locus_tag", "gene_name")
    assert result.renamed == 1
    assert result.unmapped == ["gene10"]
    assert {g.id for g in model.genes} == {"ALPHA1", "gene10"}


def test_first_nonempty_mapping_wins_for_duplicates(model):
    table = pd.DataFrame({
        "locus_tag": ["gene1", "gene1", "gene10"],
        "gene_name": ["FIRST", "SECOND", "ALPHA10"],
    })
    result = rename_model_genes(model, table, "locus_tag", "gene_name")
    assert result.renamed == 2
    assert {g.id for g in model.genes} == {"FIRST", "ALPHA10"}


def test_rows_with_empty_or_nan_target_are_skipped(model):
    table = pd.DataFrame({
        "locus_tag": ["gene1", "gene10"],
        "gene_name": ["", float("nan")],
    })
    result = rename_model_genes(model, table, "locus_tag", "gene_name")
    assert result.renamed == 0
    assert sorted(result.unmapped) == ["gene1", "gene10"]
    assert {g.id for g in model.genes} == {"gene1", "gene10"}


def test_gene_table_accepts_a_tsv_path(model, tmp_path):
    path = tmp_path / "genes.tsv"
    pd.DataFrame({"locus_tag": ["gene1", "gene10"], "gene_name": ["ALPHA1", "ALPHA10"]}).to_csv(
        path, sep="\t", index=False,
    )
    result = rename_model_genes(model, path, "locus_tag", "gene_name")
    assert result.renamed == 2
    assert {g.id for g in model.genes} == {"ALPHA1", "ALPHA10"}


def test_missing_column_raises(model):
    table = pd.DataFrame({"locus_tag": ["gene1"], "gene_name": ["ALPHA1"]})
    with pytest.raises(ValueError, match="not found in gene_table"):
        rename_model_genes(model, table, "locus_tag", "nonexistent_column")


def test_renaming_two_genes_to_an_existing_one_merges_them(model):
    """RAVEN's own regex-based renamer has no check for this and can leave
    a model with two model.genes entries sharing the same id; cobra's
    rename_genes explicitly merges the underlying Gene objects instead (the
    GPR string still references it twice, e.g. "SAME or SAME" -- renaming
    doesn't algebraically simplify a rule, only the gene identity is
    deduplicated) -- a real, disclosed improvement, not a like-for-like
    behavioural match."""
    table = pd.DataFrame({"locus_tag": ["gene1", "gene10"], "gene_name": ["SAME", "SAME"]})
    rename_model_genes(model, table, "locus_tag", "gene_name")
    assert [g.id for g in model.genes] == ["SAME"]
    r1 = model.reactions.get_by_id("r1")
    assert r1.gene_reaction_rule == "SAME or SAME"
