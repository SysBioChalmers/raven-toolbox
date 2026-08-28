"""Shared pytest fixtures for the raven_toolbox test-suite.

The linear-chain INIT model used by the tINIT/ftINIT scoring tests.
"""
from pathlib import Path

import cobra
import pytest

# cobra's FVA / sampling / gap-filling parallelise with multiprocessing by
# default. On Linux that uses fork; on macOS and Windows it uses spawn, which
# re-pickles the model into each worker and deadlocks inside pytest (observed:
# ftINIT's task-essentiality FVA hangs forever on the macOS/Windows runners).
# Force single-process for the whole test suite so it is deterministic and
# hang-free on every platform. Library behaviour for real users is unchanged.
cobra.Configuration().processes = 1


def _linear_chain_model(*, with_genes: bool = False) -> cobra.Model:
    """EX_A -> A -(r1)-> B -(r2)-> C -(r3)-> D.

    A is taken up via the reversible ``EX_A``; ``r1``/``r2`` are the productive path
    and ``r3`` the dead-end branch the INIT/scoring tests penalise. With ``with_genes``
    the three internal reactions get gene rules ``g1``/``g2``/``g3``.
    """
    m = cobra.Model("net")
    A, B, C, D = (
        cobra.Metabolite(x, name=x[:-2], compartment="c")
        for x in ("A_c", "B_c", "C_c", "D_c")
    )
    m.add_metabolites([A, B, C, D])
    exa = cobra.Reaction("EX_A", lower_bound=-1000, upper_bound=1000)
    exa.add_metabolites({A: -1})  # negative flux = uptake of A
    r1 = cobra.Reaction("r1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({A: -1, B: 1})
    r2 = cobra.Reaction("r2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({B: -1, C: 1})
    r3 = cobra.Reaction("r3", lower_bound=0, upper_bound=1000)
    r3.add_metabolites({C: -1, D: 1})
    m.add_reactions([exa, r1, r2, r3])
    if with_genes:
        for rid, rule in (("r1", "g1"), ("r2", "g2"), ("r3", "g3")):
            m.reactions.get_by_id(rid).gene_reaction_rule = rule
    return m


@pytest.fixture
def linear_chain_model() -> cobra.Model:
    """A fresh linear-chain INIT model (no gene rules)."""
    return _linear_chain_model()


@pytest.fixture
def linear_chain_model_with_genes() -> cobra.Model:
    """A fresh linear-chain INIT model with gene rules g1/g2/g3 on r1/r2/r3."""
    return _linear_chain_model(with_genes=True)


# --------------------------------------------------------------------------- #
# Synthetic KEGG dump
# --------------------------------------------------------------------------- #
# A small, entirely fictional KEGG-format dump used by the reconstruction/kegg
# tests. It mimics the flat-file *format* (so it exercises the parser's flag
# detection, overview-map skipping, InChI/formula handling, mapformula
# irreversibility, KO/gene grouping and taxonomy lineages) but contains no real
# KEGG content — identifiers, names, sequences and cross-references are all
# made up, since the project is not licensed to redistribute KEGG data. The org
# codes (aaa/bbb/ccc/zzz), gene ids (GENE0x), KO/reaction/compound numbers
# (K9xxxx/R9xxxx/C9xxxx) and protein sequences are invented; only the structural
# domain/phylum words and the comment keywords the parser keys off are generic.

_REACTION = """\
ENTRY       R90010                      Reaction
NAME        fictional glucohydrolase analogue
DEFINITION  Trehalike + Aqualike <=> 2 Glucolike
EQUATION    C91083 + C90001 <=> 2 C90031
ENZYME      9.9.9.99
PATHWAY     rn09500  Fictional sugar metabolism
            rn01199  Fictional global overview
MODULE      M90599  example module
ORTHOLOGY   K90001  fictional trehalase analogue [EC:9.9.9.99]
DBLINKS     RHEA: 99999
///
ENTRY       R90100                      Reaction
NAME        spontaneous example
COMMENT     This reaction is spontaneous.
EQUATION    C90002 <=> C90003
ORTHOLOGY   K90002  fictional dehydrogenase analogue
///
ENTRY       R90200                      Reaction
NAME        undefined stoich example
EQUATION    C90001 + n C90002 <=> C90003
///
ENTRY       R90300                      Reaction
NAME        general example
COMMENT     General reaction.
EQUATION    C90031 <=> C90006
ORTHOLOGY   K90009  lumped ortholog
///
ENTRY       R90400                      Reaction
NAME        empty after cancellation
EQUATION    C90007 <=> C90007
///
"""

_COMPOUND = """\
ENTRY       C90001                      Compound
NAME        Aqualike;
            Aquatwo
FORMULA     A2O
DBLINKS     PubChem: 9303
            ChEBI: 95377
///
ENTRY       C90002                      Compound
NAME        Atriplike
FORMULA     A10
///
ENTRY       C90003                      Compound
NAME        Endione;
            Endi
FORMULA     N21
///
ENTRY       C90006                      Compound
NAME        Endionep
FORMULA     N22
///
ENTRY       C90031                      Compound
NAME        Glucolike;
            Sweetgrain
FORMULA     G6
DBLINKS     ChEBI: 94167 97634
///
ENTRY       C91083                      Compound
NAME        Trehalike
FORMULA     T12
///
ENTRY       C90007                      Compound
NAME        Oxylike
FORMULA     O2
///
"""

# InChI for C90031 is a made-up string; only its "InChI=" prefix is asserted, and
# its presence makes the parser prefer it over the FORMULA (clearing the latter).
_COMPOUND_INCHI = "C90031\tInChI=1S/FAKE6/c1-2-3\n"

_KO = """\
ENTRY       K90001                      KO
NAME        treF, FCT
DEFINITION  fictional trehalase analogue [EC:9.9.9.99]
GENES       AAA: GENE01(alias1) GENE02
            CCC: GENE04 GENE05(alias5)
///
ENTRY       K90002                      KO
DEFINITION  fictional dehydrogenase analogue
GENES       BBB: GENE03
///
ENTRY       K90099                      KO
DEFINITION  unlinked ortholog
GENES       BBB: GENE98
///
"""

# Invented (non-biological) protein sequences; only bbb:GENE03's is asserted.
_GENES_PEP = """\
>aaa:GENE01 alias1; fictional protein one
MKKLLAVTACGGVTAGGVTAGGVTAGGVTAGGVTAGGVTAGGVTAGG
>aaa:GENE02 hypothetical protein
MKKLLAVTACGGVTAGGVTAGGVTAGGVTAGGVTAGGVTAGGVTAGX
>bbb:GENE03 fictional aspartokinase analogue
MQFKTLVIDEGHKLPSTWYNACRMQFKTLVIDEGHKLPSTWYNACR
>ccc:GENE04 fictional dehydrogenase a
MSSTTGGKKVVIICCKAAALLWWEELLKKPPFFSSIIEEVVEEVV
>ccc:GENE05 fictional dehydrogenase b
MSSTTGGKKVVIICCKAAALLWWEEVVKKPPFFSSIIEEDDVVEEVV
>zzz:GENE99 some other gene not in any KO
MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
"""

_TAXONOMY = (
    "# Prokaryotes\n"
    "## Bacteria\n"
    "### Firmicutes\n"
    "T90010\taaa\tFictobacter examplensis\tFictobacter\n"
    "### Gammaproteobacteria\n"
    "T90007\tbbb\tMockella testium\tMockella\n"
    "# Eukaryotes\n"
    "## Animals\n"
    "### Vertebrates\n"
    "T90001\tccc\tImaginaria animalis\tImaginaria\n"
)

_MAPFORMULA = (
    "R90010: 09500: C91083 => C90031\n"  # one direction in this map ...
    "R90010: 09010: C90031 => C91083\n"  # ... the opposite in another -> stays reversible
    "R90100: 09010: C90002 => C90003\n"  # only ever one direction -> irreversible
)


@pytest.fixture(scope="session")
def kegg_dump(tmp_path_factory) -> Path:
    """A synthetic, fully fictional KEGG-format dump (no real KEGG content).

    Exercises the KEGG dump parser end-to-end — reaction flags, overview-map
    skipping, InChI/formula handling, mapformula irreversibility, KO/gene
    grouping and taxonomy lineages — without redistributing any KEGG-derived
    data, which the project is not licensed to ship.
    """
    d = tmp_path_factory.mktemp("kegg_dump")
    (d / "reaction").write_text(_REACTION, encoding="utf-8")
    (d / "compound").write_text(_COMPOUND, encoding="utf-8")
    (d / "compound.inchi").write_text(_COMPOUND_INCHI, encoding="utf-8")
    (d / "ko").write_text(_KO, encoding="utf-8")
    (d / "genes.pep").write_text(_GENES_PEP, encoding="utf-8")
    (d / "reaction_mapformula.lst").write_text(_MAPFORMULA, encoding="utf-8")
    (d / "taxonomy").write_text(_TAXONOMY, encoding="utf-8")
    return d
