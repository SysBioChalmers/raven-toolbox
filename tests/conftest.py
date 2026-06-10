"""Shared pytest fixtures for the raven_python test-suite.

Currently the linear-chain INIT model used by the tINIT/ftINIT scoring tests,
which several modules previously built independently and identically.
"""
import cobra
import pytest


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
