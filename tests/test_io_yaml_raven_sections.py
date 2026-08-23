"""RAVEN MATLAB omits empty top-level sections; the reader must tolerate that.

``writeYAMLmodel`` writes no ``genes:`` block at all for a model without genes
(RAVEN's own ``tutorial/small.yml`` is such a model). cobra's
``model_from_dict`` indexes ``obj["genes"]`` directly, so those files used to
fail to load with ``KeyError: 'genes'`` -- raven-toolbox could not read a valid
model written by the toolbox it is the counterpart of.
"""
from __future__ import annotations

import textwrap

import pytest

from raven_toolbox.io import read_yaml_model, write_yaml_model

# A minimal RAVEN-shaped model with no genes, and therefore no ``genes:``
# section. Authored here rather than copied from RAVEN, which is GPL.
NO_GENES = textwrap.dedent(
    """\
    ---
    - metaData: !!omap
      - id: tiny
      - name: No-gene test model
    - metabolites:
        - !!omap
          - id: a_c
          - name: A
          - compartment: c
        - !!omap
          - id: b_c
          - name: B
          - compartment: c
    - reactions:
        - !!omap
          - id: R1
          - name: A to B
          - metabolites: !!omap
              - a_c: -1
              - b_c: 1
          - lower_bound: 0
          - upper_bound: 1000
    - compartments: !!omap
      - c: cytoplasm
    """
)


@pytest.fixture
def no_gene_file(tmp_path):
    path = tmp_path / "tiny.yml"
    path.write_text(NO_GENES, encoding="utf-8")
    return path


def test_reads_model_without_genes_section(no_gene_file):
    model = read_yaml_model(no_gene_file)

    assert [r.id for r in model.reactions] == ["R1"]
    assert {m.id for m in model.metabolites} == {"a_c", "b_c"}
    assert len(model.genes) == 0


def test_round_trips_model_without_genes_section(no_gene_file, tmp_path):
    model = read_yaml_model(no_gene_file)
    out = tmp_path / "rewritten.yml"
    write_yaml_model(model, out)

    again = read_yaml_model(out)
    assert [r.id for r in again.reactions] == [r.id for r in model.reactions]
    assert [m.id for m in again.metabolites] == [m.id for m in model.metabolites]
    assert again.reactions.R1.bounds == model.reactions.R1.bounds


def test_model_with_no_entity_sections_at_all(tmp_path):
    """A skeleton model -- metadata and compartments only -- still loads.

    Only sections that are *empty* are omitted by RAVEN. A file with reactions
    but no ``metabolites:`` is not a RAVEN file, it is a broken one, and is
    expected to fail loudly rather than load a model whose reactions reference
    metabolites that do not exist.
    """
    path = tmp_path / "skeleton.yml"
    path.write_text(
        textwrap.dedent(
            """\
            ---
            - metaData: !!omap
              - id: skeleton
            - compartments: !!omap
              - c: cytoplasm
            """
        ),
        encoding="utf-8",
    )

    model = read_yaml_model(path)

    assert len(model.reactions) == 0
    assert len(model.metabolites) == 0
    assert len(model.genes) == 0
