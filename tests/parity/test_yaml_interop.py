"""Tier 1: raven-toolbox must read what MATLAB RAVEN writes.

RAVEN's repository ships models written by its own ``writeYAMLmodel``. Reading
them here is a real cross-language check that needs no MATLAB installation --
and it is the check that would have caught raven-toolbox failing on any RAVEN
model without genes (``KeyError: 'genes'``), which no unit test covered because
every fixture in this repository happened to have genes.

Exact here means semantically exact. Both implementations write valid YAML but
differ in key order and quoting; comparing bytes would test serialiser habits
rather than the model.
"""
from __future__ import annotations

import pytest

from raven_toolbox.io import read_yaml_model, write_yaml_model

pytestmark = pytest.mark.parity


def test_reads_every_raven_authored_model(raven_models):
    """Every YAML model RAVEN ships must load."""
    failures = {}
    for name, path in raven_models.items():
        try:
            read_yaml_model(path)
        except Exception as exc:  # noqa: BLE001 - reported per file below
            failures[name] = f"{type(exc).__name__}: {exc}"

    assert not failures, "RAVEN-written models that raven-toolbox cannot read: " + repr(
        failures
    )


def test_models_are_non_trivial(raven_models):
    """Guard against the suite passing because it read nothing.

    A reader that silently produced empty models would satisfy the test above.
    """
    sizes = {}
    for name, path in raven_models.items():
        model = read_yaml_model(path)
        sizes[name] = (len(model.reactions), len(model.metabolites))

    assert any(r > 0 and m > 0 for r, m in sizes.values()), (
        f"no RAVEN model produced any content: {sizes}"
    )


def test_round_trip_is_stable(raven_models, tmp_path):
    """Writing a RAVEN model and reading it back preserves the model.

    Not byte-for-byte against RAVEN's file -- key order and quoting differ
    without meaning -- but the content a model *is* must survive intact.
    """
    for name, path in raven_models.items():
        original = read_yaml_model(path)
        out = tmp_path / name
        write_yaml_model(original, out)
        again = read_yaml_model(out)

        assert [r.id for r in again.reactions] == [
            r.id for r in original.reactions
        ], f"{name}: reaction ids changed"
        assert [m.id for m in again.metabolites] == [
            m.id for m in original.metabolites
        ], f"{name}: metabolite ids changed"
        assert [g.id for g in again.genes] == [
            g.id for g in original.genes
        ], f"{name}: gene ids changed"

        for before, after in zip(original.reactions, again.reactions, strict=True):
            assert after.bounds == before.bounds, f"{name}: {before.id} bounds changed"
            assert (
                after.gene_reaction_rule == before.gene_reaction_rule
            ), f"{name}: {before.id} GPR changed"
            assert {m.id: c for m, c in after.metabolites.items()} == {
                m.id: c for m, c in before.metabolites.items()
            }, f"{name}: {before.id} stoichiometry changed"


def test_raven_specific_fields_survive(raven_models, tmp_path):
    """RAVEN's own per-entry fields are not dropped on the way through.

    They have no cobra counterpart and live in ``notes``; losing them silently
    would make raven-toolbox a lossy step in a MATLAB workflow.
    """
    interesting = {"inchis", "deltaG", "metFrom", "rxnFrom", "confidence_score"}

    for name, path in raven_models.items():
        original = read_yaml_model(path)
        present = {
            key
            for entity in (*original.metabolites, *original.reactions)
            for key in (entity.notes or {})
        } & interesting
        if not present:
            continue

        out = tmp_path / name
        write_yaml_model(original, out)
        again = read_yaml_model(out)

        survived = {
            key
            for entity in (*again.metabolites, *again.reactions)
            for key in (entity.notes or {})
        } & interesting
        assert present <= survived, (
            f"{name}: RAVEN fields lost in round-trip: {sorted(present - survived)}"
        )
