"""Smoke tests for the ravengem scaffold.

Real functionality tests are added per-function as the port proceeds (see PLAN.md).
"""

import importlib


def test_package_imports():
    import ravengem

    assert ravengem.__version__


def test_subpackages_importable():
    for sub in (
        "manipulation",
        "io",
        "reconstruction",
        "reconstruction.kegg",
        "reconstruction.metacyc",
        "reconstruction.homology",
        "init",
        "tasks",
        "gapfilling",
        "omics",
        "localization",
        "comparison",
        "analysis",
        "plotting",
        "utils",
    ):
        assert importlib.import_module(f"ravengem.{sub}")


def test_cobra_available():
    import cobra

    assert cobra.__version__
