"""Smoke tests for the ravenpy scaffold.

Real functionality tests are added per-function as the port proceeds (see PLAN.md).
"""

import importlib


def test_package_imports():
    import ravenpy

    assert ravenpy.__version__


def test_subpackages_importable():
    for sub in (
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
        assert importlib.import_module(f"ravenpy.{sub}")


def test_cobra_available():
    import cobra

    assert cobra.__version__
