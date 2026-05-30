"""Sphinx configuration for the raven-python documentation.

Markdown sources are rendered by MyST-Parser, so the docs stay authored in the
same Markdown used on GitHub. The API reference is generated from the NumPy-style
docstrings via autodoc + autosummary (``:recursive:``).
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

# -- Make the package importable for autodoc --------------------------------
# On ReadTheDocs the package is pip-installed (see .readthedocs.yaml); this also
# lets a local `make html` work straight from a source checkout.
_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT / "src"))


def _get_version() -> str:
    try:
        from raven_python import __version__

        return __version__
    except Exception:  # pragma: no cover - docs must still build
        return "0.0.1"


# -- Project information -----------------------------------------------------
project = "raven-python"
author = "Eduard Kerkhoven"
copyright = f"{date.today().year}, SysBioChalmers"
release = _get_version()
version = release

# -- General configuration ---------------------------------------------------
extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.todo",
    "sphinx_copybutton",
    "sphinx_design",
]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "README.md"]

# Markdown + reStructuredText both allowed; .md goes through MyST.
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# -- MyST configuration ------------------------------------------------------
myst_enable_extensions = [
    "colon_fence",   # ::: fenced directives
    "deflist",
    "fieldlist",
    "tasklist",
    "smartquotes",
    "substitution",
]
myst_heading_anchors = 3  # auto-generate anchors for ##, ###, #### headings

# -- autodoc / autosummary ---------------------------------------------------
# The API reference is built from explicit per-subpackage ``automodule`` pages
# (docs/reference/api/*), so no autosummary stub generation is needed.
autosummary_generate = False
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
    "member-order": "bysource",
}
autodoc_typehints = "description"
autoclass_content = "both"
# Optional / heavy imports that need not be present for the docs to build.
autodoc_mock_imports = []

# -- napoleon (NumPy-style docstrings) --------------------------------------
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = False
napoleon_use_param = True
# Render docstring "Attributes:" sections as :ivar: fields rather than separate
# attribute directives — avoids duplicate-description warnings on the dataclasses
# (whose fields autodoc also documents).
napoleon_use_ivar = True

# -- intersphinx -------------------------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "cobra": ("https://cobrapy.readthedocs.io/en/latest/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

# -- todo extension ----------------------------------------------------------
todo_include_todos = True

# -- HTML output -------------------------------------------------------------
html_theme = "furo"
html_title = f"raven-python {release}"
html_static_path = ["_static"]
html_theme_options = {
    "source_repository": "https://github.com/SysBioChalmers/raven-python/",
    "source_branch": "develop",
    "source_directory": "docs/",
}

# Suppress noisy-but-harmless warnings (e.g. cross-references to GitHub source
# files that intentionally live outside the doc tree).
suppress_warnings = ["myst.xref_missing"]
