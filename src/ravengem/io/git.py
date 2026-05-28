"""Export a model into a Standard-GEM versioned-repository layout.

Writes the model in several formats into the Standard-GEM folder structure (a
``model/`` directory with one subfolder per format), ready to commit to a
Git-maintained model repository (Metabolic Atlas / Human-GEM / yeast-GEM style),
plus a ``dependencies.txt`` recording tool versions.

Thin orchestration over the writers ravengem already exposes: ``write_yaml_model``,
cobra's ``write_sbml_model`` and ``save_matlab_model``, ``export_to_excel``, plus a
single-file reaction table (txt).
"""
from __future__ import annotations

import importlib.metadata as _md
import platform
from collections.abc import Iterable
from pathlib import Path

import cobra

from ravengem.io.excel import _equation, export_to_excel
from ravengem.io.yaml import write_yaml_model
from ravengem.utils.sort import sort_identifiers

_ALL_FORMATS = ("yml", "xml", "mat", "xlsx", "txt")


def _version(package: str) -> str:
    try:
        return _md.version(package)
    except _md.PackageNotFoundError:
        return "unknown"


def _write_txt(model: cobra.Model, path: Path) -> None:
    """Single-file, human-readable reaction table (RAVEN exportForGit txt)."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("Rxn name\tFormula\tGene-reaction association\tLB\tUB\tObjective\n")
        for r in model.reactions:
            fh.write(
                f"{r.id}\t{_equation(r)}\t{r.gene_reaction_rule}\t"
                f"{r.lower_bound:g}\t{r.upper_bound:g}\t{r.objective_coefficient:g}\n"
            )


def export_for_git(
    model: cobra.Model,
    path: str | Path = ".",
    *,
    prefix: str = "model",
    formats: Iterable[str] = ("yml", "xml", "mat", "xlsx"),
    sub_dirs: bool = True,
) -> Path:
    """Write ``model`` into a Standard-GEM repository layout.

    Parameters
    ----------
    path
        Directory to populate.
    prefix
        Base filename for every format (default ``"model"``).
    formats
        Which formats to write; any of ``"yml"``, ``"xml"``, ``"mat"``,
        ``"xlsx"``, ``"txt"`` (default ``yml``/``xml``/``mat``/``xlsx``).
    sub_dirs
        If True (default), write ``model/<fmt>/<prefix>.<fmt>`` (standard-GEM
        layout); otherwise all files go directly in ``path``.

    Returns
    -------
    pathlib.Path
        The root directory written to.
    """
    formats = list(formats)
    unknown = set(formats) - set(_ALL_FORMATS)
    if unknown:
        raise ValueError(f"Unknown format(s): {sorted(unknown)}; allowed: {_ALL_FORMATS}")

    # Sort a copy so the caller's model is untouched.
    model = sort_identifiers(model.copy())

    root = Path(path) / "model" if sub_dirs else Path(path)
    root.mkdir(parents=True, exist_ok=True)

    def target(fmt: str) -> Path:
        folder = root / fmt if sub_dirs else root
        folder.mkdir(parents=True, exist_ok=True)
        return folder / f"{prefix}.{fmt}"

    if "yml" in formats:
        write_yaml_model(model, target("yml"))
    if "xml" in formats:
        cobra.io.write_sbml_model(model, str(target("xml")))
    if "mat" in formats:
        cobra.io.save_matlab_model(model, str(target("mat")))
    if "xlsx" in formats:
        export_to_excel(model, target("xlsx"))
    if "txt" in formats:
        _write_txt(model, target("txt"))

    with open(root / "dependencies.txt", "w", encoding="utf-8") as fh:
        fh.write(f"python\t{platform.python_version()}\n")
        fh.write(f"cobra\t{_version('cobra')}\n")
        fh.write(f"ravengem\t{_version('ravengem')}\n")

    return root
