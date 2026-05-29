"""Fetch and cache published data artefacts (KEGG reference model, tables, HMMs).

The mirror of :mod:`raven_python.binaries` for *data*: a version-pinned registry of
downloadable artefacts, fetched on first use, SHA256-verified, and cached under
platformdirs so end users never rebuild them from a KEGG dump (that is the
maintainer's job — see docs/maintaining_kegg_data.md).

Resolution for any artefact file:

    explicit local dir  →  cached copy  →  download from the registry (verify,
        cache)  →  FileNotFoundError with guidance

The registry is **empty until the artefacts are published** (same as
``binaries._REGISTRY``); until then ``ensure_data_file`` raises an actionable
error. Cache layout::

    $XDG_CACHE_HOME/raven_python/data/<dataset>-<version>/<filename>
    (or ~/.cache/raven_python/data/... if XDG_CACHE_HOME is unset)
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path
from urllib.request import urlopen

from raven_python.binaries import _sha256

# dataset -> {"version": str, "files": {filename: {"url": str, "sha256": str}}}
# Populated when raven_python publishes the KEGG artefacts as release assets.
_DATA_REGISTRY: dict = {}

# The core KEGG artefacts needed to build a model (no HMM libraries).
CORE_KEGG_FILES = (
    "reference_model.yml.gz",
    "ko_reaction.tsv.gz",
    "ko_names.tsv.gz",
    "organism_gene_ko.tsv.xz",
    "rxn_flags.tsv.gz",
)


def _data_cache_dir() -> Path:
    base = os.environ.get("XDG_CACHE_HOME") or (Path.home() / ".cache")
    return Path(base) / "raven_python" / "data"


def _bundle(dataset: str, registry: dict) -> dict:
    bundle = registry.get(dataset)
    if bundle is None:
        raise FileNotFoundError(
            f"No data artefacts registered for {dataset!r}. Either pass a local "
            f"directory of artefacts, or build them per docs/maintaining_kegg_data.md."
        )
    return bundle


def ensure_data_file(
    dataset: str,
    filename: str,
    *,
    version: str | None = None,
    registry: dict | None = None,
) -> Path:
    """Download (if needed) and return the cached path to one artefact file.

    Looks the file up in the registry for ``dataset`` (at ``version`` or the
    registry's default), downloads it to the version-pinned cache directory,
    verifies its SHA256, and returns the path. Re-uses an already-cached copy.
    """
    registry = _DATA_REGISTRY if registry is None else registry
    bundle = _bundle(dataset, registry)
    ver = version or bundle["version"]
    entry = bundle.get("files", {}).get(filename)
    if entry is None:
        raise FileNotFoundError(
            f"{filename!r} is not registered for {dataset!r} {ver}. "
            f"Available: {sorted(bundle.get('files', {}))}."
        )

    dest_dir = _data_cache_dir() / f"{dataset}-{ver}"
    dest = dest_dir / filename
    if dest.exists():
        return dest

    dest_dir.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    with urlopen(entry["url"]) as resp, open(tmp, "wb") as out:  # noqa: S310 (trusted registry URLs)
        shutil.copyfileobj(resp, out)
    digest = _sha256(tmp)
    if digest != entry["sha256"]:
        tmp.unlink(missing_ok=True)
        raise ValueError(
            f"SHA256 mismatch for {dataset}/{filename} ({ver}): "
            f"expected {entry['sha256']}, got {digest}."
        )
    tmp.replace(dest)
    return dest


def ensure_kegg_data(
    *,
    version: str | None = None,
    files: tuple[str, ...] = CORE_KEGG_FILES,
    registry: dict | None = None,
) -> Path:
    """Ensure the core KEGG artefacts are cached; return their directory.

    Fetches each of ``files`` (default :data:`CORE_KEGG_FILES`) for the ``kegg``
    dataset and returns the cache directory holding them — ready to pass as the
    ``artefact_dir`` of :func:`get_kegg_model_for_organism_from_artefacts`.
    """
    registry = _DATA_REGISTRY if registry is None else registry
    ver = version or _bundle("kegg", registry)["version"]
    for filename in files:
        ensure_data_file("kegg", filename, version=ver, registry=registry)
    return _data_cache_dir() / f"kegg-{ver}"


def ensure_kegg_hmm_library(
    domain: str, *, version: str | None = None, registry: dict | None = None
) -> Path:
    """Ensure a domain HMM library (and its hmmpress index) is cached; return its path.

    ``domain`` is ``"prokaryotes"`` or ``"eukaryotes"``. Fetches ``<domain>.hmm``
    plus the ``hmmpress`` sidecar files (``.h3f/.h3i/.h3m/.h3p``) and returns the
    path to the ``.hmm`` (the argument for :func:`run_hmmscan`).
    """
    registry = _DATA_REGISTRY if registry is None else registry
    ver = version or _bundle("kegg", registry)["version"]
    base = f"{domain}.hmm"
    library = ensure_data_file("kegg", base, version=ver, registry=registry)
    for suffix in (".h3f", ".h3i", ".h3m", ".h3p"):
        ensure_data_file("kegg", base + suffix, version=ver, registry=registry)
    return library
