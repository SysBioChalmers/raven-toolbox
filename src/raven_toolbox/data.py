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

import gzip
import os
import shutil
import tarfile
from pathlib import Path
from urllib.request import urlopen

from raven_python.binaries import _sha256

# dataset -> {"version": str, "files": {filename: {"url": str, "sha256": str}}}
# Mirrors data/manifest.json (the cross-language source of truth); regenerate the
# block with scripts/make_registry_snippet.py when publishing a new KEGG release.
_DATA_REGISTRY: dict = {
    "kegg": {
        "version": "kegg116",
        "files": {
            "kegg116_core.tar.gz": {
                "url": "https://github.com/SysBioChalmers/raven-python/releases/download/v0.1.0/kegg116_core.tar.gz",
                "sha256": "155d5806d43db2fde5783fb124f8782bbcad390a1dd80879c520d2eac9d780e7",
            },
            "kegg116_eukaryotes.hmm.gz": {
                "url": "https://github.com/SysBioChalmers/raven-python/releases/download/v0.1.0/kegg116_eukaryotes.hmm.gz",
                "sha256": "2d48bc9935575d0f9ba4178bf2df19279bff866b49c1bf83a8e15787b11d6708",
            },
            "kegg116_prokaryotes.hmm.gz": {
                "url": "https://github.com/SysBioChalmers/raven-python/releases/download/v0.1.0/kegg116_prokaryotes.hmm.gz",
                "sha256": "d80cb2a22dec9fd8336b3998e3b96ee121672f63f4041cddaf09624fe739f1af",
            },
            "kegg116_taxonomy.gz": {
                "url": "https://github.com/SysBioChalmers/raven-python/releases/download/v0.1.0/kegg116_taxonomy.gz",
                "sha256": "1edc56da94d71433e5f08c133600292c311baaf33279a959518ab08389b0e538",
            },
        },
    },
}

# The core KEGG artefacts needed to build a model (no HMM libraries). These are
# the *base* names of the files bundled into the published ``<version>_core.tar.gz``
# (each stored version-prefixed inside the archive); ``ensure_kegg_data`` fetches the
# bundle and extracts these, and the build groups exactly this set.
CORE_KEGG_FILES = (
    "reference_model.yml.gz",
    "ko_reaction.tsv.gz",
    "ko_names.tsv.gz",
    "organism_gene_ko.tsv.gz",
    "rxn_flags.tsv.gz",
)


def _data_cache_dir() -> Path:
    base = os.environ.get("XDG_CACHE_HOME") or (Path.home() / ".cache")
    return Path(base) / "raven_python" / "data"


def _maybe_autoload(registry: dict) -> None:
    """Populate the default registry from ``$RAVEN_PYTHON_MANIFEST`` on first use, if set.

    Fires only when the caller relies on the default (still-empty) ``_DATA_REGISTRY`` and
    the environment variable points at a manifest. Local import avoids an import cycle with
    :mod:`raven_python.manifest`.
    """
    if registry is _DATA_REGISTRY and not registry and os.environ.get("RAVEN_PYTHON_MANIFEST"):
        from raven_python import manifest as _manifest

        _manifest.load_into_registries()


def _bundle(dataset: str, registry: dict) -> dict:
    _maybe_autoload(registry)
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
    verify: bool = False,
) -> Path:
    """Download (if needed) and return the cached path to one artefact file.

    Looks the file up in the registry for ``dataset`` (at ``version`` or the
    registry's default), downloads it to the version-pinned cache directory,
    verifies its SHA256, and returns the path. Re-uses an already-cached copy.

    A freshly downloaded file is always SHA256-checked. ``verify`` additionally
    re-checks an *already-cached* file's SHA256 (a mismatch — i.e. a corrupted
    cache — discards it and re-downloads); it is off by default so the common
    cache-hit path stays fast.
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
        if not verify or _sha256(dest) == entry["sha256"]:
            return dest
        dest.unlink()  # corrupted cache → fall through and re-download

    dest_dir.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    with urlopen(entry["url"], timeout=60) as resp, open(tmp, "wb") as out:  # noqa: S310 (trusted registry URLs)
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


def ensure_kegg_data(*, version: str | None = None, registry: dict | None = None) -> Path:
    """Ensure the core KEGG artefacts are cached; return their directory.

    Fetches the single ``<version>_core.tar.gz`` bundle (the gene-free reference
    model + the KO/reaction/organism-gene tables of :data:`CORE_KEGG_FILES`),
    SHA256-verifies it, and extracts the version-prefixed members into the cache
    directory on first use — ready to pass as the ``artefact_dir`` of
    :func:`get_kegg_model_for_organism_from_artefacts`. The HMM libraries and the
    taxonomy file are *separate* artefacts (see :func:`ensure_kegg_hmm_library`,
    :func:`ensure_kegg_taxonomy`).
    """
    registry = _DATA_REGISTRY if registry is None else registry
    ver = version or _bundle("kegg", registry)["version"]
    dest_dir = _data_cache_dir() / f"kegg-{ver}"
    archive = ensure_data_file("kegg", f"{ver}_core.tar.gz", version=ver, registry=registry)
    # Extract once; a marker avoids re-extracting (and re-reading the archive) per call.
    marker = dest_dir / ".core-extracted"
    if not marker.exists():
        with tarfile.open(archive, "r:gz") as tar:
            tar.extractall(dest_dir, filter="data")  # safe extraction (matches download.py)
        marker.touch()
    return dest_dir


def ensure_kegg_hmm_library(
    domain: str,
    *,
    version: str | None = None,
    registry: dict | None = None,
) -> Path:
    """Ensure a domain HMM library is cached and decompressed; return the ``.hmm`` path.

    ``domain`` is ``"prokaryotes"`` or ``"eukaryotes"``. Fetches the gzipped
    concatenated library ``<version>_<domain>.hmm.gz`` and decompresses it once
    (cached). Returns the path to the ``.hmm`` flatfile — the argument for
    :func:`run_hmmsearch`, which searches it directly (no ``hmmpress`` needed).

    Shipping the gzip flatfile keeps the download ~10x smaller than a binary index,
    stays portable across HMMER versions/platforms, and lets the same artefact serve
    MATLAB RAVEN.
    """
    registry = _DATA_REGISTRY if registry is None else registry
    ver = version or _bundle("kegg", registry)["version"]
    archive = ensure_data_file("kegg", f"{ver}_{domain}.hmm.gz", version=ver, registry=registry)
    library = archive.with_suffix("")  # strip ".gz" -> <version>_<domain>.hmm
    if not library.exists():
        tmp = library.with_name(library.name + ".part")
        with gzip.open(archive, "rb") as src, open(tmp, "wb") as out:
            shutil.copyfileobj(src, out)
        tmp.replace(library)
    return library


def ensure_kegg_taxonomy(*, version: str | None = None, registry: dict | None = None) -> Path:
    """Ensure the KEGG ``taxonomy`` artefact is cached; return its (gzipped) path.

    The gzipped KEGG ``taxonomy`` file is the source for domain classification and for
    regenerating the phylogenetic distance matrix — RAVEN's ``keggPhylDist``, which GECKO
    uses to pick the closest organism for kcat assignment — via
    :func:`raven_python.reconstruction.kegg.phyl_dist` (which reads ``.gz`` directly). So
    that capability needs only this published artefact, no MATLAB ``.mat`` file.
    """
    registry = _DATA_REGISTRY if registry is None else registry
    ver = version or _bundle("kegg", registry)["version"]
    return ensure_data_file("kegg", f"{ver}_taxonomy.gz", version=ver, registry=registry)
