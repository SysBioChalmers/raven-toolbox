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
import subprocess
from pathlib import Path
from urllib.request import urlopen

from raven_python.binaries import _sha256, resolve_binary

# dataset -> {"version": str, "files": {filename: {"url": str, "sha256": str}}}
# Mirrors data/manifest.json (the cross-language source of truth); regenerate the
# block with scripts/make_registry_snippet.py when publishing a new KEGG release.
_DATA_REGISTRY: dict = {
    "kegg": {
        "version": "kegg116",
        "files": {
            "kegg116_eukaryotes.hmm.gz": {
                "url": "https://github.com/SysBioChalmers/raven-python/releases/download/v0.1.0/kegg116_eukaryotes.hmm.gz",
                "sha256": "2d48bc9935575d0f9ba4178bf2df19279bff866b49c1bf83a8e15787b11d6708",
            },
            "kegg116_ko_names.tsv.gz": {
                "url": "https://github.com/SysBioChalmers/raven-python/releases/download/v0.1.0/kegg116_ko_names.tsv.gz",
                "sha256": "84f9c7150172d948f794d91a6608d55f7140f31e53249c705057ae49b11c93b3",
            },
            "kegg116_ko_reaction.tsv.gz": {
                "url": "https://github.com/SysBioChalmers/raven-python/releases/download/v0.1.0/kegg116_ko_reaction.tsv.gz",
                "sha256": "e1a4ac22875bd3030d03b78368b0153b6d99000acb2ee0f474340a03c180323c",
            },
            "kegg116_organism_gene_ko.tsv.gz": {
                "url": "https://github.com/SysBioChalmers/raven-python/releases/download/v0.1.0/kegg116_organism_gene_ko.tsv.gz",
                "sha256": "27bf7dd58eb1acd5904990dc2be187aae4d8d9b9f7421375618e7c8d6ff7253d",
            },
            "kegg116_prokaryotes.hmm.gz": {
                "url": "https://github.com/SysBioChalmers/raven-python/releases/download/v0.1.0/kegg116_prokaryotes.hmm.gz",
                "sha256": "d80cb2a22dec9fd8336b3998e3b96ee121672f63f4041cddaf09624fe739f1af",
            },
            "kegg116_reference_model.yml.gz": {
                "url": "https://github.com/SysBioChalmers/raven-python/releases/download/v0.1.0/kegg116_reference_model.yml.gz",
                "sha256": "73ff313fe2aa2830ec511f4e522226c98c5714c2d5c4632844544e5a409c7f0c",
            },
            "kegg116_rxn_flags.tsv.gz": {
                "url": "https://github.com/SysBioChalmers/raven-python/releases/download/v0.1.0/kegg116_rxn_flags.tsv.gz",
                "sha256": "c4c134effc9edeeb74b925ae8616320af162edbaad3a9b44dcc29d2c4d12db9b",
            },
        },
    },
}

# The core KEGG artefacts needed to build a model (no HMM libraries). These are
# the *base* names; published assets are version-prefixed (``<version>_<base>``),
# which is what the resolvers below construct and what the registry keys hold.
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


def ensure_kegg_data(
    *,
    version: str | None = None,
    files: tuple[str, ...] = CORE_KEGG_FILES,
    registry: dict | None = None,
) -> Path:
    """Ensure the core KEGG artefacts are cached; return their directory.

    Fetches each of ``files`` (default :data:`CORE_KEGG_FILES`, given as *base*
    names) for the ``kegg`` dataset and returns the cache directory holding them —
    ready to pass as the ``artefact_dir`` of
    :func:`get_kegg_model_for_organism_from_artefacts`. Each file is fetched under
    its version-prefixed published name (``<version>_<base>``).
    """
    registry = _DATA_REGISTRY if registry is None else registry
    ver = version or _bundle("kegg", registry)["version"]
    for base in files:
        ensure_data_file("kegg", f"{ver}_{base}", version=ver, registry=registry)
    return _data_cache_dir() / f"kegg-{ver}"


def ensure_kegg_hmm_library(
    domain: str,
    *,
    version: str | None = None,
    registry: dict | None = None,
    hmmpress: str | os.PathLike | None = None,
) -> Path:
    """Ensure a domain HMM library is cached and pressed; return the ``.hmm`` path.

    ``domain`` is ``"prokaryotes"`` or ``"eukaryotes"``. Fetches the gzipped
    concatenated library ``<version>_<domain>.hmm.gz``, decompresses it once, and
    runs ``hmmpress`` to build the ``.h3f/.h3i/.h3m/.h3p`` index ``hmmscan`` needs
    (HMMER is already a requirement of the de-novo query path). Both steps are
    cached, so they run only on first use. Returns the path to the decompressed
    ``.hmm`` (the argument for :func:`run_hmmscan`).

    Shipping the gzip flatfile and pressing on the client keeps the download ~10x
    smaller than the binary index, stays portable across HMMER versions/platforms,
    and lets the same artefact serve MATLAB RAVEN.
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
    sidecars = (".h3f", ".h3i", ".h3m", ".h3p")
    if not all(library.with_name(library.name + s).exists() for s in sidecars):
        exe = resolve_binary("hmmpress", binary=hmmpress)
        proc = subprocess.run([exe, "-f", str(library)], capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"hmmpress failed:\n{(proc.stderr or '').strip()}")
    return library
