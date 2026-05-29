"""Locate and provision external command-line binaries (BLAST+, DIAMOND, …).

Shared across tools (not homology-specific). Resolution order for any executable:

    explicit path arg  →  env var (RAVENGEM_<TOOL>)  →  shutil.which (PATH)
      →  ensure_binary  (download the version-pinned ZIP from a ravengem release,
                         verify SHA256, cache, return the path)
      →  FileNotFoundError with install guidance

So a pre-installed/conda binary always wins; the bundled ZIP is the zero-setup
fallback. See docs/maintaining_binaries.md for how the release ZIPs and the
registry are produced and updated.
"""
from __future__ import annotations

import hashlib
import os
import platform
import shutil
import zipfile
from pathlib import Path
from urllib.request import urlopen

# Registry of bundled binaries. Empty until release ZIPs are published; populated
# per docs/maintaining_binaries.md. Keyed by *bundle*; one bundle can provide
# several executables (e.g. "blast" -> blastp + makeblastdb).
#   bundle -> {version, provides:[exe...], platforms:{"<os>-<arch>": {url, sha256}}}
_REGISTRY: dict = {}

# Environment variable overrides per executable.
_ENV_VARS = {
    "diamond": "RAVENGEM_DIAMOND",
    "blastp": "RAVENGEM_BLASTP",
    "makeblastdb": "RAVENGEM_MAKEBLASTDB",
    "hmmbuild": "RAVENGEM_HMMBUILD",
    "hmmpress": "RAVENGEM_HMMPRESS",
    "hmmsearch": "RAVENGEM_HMMSEARCH",
    "hmmscan": "RAVENGEM_HMMSCAN",
    "mafft": "RAVENGEM_MAFFT",
    "cd-hit": "RAVENGEM_CDHIT",
}


def platform_key() -> str:
    """Return the ``<os>-<arch>`` key used in the registry (e.g. ``linux-x86_64``)."""
    system = {"linux": "linux", "darwin": "macos", "windows": "windows"}.get(
        platform.system().lower(), platform.system().lower()
    )
    machine = platform.machine().lower()
    arch = {"x86_64": "x86_64", "amd64": "x86_64", "arm64": "arm64", "aarch64": "arm64"}.get(
        machine, machine
    )
    return f"{system}-{arch}"


def _cache_dir() -> Path:
    base = os.environ.get("XDG_CACHE_HOME") or (Path.home() / ".cache")
    return Path(base) / "ravengem" / "binaries"


def _bundle_for(executable: str, registry: dict):
    for name, bundle in registry.items():
        if executable in bundle.get("provides", []):
            return name, bundle
    return None, None


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def ensure_binary(executable: str, *, registry: dict | None = None) -> Path:
    """Download (if needed) and return the path to a bundled ``executable``.

    Consults the registry for the current platform, downloads the pinned ZIP,
    verifies its SHA256, extracts it into the cache, and returns the executable
    path. Raises ``FileNotFoundError`` if no bundle for this platform is hosted.
    """
    registry = _REGISTRY if registry is None else registry
    bundle_name, bundle = _bundle_for(executable, registry)
    if bundle is None:
        raise FileNotFoundError(
            f"No bundled binary registered for {executable!r}. Install it (e.g. "
            f"`conda install -c bioconda {executable}`) or pass an explicit path."
        )
    key = platform_key()
    entry = bundle.get("platforms", {}).get(key)
    if entry is None:
        raise FileNotFoundError(
            f"No bundled {executable!r} for platform {key!r}. Install it "
            f"(e.g. `conda install -c bioconda {executable}`), set "
            f"{_ENV_VARS.get(executable, 'the binary path')}, or pass binary=."
        )

    dest_dir = _cache_dir() / f"{bundle_name}-{bundle['version']}-{key}"
    exe = dest_dir / executable
    if exe.exists():
        return exe

    dest_dir.mkdir(parents=True, exist_ok=True)
    archive = dest_dir / "_download.zip"
    # Download into a sibling .part file and rename on success — an interrupted
    # download leaves the partial behind .part, never as a half-complete .zip
    # that a later run might mistake for a finished one. Mirrors data.py.
    part = archive.with_suffix(archive.suffix + ".part")
    try:
        with urlopen(entry["url"]) as resp, open(part, "wb") as out:  # noqa: S310
            shutil.copyfileobj(resp, out)
        digest = _sha256(part)
        if digest != entry["sha256"]:
            raise ValueError(
                f"SHA256 mismatch for {executable!r} ({key}): "
                f"expected {entry['sha256']}, got {digest}."
            )
        os.replace(part, archive)
    finally:
        part.unlink(missing_ok=True)
    with zipfile.ZipFile(archive) as zf:
        zf.extractall(dest_dir)
    archive.unlink(missing_ok=True)
    if not exe.exists():
        raise FileNotFoundError(f"{executable!r} not found in the extracted bundle at {dest_dir}.")
    exe.chmod(0o755)
    return exe


def resolve_binary(executable: str, *, binary: str | os.PathLike | None = None) -> str:
    """Resolve an executable to a path: arg → env var → PATH → bundled ZIP → error."""
    if binary is not None:
        return os.fspath(binary)
    env_var = _ENV_VARS.get(executable)
    if env_var and os.environ.get(env_var):
        return os.environ[env_var]
    found = shutil.which(executable)
    if found:
        return found
    try:
        return os.fspath(ensure_binary(executable))
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            f"Could not find {executable!r}. Install it (e.g. "
            f"`conda install -c bioconda {executable}`), put it on PATH, set "
            f"{env_var or 'the binary path'}, or pass binary=. ({exc})"
        ) from exc
