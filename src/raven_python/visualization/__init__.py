"""Pathway-map and omics-overlay visualisation (stub — not yet implemented)."""
from __future__ import annotations


def __getattr__(name: str):
    # PEP 562: accessing any attribute of this stub package fails with a clear,
    # actionable message rather than a bare AttributeError. Dunder lookups
    # (``__path__`` etc.) fall through to normal handling.
    if name.startswith("__"):
        raise AttributeError(name)
    raise NotImplementedError(
        f"raven_python.visualization.{name!r} is not implemented yet. Pathway-map "
        "and omics-overlay visualisation is on the roadmap; until then export the "
        "model via raven_python.io and use an external viewer such as Escher."
    )
