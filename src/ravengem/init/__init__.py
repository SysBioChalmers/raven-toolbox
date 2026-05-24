"""Context-specific model extraction (tINIT / ftINIT).

Phase 4c (tINIT): :func:`run_init` — the INIT MILP (port of RAVEN ``runINIT``).
The expression-scoring wrapper (``getINITModel``) and ftINIT (Phase 4d) follow.
"""
from ravengem.init.init import InitResult, run_init

__all__ = ["InitResult", "run_init"]
