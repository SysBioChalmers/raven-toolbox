"""ftINIT step schedule.

ftINIT runs as a short sequence of MILP steps instead of one big MILP. Each step
(:class:`InitStep`) chooses which reaction categories to hold out of the problem
(``ignore_mask``, an 8-bit pattern over :class:`raven_toolbox.init.ReactionMasks`), whether
to drop positive reversibles and allow metabolite secretion, and how to treat the
reactions turned on by previous steps (``'ignore'`` for the first step, ``'essential'``
to fix them on). :func:`get_init_steps` builds the standard schedules.

The default ``'1+1'`` is two steps: step 1 decides only the GPR-associated reactions
(everything GPR-less is held out); step 2 brings the GPR-less transport / extracellular
reactions in with step-1 reactions fixed as essential. ``'full'`` is the single-MILP
classic-tINIT variant (nothing held out).
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

# 8-bit ignore patterns (exchange, import, simple-transp, adv-transp, spontaneous,
# extracellular, custom, no-GPR) — see ReactionMasks.
_ALL_NO_GPR_KEPT = (1, 1, 1, 1, 1, 1, 1, 0)  # hold out every GPR-less category but "all no-GPR"
_EXCH_SPONT = (1, 0, 0, 0, 1, 0, 0, 0)        # hold out only exchange + spontaneous
_NONE = (0, 0, 0, 0, 0, 0, 0, 0)


@dataclass
class InitStep:
    """One ftINIT MILP step."""

    how_to_use_prev: str = "essential"          # 'ignore' | 'essential'
    ignore_mask: tuple[int, ...] = _ALL_NO_GPR_KEPT
    pos_rev_off: bool = False                    # drop positive reversibles from the problem
    allow_met_secr: bool = False                 # relax S·v = 0 to ≥ 0
    mets_to_ignore: Sequence[str] = field(default_factory=tuple)  # met names zeroed from S (e.g. H2O)


def get_init_steps(series: str = "1+1", *, mets_to_ignore: Sequence[str] = ()) -> list[InitStep]:
    """Return the step schedule for a named ftINIT ``series`` (RAVEN ``getINITSteps``).

    ``'1+1'`` (default, step 1+2 merged), ``'2+1'`` (3-step), ``'1+0'``/``'2+0'``
    (skip the final GPR-less step), ``'full'`` (single MILP). ``mets_to_ignore`` are
    metabolite names removed from the stoichiometry in each step (e.g. H2O, H+).
    """
    m = tuple(mets_to_ignore)
    s1 = InitStep("ignore", _ALL_NO_GPR_KEPT, mets_to_ignore=m)
    s1_posrev = InitStep("ignore", _ALL_NO_GPR_KEPT, pos_rev_off=True, allow_met_secr=True,
                         mets_to_ignore=m)
    s2_all = InitStep("essential", _ALL_NO_GPR_KEPT, mets_to_ignore=m)
    s_final = InitStep("essential", _EXCH_SPONT, mets_to_ignore=m)

    if series == "1+1":
        return [s1, s_final]
    if series == "2+1":
        return [s1_posrev, s2_all, s_final]
    if series == "1+0":
        return [s1]
    if series == "2+0":
        return [s1_posrev, s2_all]
    if series == "full":
        return [InitStep("ignore", _NONE, mets_to_ignore=m)]
    raise ValueError(f"Unknown ftINIT series {series!r}; expected 1+1, 2+1, 1+0, 2+0, full.")
