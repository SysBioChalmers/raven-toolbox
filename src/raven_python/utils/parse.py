"""Small parsing helpers shared across raven_python."""
from __future__ import annotations

import re

# A metabolite written as ``name[comp]``. The name is greedy so that, for a
# pathological name that itself contains brackets, the *last* ``[...]`` is taken
# as the compartment (matching RAVEN getIndexes' ``max(strfind('['))`` rule).
_NAME_COMP_RE = re.compile(r"^(?P<name>.+)\[(?P<comp>[^\[\]]+)\]$")


def parse_name_comp(token: str) -> tuple[str, str | None]:
    """Split a ``name[comp]`` token into ``(name, compartment)``.

    This is the one genuinely cobra-absent sliver of RAVEN ``getIndexes``'
    ``metcomps`` mode and ``addRxns`` eqnType 3: resolving a metabolite written
    as its *name* plus a compartment in square brackets, e.g. ``"ATP[c]"``.

    Returns ``(name, None)`` when there is no trailing ``[...]``.

    Examples
    --------
    >>> parse_name_comp("ATP[c]")
    ('ATP', 'c')
    >>> parse_name_comp("ATP")
    ('ATP', None)
    >>> parse_name_comp("weird[name][m]")
    ('weird[name]', 'm')
    """
    match = _NAME_COMP_RE.match(token.strip())
    if match:
        return match.group("name").strip(), match.group("comp").strip()
    return token.strip(), None


def subsystem_to_str(value) -> str:
    """Coerce a reaction ``subsystem`` to cobra's canonical ``str`` form.

    cobra defines ``Reaction.subsystem`` as a plain string, but RAVEN/MATLAB and
    "messy" YAML sometimes deliver a list/tuple of subsystem names (a reaction can
    belong to several). This collapses any such value to a single ``;``-joined
    string of its parts — losing no data, unlike taking only the first — and returns
    ``""`` for an empty/None subsystem. Used wherever a subsystem is rendered or
    compared so the handling is identical across modules.

    Examples
    --------
    >>> subsystem_to_str("glycolysis")
    'glycolysis'
    >>> subsystem_to_str(["glycolysis", "TCA cycle"])
    'glycolysis;TCA cycle'
    >>> subsystem_to_str(None)
    ''
    """
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return ";".join(str(s) for s in value if str(s).strip())
    return str(value)
