"""ftINIT preprocessing — reaction classification + prepData (port of ``prepINITModel``).

ftINIT does all omics-independent work once per template model: classify reactions
into the categories the staged MILP may *ignore* (leave in, never remove), discover
task-essential reactions, linearly merge, and scale. The result (:class:`PrepData`)
is reused across every sample.

:func:`classify_reactions` is the reaction taxonomy (RAVEN's ``toIgnore*`` masks):
exchange, GPR-less import/simple/advanced transport, spontaneous, GPR-less
extracellular, custom, and "any without a GPR". The staged schedule
(:func:`ravengem.init.get_init_steps`) selects which categories to keep out of each
MILP step via an 8-bit pattern.
"""
from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field

import cobra

from ravengem.init.merge import merge_linear
from ravengem.tasks import Task, find_task_essential_reactions


@dataclass
class ReactionMasks:
    """Reaction-category id sets (RAVEN's ``toIgnore*``), in 8-bit-pattern order.

    ``ignored(pattern)`` returns the union of the categories whose bit is set — the
    reactions held out of (left untouched by) that MILP step.
    """

    exchange: set[str] = field(default_factory=set)            # b1
    import_rxns: set[str] = field(default_factory=set)         # b2
    simple_transport: set[str] = field(default_factory=set)    # b3
    advanced_transport: set[str] = field(default_factory=set)  # b4
    spontaneous: set[str] = field(default_factory=set)         # b5
    extracellular: set[str] = field(default_factory=set)       # b6 (no-GPR, all mets in ext comp)
    custom: set[str] = field(default_factory=set)              # b7
    no_gpr: set[str] = field(default_factory=set)              # b8

    def _ordered(self) -> list[set[str]]:
        return [self.exchange, self.import_rxns, self.simple_transport,
                self.advanced_transport, self.spontaneous, self.extracellular,
                self.custom, self.no_gpr]

    def ignored(self, pattern: Iterable[int]) -> set[str]:
        out: set[str] = set()
        for bit, group in zip(pattern, self._ordered(), strict=True):
            if bit:
                out |= group
        return out


def _is_advanced_transport(rxn: cobra.Reaction) -> bool:
    """Even number (>2) of mets pairing up by name across compartments with canceling stoich."""
    mets = list(rxn.metabolites.items())
    if len(mets) <= 2 or len(mets) % 2 != 0:
        return False
    remaining = [(m.name, m.compartment, c) for m, c in mets]
    while remaining:
        name, comp, coeff = remaining[0]
        matches = [i for i in range(1, len(remaining)) if remaining[i][0] == name]
        if len(matches) != 1:
            return False
        j = matches[0]
        if coeff + remaining[j][2] != 0 or comp == remaining[j][1]:
            return False
        remaining = [r for k, r in enumerate(remaining) if k not in (0, j)]
    return True


def classify_reactions(
    model: cobra.Model,
    *,
    ext_comp: str = "e",
    spontaneous: Iterable[str] = (),
    custom: Iterable[str] = (),
) -> ReactionMasks:
    """Classify reactions into the ftINIT ``toIgnore`` categories (``prepINITModel``).

    ``ext_comp`` is the extracellular compartment. ``spontaneous``/``custom`` are
    reaction-id lists. A reaction is "GPR-less" when its gene rule is empty.
    """
    spont, cust = set(spontaneous), set(custom)
    masks = ReactionMasks(
        exchange={r.id for r in model.boundary},
        spontaneous={r.id for r in model.reactions if r.id in spont},
        custom={r.id for r in model.reactions if r.id in cust},
        no_gpr={r.id for r in model.reactions if not r.gene_reaction_rule.strip()},
    )
    for rxn in model.reactions:
        if rxn.gene_reaction_rule.strip():
            continue  # transport categories are GPR-less only
        mets = list(rxn.metabolites)
        if len(mets) == 2:
            (m1, m2) = mets
            if m1.compartment != m2.compartment and m1.name == m2.name:
                if ext_comp in (m1.compartment, m2.compartment):
                    masks.import_rxns.add(rxn.id)
                else:
                    masks.simple_transport.add(rxn.id)
        elif _is_advanced_transport(rxn):
            masks.advanced_transport.add(rxn.id)
        if len(mets) > 1 and all(m.compartment == ext_comp for m in mets):
            masks.extracellular.add(rxn.id)
    return masks


@dataclass
class PrepData:
    """One-time ftINIT preprocessing of a template model (RAVEN ``prepData``).

    Built once per template, reused across samples. ``min_model`` is the merged model
    the MILP runs on; ``orig_rxn_ids``/``group_ids`` map its reactions back to the
    ``ref_model`` (the simplified, pre-merge reference). ``essential_rxns`` are in
    **merged** ids and pre-oriented irreversibly (so the MILP forces flux *forward*).
    ``masks`` is on ``ref_model`` (= original) ids.
    """

    ref_model: cobra.Model
    min_model: cobra.Model
    orig_rxn_ids: list[str]
    group_ids: list[int]
    reversed_rxns: list[bool]
    masks: ReactionMasks
    essential_rxns: set[str] = field(default_factory=set)
    essential_mets_for_tasks: set[str] = field(default_factory=set)
    tasks: list[Task] = field(default_factory=list)

    @property
    def group_of(self) -> dict[str, int]:
        return dict(zip(self.orig_rxn_ids, self.group_ids, strict=True))


def _orient_forward(rxn: cobra.Reaction, direction: int) -> None:
    """Make ``rxn`` carry flux only in its forced direction (irreversible forward)."""
    if direction < 0:  # flip so the forced (reverse) direction becomes forward
        rxn.add_metabolites({m: -2 * c for m, c in rxn.metabolites.items()})
        rxn.bounds = (-rxn.upper_bound, -rxn.lower_bound)
    rxn.lower_bound = max(rxn.lower_bound, 0.0)


def prep_init_model(
    template: cobra.Model,
    tasks: Iterable[Task] | None = None,
    *,
    ext_comp: str = "e",
    spontaneous: Iterable[str] = (),
    custom: Iterable[str] = (),
) -> PrepData:
    """Build :class:`PrepData` from a template model (port of ``prepINITModel``).

    With ``tasks``, discovers the task-essential reactions (kept regardless of score),
    orients them irreversibly in their required direction, and drops tasks that are
    infeasible. Then classifies reactions into the ``toIgnore`` categories and linearly
    merges. (Dead-end simplification and numerical scaling are deferred to 4d.7; on
    genome-scale templates they should be applied here.)
    """
    ref_model = template.copy()

    essential_pre: dict[str, int] = {}
    task_mets: set[str] = set()
    kept_tasks: list[Task] = []
    if tasks is not None:
        tasks = list(tasks)
        ess = find_task_essential_reactions(ref_model, tasks)
        essential_pre = ess.reactions
        task_mets = ess.task_metabolites
        kept_tasks = [t for t in tasks if t.id not in ess.failed_tasks]

    # Orient essentials irreversibly (forced direction → forward) before merging, so
    # the merge keeps them forward and the MILP forces them with a simple lower bound.
    for rid, direction in essential_pre.items():
        _orient_forward(ref_model.reactions.get_by_id(rid), direction)

    masks = classify_reactions(ref_model, ext_comp=ext_comp,
                               spontaneous=spontaneous, custom=custom)

    min_model, orig_ids, group_ids, reversed_rxns = merge_linear(ref_model)
    group_of = dict(zip(orig_ids, group_ids, strict=True))

    # Map essentials to the merged model: the survivor of each group containing an
    # essential (or the reaction itself if unmerged). All are forward after orientation.
    # An essential that merged into a group which collapsed away (e.g. a trivial
    # source→sink chain) has no survivor and imposes no constraint — skip it.
    survivor_by_group = {group_of[r.id]: r.id for r in min_model.reactions if group_of[r.id]}
    essential_merged: set[str] = set()
    for rid in essential_pre:
        gid = group_of[rid]
        if gid == 0:
            essential_merged.add(rid)
        elif gid in survivor_by_group:
            essential_merged.add(survivor_by_group[gid])

    return PrepData(
        ref_model=ref_model,
        min_model=min_model,
        orig_rxn_ids=orig_ids,
        group_ids=group_ids,
        reversed_rxns=reversed_rxns,
        masks=masks,
        essential_rxns=essential_merged,
        essential_mets_for_tasks=task_mets,
        tasks=kept_tasks,
    )
