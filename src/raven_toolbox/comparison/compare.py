"""N-model structural and functional comparison.

Compare two or more models — typically context-specific models extracted from the same
template — on their reactions, metabolites, genes, subsystems, and (optionally) which
metabolic tasks they perform. Returns tidy :class:`pandas.DataFrame`\\ s suitable for
downstream plotting (heatmaps, tSNE/MDS, …) in seaborn / scikit-learn; plotting is
intentionally not in this function so it stays usable inside pipelines.

All matrices use the union of ids across the input models as the row index, so missing
entries are unambiguously ``0`` / ``False`` rather than ``NaN``.
"""
from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import cast

import cobra
import pandas as pd

from raven_toolbox.tasks import Task, check_tasks
from raven_toolbox.utils.parse import subsystem_to_str


@dataclass
class ModelComparison:
    """Tabular result of :func:`compare_models`.

    All matrices are indexed by id (reactions/metabolites/genes/subsystems) with one
    column per model. ``presence`` matrices are 0/1; ``subsystems`` is the per-model
    reaction count per subsystem. ``similarity`` is the model × model Jaccard on the
    reaction set (1 = identical, 0 = disjoint).
    """

    model_ids: list[str]
    reactions: pd.DataFrame
    metabolites: pd.DataFrame
    genes: pd.DataFrame
    subsystems: pd.DataFrame
    similarity: pd.DataFrame
    tasks: pd.DataFrame | None = None  # filled iff tasks were supplied
    failed_tasks: dict[str, list[str]] = field(default_factory=dict)


def _presence_matrix(items_per_model: list[list[str]], model_ids: list[str]) -> pd.DataFrame:
    """Build a 0/1 DataFrame: union of items as index × one column per model."""
    ordered: list[str] = []
    seen: set[str] = set()
    for items in items_per_model:
        for it in items:
            if it not in seen:
                seen.add(it)
                ordered.append(it)
    df = pd.DataFrame(0, index=ordered, columns=model_ids, dtype="int8")
    for mid, items in zip(model_ids, items_per_model, strict=True):
        if items:  # avoid empty-list edge case
            df.loc[list(set(items) & seen), mid] = 1
    return df


def _subsystem_counts(model: cobra.Model) -> dict[str, int]:
    """{subsystem_name: reaction_count}. Reactions with empty subsystem fall under '(none)'."""
    counts: dict[str, int] = {}
    for r in model.reactions:
        # cobra stores subsystem as a string; RAVEN sometimes uses cell-of-cells.
        # Coerce to the same ;-joined string used everywhere else (no data lost).
        sub = subsystem_to_str(r.subsystem).strip() or "(none)"
        counts[sub] = counts.get(sub, 0) + 1
    return counts


def _jaccard_matrix(presence: pd.DataFrame) -> pd.DataFrame:
    """Pairwise Jaccard similarity from a 0/1 presence matrix (rows = items, cols = models)."""
    arr = presence.values.astype(bool)
    out = pd.DataFrame(0.0, index=presence.columns, columns=presence.columns)
    for i, a in enumerate(presence.columns):
        ai = arr[:, i]
        for j, b in enumerate(presence.columns):
            bj = arr[:, j]
            inter = int((ai & bj).sum())
            union = int((ai | bj).sum())
            out.loc[a, b] = inter / union if union else 1.0
    return out


def compare_models(
    models: Iterable[cobra.Model],
    *,
    tasks: str | Iterable[Task] | None = None,
) -> ModelComparison:
    """Compare N cobra models on their reactions / metabolites / genes / subsystems
    (and tasks, if provided).

    ``tasks`` is forwarded to :func:`raven_toolbox.tasks.check_tasks` on each model; pass a
    file path or a parsed task list. When omitted, ``ModelComparison.tasks`` is ``None``.

    Models are identified by ``model.id`` (with a fallback to ``model_<i>`` if missing
    or duplicated).
    """
    models_list = list(models)
    if len(models_list) < 2:
        raise ValueError(f"compare_models needs ≥2 models; got {len(models_list)}")

    # Unique, stable model ids.
    model_ids: list[str] = []
    seen: set[str] = set()
    for i, m in enumerate(models_list):
        mid = (m.id or "").strip() or f"model_{i}"
        base, n = mid, 2
        while mid in seen:
            mid, n = f"{base}__{n}", n + 1
        seen.add(mid)
        model_ids.append(mid)

    reactions = _presence_matrix([[r.id for r in m.reactions] for m in models_list], model_ids)
    metabolites = _presence_matrix([[x.id for x in m.metabolites] for m in models_list], model_ids)
    genes = _presence_matrix([[g.id for g in m.genes] for m in models_list], model_ids)

    # Subsystems: union of names, per-model reaction counts.
    sub_counts = [_subsystem_counts(m) for m in models_list]
    sub_ids = sorted({s for c in sub_counts for s in c})
    subsystems = pd.DataFrame(0, index=sub_ids, columns=model_ids, dtype="int32")
    for mid, c in zip(model_ids, sub_counts, strict=True):
        for s, n in c.items():
            subsystems.at[s, mid] = n

    similarity = _jaccard_matrix(reactions)

    task_df: pd.DataFrame | None = None
    failed: dict[str, list[str]] = {}
    if tasks is not None:
        # raven_toolbox.tasks.check_tasks accepts a path or an iterable of Task; preserve task
        # ids for the index. Capture the list once so all models test the same set.
        from raven_toolbox.tasks.tasklist import parse_task_list
        task_list = (parse_task_list(cast("str | Path", tasks))
                     if isinstance(tasks, (str, bytes)) or hasattr(tasks, "__fspath__")
                     else list(tasks))
        task_ids = [t.id for t in task_list]
        task_df = pd.DataFrame(False, index=task_ids, columns=model_ids, dtype=bool)
        for mid, m in zip(model_ids, models_list, strict=True):
            results = check_tasks(m, task_list)
            for r in results:
                task_df.at[r.id, mid] = bool(r.passed)
                if not r.passed and r.error:
                    failed.setdefault(mid, []).append(f"{r.id}: {r.error}")

    return ModelComparison(model_ids=model_ids, reactions=reactions, metabolites=metabolites,
                           genes=genes, subsystems=subsystems, similarity=similarity,
                           tasks=task_df, failed_tasks=failed)
