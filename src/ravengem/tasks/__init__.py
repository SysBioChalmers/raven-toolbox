"""Metabolic task definition, parsing, and checking (Phase 4a).

* :class:`Task` + :func:`parse_task_list` — the task-list file format (``parseTaskList``).
* :func:`check_tasks` + :class:`TaskResult` — run tasks against a model (``checkTasks``).
* :func:`find_task_essential_reactions` + :class:`EssentialReactionsResult` — reactions
  a model must use to satisfy a task list (the ``prepINITModel`` input for (ft)INIT).
"""
from ravengem.tasks.check import (
    EssentialReactionsResult,
    TaskResult,
    check_tasks,
    find_task_essential_reactions,
)
from ravengem.tasks.tasklist import Task, parse_task_list

__all__ = [
    "EssentialReactionsResult",
    "Task",
    "TaskResult",
    "check_tasks",
    "find_task_essential_reactions",
    "parse_task_list",
]
