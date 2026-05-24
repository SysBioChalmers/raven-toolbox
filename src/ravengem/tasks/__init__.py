"""Metabolic task definition, parsing, and checking (Phase 4a).

* :class:`Task` + :func:`parse_task_list` — the task-list file format (``parseTaskList``).
* :func:`check_tasks` + :class:`TaskResult` — run tasks against a model (``checkTasks``).
"""
from ravengem.tasks.check import TaskResult, check_tasks
from ravengem.tasks.tasklist import Task, parse_task_list

__all__ = ["Task", "TaskResult", "check_tasks", "parse_task_list"]
