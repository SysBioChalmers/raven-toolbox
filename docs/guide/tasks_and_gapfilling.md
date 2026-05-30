# Metabolic tasks & gap-filling

## Tasks — {mod}`raven_python.tasks`

A metabolic task asserts that a model can (or cannot) produce/consume a set of metabolites
under given constraints.

- {class}`raven_python.tasks.Task` + {func}`raven_python.tasks.parse_task_list` — the
  task-list file format.
- {func}`raven_python.tasks.check_tasks` + {class}`raven_python.tasks.TaskResult` — run tasks
  against a model. `check_tasks` reuses one model across the whole list (no per-task model
  copy) — at genome scale ~12× faster than a copy-per-task implementation.
- {func}`raven_python.tasks.find_task_essential_reactions` — the reactions a model must use to
  satisfy a task list. This is the input to (f)tINIT's task layer (see the
  [context-specific modeling guide](context_specific.md)).

## Connectivity gap-filling — {mod}`raven_python.gapfilling`

{func}`raven_python.gapfilling.connect_blocked_reactions` adds the fewest (lowest-penalty)
template reactions so reactions that are blocked in a draft can carry flux — the connectivity
flavour of RAVEN's `fillGaps`.

:::{tip}
For the other flavour — add reactions until a specific objective becomes *feasible* — use
cobra's {func}`cobra.flux_analysis.gapfill` directly.
:::
