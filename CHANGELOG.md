# Changelog

Milestones in the ravengem port. For function-level status see
[docs/raven_migration.md](docs/raven_migration.md); for open work see
[docs/todo.md](docs/todo.md).

## Infrastructure

* **GitHub Actions CI** ([.github/workflows/ci.yml](.github/workflows/ci.yml)) —
  ruff + pytest matrix over Python 3.11/3.12/3.13. Tests that require Gurobi
  auto-skip (no Gurobi on free runners); the known HiGHS upstream blocker
  (`hybrid_interface.Configuration` rejects `lp_method='primal'`) is marked
  `xfail(strict=True)` so CI flips red when optlang fixes it.

## Quality sweep — known-issues sections C / D / E

Closed all the robustness, efficiency, and dead-code items in one pass.

**Robustness (C):**
* `constrain_reversible_reactions` wraps FVA in try/except + NaN check; both
  backend-raised `OptimizationError` and silent-NaN returns now surface as one
  clear `RuntimeError` (the original `abs(NaN) < eps` silently no-op'd).
* `ensure_binary` downloads through `.part` + `os.replace`, matching `data.py` —
  an interrupted download leaves a `.part`, never a half-complete `.zip`.
* `parse_task_list` (.xlsx) checks `wb.sheetnames` before lookup; missing
  `TASKS` sheet now raises a clear `ValueError` instead of a bare `KeyError`.
* `parse_taxonomy` pads with explicit `""` when a depth level is skipped and
  warns once.

**Efficiency (D):**
* `group_linear_reactions` rewritten with a metabolite worklist (re-enqueue
  the mets touched by each merge); same observable result, O(n+m) work per
  pass instead of restarting the full scan after every merge.
* `parse_kegg_reactions` now caches the parsed stoichiometry on each
  `KeggReaction.stoichiometry`; `build_reference_model` reuses it instead of
  re-parsing.

**Dead code (E):**
* Dropped `KeggReaction.modules` and `.rhea` (parsed but never consumed).
* Dropped the vestigial `only_genes_in_models` parameter from `_ortholog_map`.

Six new regression tests; the only one without a test is the `.part` atomic
download (defensive, needs urlopen mocking).

## Quality sweep — known-issues section B

Closed all four "silent misbehaviour" items from [docs/known_issues.md](docs/known_issues.md):
* `merge_models` warns on `formula` / `charge` conflicts when two source models
  share a name[comp] but disagree (used to silently keep the first-seen).
* `add_reactions_from_equations` warns when creating a metabolite in an
  unregistered compartment — both the `mets_by="id"` and `mets_by="name"` paths
  (id-mode used to skip the check entirely, an asymmetry).
* `parse_task_list` warns when continuation data appears before any task ID
  has been seen (used to silently drop the orphan row).
* `export_model_to_sif` warns up front when a custom label map sends two
  distinct ids to the same label (used to silently collapse nodes).
Four new regression tests cover them.

## Quality sweep — known-issues section A

Closed all six "latent edge-case bug" items from [docs/known_issues.md](docs/known_issues.md):
* `add_reactions_from_equations` no longer misparses `"2 oxoglutarate"` (or any
  leading-number metabolite name) — the resolver tries the full token before
  splitting off a coefficient.
* `add_reactions_from_equations` warns when an equation's terms cancel to a
  zero-metabolite reaction.
* `add_reactions_from_model` tracks ids minted within the batch so two source
  metabolites whose ids both collide with the draft don't collapse onto the
  same generated id.
* `add_transport_reactions` warns on duplicate metabolite names in the source
  or target compartment instead of silently dropping all but one.
* `connect_blocked_reactions` membership-guards the FVA result before
  `.at[]` lookup.
* `assign_kos` rejects `cutoff >= 1` up front — would have crashed inside the
  ratio filter at `log(best_evalue) == 0`.
Six new regression tests cover the user-reachable cases.

## Phase 7 — Localization

* **Sub-cellular localisation by MILP.** [`localization.predict_localization`](src/ravengem/localization/predict.py)
  + [`apply_localization`](src/ravengem/localization/predict.py). Deterministic (not simulated
  annealing); caller-passed `reactions_to_relocate` set with everything else pinned;
  incomplete-model tolerant (no silent reaction removal); `apply=False` returns a diff
  preview; multi-compartment by default with primary-free, extras-penalised scoring.
* **Predictor loaders.** [`load_wolfpsort`, `load_deeploc`](src/ravengem/localization/scores.py),
  with the `gene × compartment` DataFrame contract open for any predictor.
* **Compartment helpers** ([`manipulation/compartments.py`](src/ravengem/manipulation/compartments.py)):
  `merge_compartments`, `copy_to_compartment` — useful standalone for model curation.
* **Real-data validation on yeast-GEM** ([docs/yeast_localization_benchmark.md](docs/yeast_localization_benchmark.md))
  — accuracy 0.72 → 0.39 on 298 GPR'd reactions as confident predictor mis-scoring rises
  from 0 % to 50 %; perfect on compartments with disjoint gene sets (c/g/lp/p/v/vm), and
  surfaces a `transport_cost` calibration insight for soft-probability score tables.

## Phase 5 — Data integration & analysis

* **Reporter metabolites, FSEOF, random sampling** ([`analysis/`](src/ravengem/analysis/)).
* **HPA omics ingestion** ([`omics.parse_hpa`, `parse_hpa_rna`, `hpa_gene_scores`, `rna_gene_scores`](src/ravengem/omics/hpa.py))
  — pandas-tidy DataFrames replace RAVEN's sparse-matrix layout; scoring adapters reuse the
  existing GPR walk.
* **N-model comparison** ([`comparison.compare_models`](src/ravengem/comparison/compare.py)).
* **Dynamic FBA** is **not ported** — established Python packages cover it (`dfba`,
  `reframed`, `mewpy`).

## Phase 4d — ftINIT

* **ftINIT pipeline** ([`init.ftinit`](src/ravengem/init/ftinit.py)) — staged MILP, linear merge,
  task-aware gap-filling, gene pruning.
* **Validated against MATLAB RAVEN on Human-GEM.** 5 Hart2015 cell-line models;
  Jaccard 0.973–0.977 (no-task) and 0.978–0.980 (task-constrained). See
  [docs/humangem_validation.md](docs/humangem_validation.md).
* **Parameter calibration & input-robustness study** ([docs/init_param_calibration.md](docs/init_param_calibration.md))
  — `mip_gap=0.01` is the genome-scale full-pipeline sweet spot (~37% faster than 0.001 at
  Jaccard 0.995); pipeline is robust to expression noise (Jaccard 0.92–0.95) but sensitive
  to sparsity (50–70% dropout → Jaccard 0.59–0.71); the task + gap-fill layer keeps the
  essential-task pass-rate at 67–69/69 across the gradient, whereas tINIT-without-it passes
  only 35/69 even on clean data.
* **Cross-solver portability** ([docs/init_solver_benchmark.md](docs/init_solver_benchmark.md))
  + [`tests/test_init_solvers.py`](tests/test_init_solvers.py): Gurobi and GLPK pass at toy
  scale; only Gurobi is viable at genome scale today (HiGHS hits an upstream optlang
  `clone()` bug; GLPK ignores `configuration.timeout` on MIP).
* **Engineering wins surfaced by the genome-scale work:** `check_tasks` and
  `fill_tasks._feasible` rewritten in-place (~12× each); `optlang.symbolics.add` builds
  in the MILP construction (the O(n²) sympy `sum()` blow-up was the original genome-scale
  blocker); bounded gap-fill MILP; `rescaleModelForINIT` ported.

## Phase 4c — tINIT

* **INIT MILP and the tINIT pipeline** ([`init.run_init`](src/ravengem/init/init.py),
  [`init.get_init_model`](src/ravengem/init/build.py)). Clean optlang reformulation;
  RNA-seq scoring via `5·ln(level/ref)`-clamped.

## Phase 4b — Gap-filling

* **Connectivity gap-filling** ([`gapfilling.connect_blocked_reactions`](src/ravengem/gapfilling/fill.py))
  — MILP. Targeted (toward objective) mode delegates to `cobra.gapfill`.

## Phase 4a — Metabolic tasks

* **Task list parsing + `check_tasks`** ([`tasks/`](src/ravengem/tasks/)).

## Phase 3 — Reconstruction

* **Homology-based draft** from a template GEM + BLAST/DIAMOND wrappers
  ([`reconstruction/homology/`](src/ravengem/reconstruction/homology/)) — with structured
  improvements over RAVEN's `getModelFromHomology` (see IMPROVEMENTS H1–H6).
* **KEGG five-step pipeline** ([`reconstruction/kegg/`](src/ravengem/reconstruction/kegg/)):
  dump → parser → HMM library builder → species model → HMM-query draft.
* **MetaCyc reconstruction** **not ported** (and flagged for removal from MATLAB RAVEN —
  see IMPROVEMENTS R-MetaCyc).

## Phase 2 — I/O

* **YAML** aligned to cobra's `!!omap` writer + RAVEN-only fields preserved into `.notes`,
  plus geckopy `ec-*` for enzyme-constrained models
  ([`io/yaml.py`](src/ravengem/io/yaml.py)).
* **SIF**, **Excel export**, and **Standard-GEM `model/<fmt>/…` git layout**
  ([`io/`](src/ravengem/io/)). Excel import intentionally excluded.

## Phase 1 — Foundation

* **GPR / balance / validation / parsing helpers** ([`utils/`](src/ravengem/utils/)) —
  cobra-absent bits only; the rest are cheatsheeted.
* **Manipulation ergonomic layer** ([`manipulation/`](src/ravengem/manipulation/)) —
  add/change/remove/transport/transfer/merge/simplify/variance + adopted transforms.
* **External-binary resolver** ([`binaries.py`](src/ravengem/binaries.py)) — version-pinned
  release-ZIP registry, SHA256-verified cache.

## Phase 0 — Scaffold

* Project structure, packaging, pytest skeleton, license alignment with MATLAB RAVEN
  (GPL-3.0-or-later).
