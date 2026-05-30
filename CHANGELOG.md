# Changelog

Milestones in the raven-python port. For function-level status see
[docs/raven_migration.md](https://github.com/SysBioChalmers/raven-python/blob/develop/docs/reference/migration.md); for open work see
[docs/todo.md](https://github.com/SysBioChalmers/raven-python/blob/develop/docs/reference/todo.md).

## 0.1.0a1 — 2026-05-30

First alpha release. Covers the functional scope of RAVEN built on cobrapy:
de-novo reconstruction (KEGG / homology), context-specific modeling (tINIT / ftINIT),
metabolic-task validation, connectivity gap-filling, HPA omics ingestion, sub-cellular
localisation, N-model comparison, reporter metabolites, FSEOF, flux sampling, and the
RAVEN-style I/O formats (YAML / SIF / Excel). Validated against MATLAB RAVEN on Human-GEM
(Jaccard 0.975–0.980).

* **Licensing:** released under the **MIT** license (previously GPL-3.0-or-later).
* **Docs:** Sphinx + MyST documentation site (sources under `docs/`).
* Not yet implemented: visualization (`plotting/`), metabolomics-based (f)tINIT scoring,
  and published binary / KEGG-artefact release bundles. See the README and
  [docs/todo.md](https://github.com/SysBioChalmers/raven-python/blob/develop/docs/reference/todo.md).

The milestone sections below record the incremental development history leading to this release.

## Infrastructure

* **GitHub Actions CI** ([.github/workflows/ci.yml](https://github.com/SysBioChalmers/raven-python/blob/develop/.github/workflows/ci.yml)) —
  ruff + pytest matrix over Python 3.11/3.12/3.13. Tests that require Gurobi
  auto-skip (no Gurobi on free runners); the known HiGHS upstream blocker
  (`hybrid_interface.Configuration` rejects `lp_method='primal'`) is marked
  `xfail(strict=True)` so CI flips red when optlang fixes it.

## Quality sweep — known-issues section F (design-choice divergences)

Closed the five items in section F (the "design choices that differ from RAVEN"
backlog from the original review). Three docstring/comment fixes; two code
fixes with matching MATLAB back-port proposals in IMPROVEMENTS.md (FS4, B2).

* `run_init` docstring spells out the score-0 semantics divergence between
  classic INIT and ftINIT.
* `get_init_model` inaccurate "same regime" comment replaced with an accurate
  description of the conservative pre-filter.
* `fseof` classifier now uses the slope of `|flux|` (`linregress(enforced, |flux|)`)
  instead of first-vs-last endpoints. A track whose endpoints straddle a
  peak/trough no longer ends up mislabelled.
* `reporter_metabolites` docstring documents the one-sided p-value + z-score
  ordering vs RAVEN's two-tailed sort, and points at the up/down split via
  `gene_fold_changes`.
* `get_elemental_balance` now reports `unknown` for empty-stoichiometry
  reactions (previously vacuously `balanced`). Original review attributed the
  bug to `check_model`; the actual code is in `balance.py`.

Two new regression tests (F3 in `test_analysis_fseof.py`, F5 in
`test_utils_balance.py`). [docs/known_issues.md](https://github.com/SysBioChalmers/raven-python/blob/develop/docs/reference/known_issues.md) now
fully closed (all sections A–F).

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

Closed all four "silent misbehaviour" items from [docs/known_issues.md](https://github.com/SysBioChalmers/raven-python/blob/develop/docs/reference/known_issues.md):
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

Closed all six "latent edge-case bug" items from [docs/known_issues.md](https://github.com/SysBioChalmers/raven-python/blob/develop/docs/reference/known_issues.md):
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

* **Sub-cellular localisation by MILP.** [`localization.predict_localization`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/localization/predict.py)
  + [`apply_localization`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/localization/predict.py). Deterministic (not simulated
  annealing); caller-passed `reactions_to_relocate` set with everything else pinned;
  incomplete-model tolerant (no silent reaction removal); `apply=False` returns a diff
  preview; multi-compartment by default with primary-free, extras-penalised scoring.
* **Predictor loaders.** [`load_wolfpsort`, `load_deeploc`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/localization/scores.py),
  with the `gene × compartment` DataFrame contract open for any predictor.
* **Compartment helpers** ([`manipulation/compartments.py`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/manipulation/compartments.py)):
  `merge_compartments`, `copy_to_compartment` — useful standalone for model curation.
* **Real-data validation on yeast-GEM** ([docs/yeast_localization_benchmark.md](https://github.com/SysBioChalmers/raven-python/blob/develop/docs/studies/yeast_localization_benchmark.md))
  — accuracy 0.72 → 0.39 on 298 GPR'd reactions as confident predictor mis-scoring rises
  from 0 % to 50 %; perfect on compartments with disjoint gene sets (c/g/lp/p/v/vm), and
  surfaces a `transport_cost` calibration insight for soft-probability score tables.

## Phase 5 — Data integration & analysis

* **Reporter metabolites, FSEOF, random sampling** ([`analysis/`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/analysis/)).
* **HPA omics ingestion** ([`omics.parse_hpa`, `parse_hpa_rna`, `hpa_gene_scores`, `rna_gene_scores`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/omics/hpa.py))
  — pandas-tidy DataFrames replace RAVEN's sparse-matrix layout; scoring adapters reuse the
  existing GPR walk.
* **N-model comparison** ([`comparison.compare_models`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/comparison/compare.py)).
* **Dynamic FBA** is **not ported** — established Python packages cover it (`dfba`,
  `reframed`, `mewpy`).

## Phase 4d — ftINIT

* **ftINIT pipeline** ([`init.ftinit`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/init/ftinit.py)) — staged MILP, linear merge,
  task-aware gap-filling, gene pruning.
* **Validated against MATLAB RAVEN on Human-GEM.** 5 Hart2015 cell-line models;
  Jaccard 0.973–0.977 (no-task) and 0.978–0.980 (task-constrained). See
  [docs/humangem_validation.md](https://github.com/SysBioChalmers/raven-python/blob/develop/docs/studies/humangem_validation.md).
* **Parameter calibration & input-robustness study** ([docs/init_param_calibration.md](https://github.com/SysBioChalmers/raven-python/blob/develop/docs/studies/init_param_calibration.md))
  — `mip_gap=0.01` is the genome-scale full-pipeline sweet spot (~37% faster than 0.001 at
  Jaccard 0.995); pipeline is robust to expression noise (Jaccard 0.92–0.95) but sensitive
  to sparsity (50–70% dropout → Jaccard 0.59–0.71); the task + gap-fill layer keeps the
  essential-task pass-rate at 67–69/69 across the gradient, whereas tINIT-without-it passes
  only 35/69 even on clean data.
* **Cross-solver portability** ([docs/init_solver_benchmark.md](https://github.com/SysBioChalmers/raven-python/blob/develop/docs/studies/init_solver_benchmark.md))
  + [`tests/test_init_solvers.py`](https://github.com/SysBioChalmers/raven-python/blob/develop/tests/test_init_solvers.py): Gurobi and GLPK pass at toy
  scale; only Gurobi is viable at genome scale today (HiGHS hits an upstream optlang
  `clone()` bug; GLPK ignores `configuration.timeout` on MIP).
* **Engineering wins surfaced by the genome-scale work:** `check_tasks` and
  `fill_tasks._feasible` rewritten in-place (~12× each); `optlang.symbolics.add` builds
  in the MILP construction (the O(n²) sympy `sum()` blow-up was the original genome-scale
  blocker); bounded gap-fill MILP; `rescaleModelForINIT` ported.

## Phase 4c — tINIT

* **INIT MILP and the tINIT pipeline** ([`init.run_init`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/init/init.py),
  [`init.get_init_model`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/init/build.py)). Clean optlang reformulation;
  RNA-seq scoring via `5·ln(level/ref)`-clamped.

## Phase 4b — Gap-filling

* **Connectivity gap-filling** ([`gapfilling.connect_blocked_reactions`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/gapfilling/fill.py))
  — MILP. Targeted (toward objective) mode delegates to `cobra.gapfill`.

## Phase 4a — Metabolic tasks

* **Task list parsing + `check_tasks`** ([`tasks/`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/tasks/)).

## Phase 3 — Reconstruction

* **Homology-based draft** from a template GEM + BLAST/DIAMOND wrappers
  ([`reconstruction/homology/`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/reconstruction/homology/)) — with structured
  improvements over RAVEN's `getModelFromHomology` (see IMPROVEMENTS H1–H6).
* **KEGG five-step pipeline** ([`reconstruction/kegg/`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/reconstruction/kegg/)):
  dump → parser → HMM library builder → species model → HMM-query draft.
* **MetaCyc reconstruction** **not ported** (and flagged for removal from MATLAB RAVEN —
  see IMPROVEMENTS R-MetaCyc).

## Phase 2 — I/O

* **YAML** aligned to cobra's `!!omap` writer + RAVEN-only fields preserved into `.notes`,
  plus geckopy `ec-*` for enzyme-constrained models
  ([`io/yaml.py`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/io/yaml.py)).
* **SIF**, **Excel export**, and **Standard-GEM `model/<fmt>/…` git layout**
  ([`io/`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/io/)). Excel import intentionally excluded.

## Phase 1 — Foundation

* **GPR / balance / validation / parsing helpers** ([`utils/`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/utils/)) —
  cobra-absent bits only; the rest are cheatsheeted.
* **Manipulation ergonomic layer** ([`manipulation/`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/manipulation/)) —
  add/change/remove/transport/transfer/merge/simplify/variance + adopted transforms.
* **External-binary resolver** ([`binaries.py`](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/binaries.py)) — version-pinned
  release-ZIP registry, SHA256-verified cache.

## Phase 0 — Scaffold

* Project structure, packaging, pytest skeleton, license alignment with MATLAB RAVEN
  (GPL-3.0-or-later).
