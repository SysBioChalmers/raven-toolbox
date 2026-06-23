# Changelog

Milestones in the raven-toolbox port. For function-level status see
[docs/raven_migration.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/reference/migration.md); for open work see
[docs/todo.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/reference/todo.md).

## Unreleased

* **Curation triage for localisation.** Added `triage_localization` — an optional companion to
  compartment assignment that ranks the genes/reactions whose localisation is shakiest (low DeepLoc
  confidence, borderline top-two margin, multi-source disagreement, no evidence, low-trust
  compartment, multi-localised), each with a plain-English reason, so a curator knows where to look.
  Returns a `ReviewReport`. `load_deeploc` gained `keep_raw_confidence=True` and `LocalizationScores`
  a `raw_confidence` field (per-gene normalisation otherwise discards the confidence the triage needs).
* **Cross-species DeepLoc benchmark.** Generalised the predictor benchmark to any curated model
  (`scripts/benchmark_deeploc.py --species {yeast,aracore}`, a per-species compartment config) and
  added an independent *Arabidopsis* test ([AraCore study](docs/studies/deeploc_aracore_benchmark.md)):
  DeepLoc 2.1 generalises across kingdoms — **80.3%** overall with the **chloroplast at 89.9%** (the
  organelle yeast could not exercise). The yeast run was refreshed to DeepLoc's slow (ProtT5) model,
  lifting organelle-collapsed accuracy 54.6% → 64.6% and mitochondrial-membrane recall 47% → 86%.
* **Optional raw DeepLoc probabilities.** `load_deeploc` gained `normalise=False` to keep DeepLoc's
  calibrated probabilities instead of rescaling each gene's best compartment to 1.0. A whole-model
  yeast-GEM benchmark ([study](docs/studies/deeploc_normalisation_benchmark.md)) finds normalisation
  is **accuracy-neutral** for compartment assignment (raw does not rescue the contested or
  high-confidence calls); the only reproducible difference is that raw assigns fewer genes to
  multiple compartments — a re-scaling of the existing `transport_cost`/`multi_compartment_penalty`
  knobs, not new signal. So normalisation stays the **default** and `normalise=False` is an opt-in
  for callers wanting the calibrated magnitudes (e.g. the `triage_localization` confidence signal).
* **Fuse and tune localisation evidence.** Added `combine_scores` (weighted-sum consensus of several
  `LocalizationScores`, so agreement across DeepLoc / UniProt / COMPARTMENTS is reinforced), and gave
  `load_deeploc` / `load_mulocdeep` a `min_confidence=` gate (drop unreliable low-confidence genes)
  plus, for `load_deeploc`, `membrane_split={"m":"mm"}` (route mitochondrion to its membrane
  sub-compartment using the transmembrane signal — mito only; ER is not separable). Motivated and
  validated by the [DeepLoc 2.1 yeast-GEM benchmark](docs/studies/deeploc_yeast_benchmark.md).
* **Prepare sequence-predictor input.** Added `prepare_deeploc_input` (plus `fetch_protein_sequences`
  and `write_fasta`) to write a DeepLoc-2.1-ready protein FASTA for a model's genes — sequences
  fetched from UniProtKB, headers set to the gene ids so the predictor output lines up with the model
  and `load_deeploc`. DeepLoc 2.1 has no batch API; the FASTA is chunked at the web server's
  500-sequence limit, and genes without a reviewed sequence are reported. Script:
  `scripts/prepare_deeploc_yeast.py`.
* **Localisation loaders modernised.** Added `load_mulocdeep` (MULocDeep wide tables),
  `load_compartments` (the COMPARTMENTS evidence database), `load_uniprot` (curated UniProtKB
  `Subcellular location` exports) and `fetch_uniprot_localization` (the same via the UniProt REST
  API by organism id), plus `DEFAULT_COMPARTMENT_MAP` to rename predictor labels to
  model compartment ids and collapse synonyms. `load_deeploc` gained a `compartment_map`
  argument. **Removed `load_wolfpsort`** — modern multi-label predictors, the COMPARTMENTS
  database and UniProt supersede the single-label WoLF PSORT caller.

## 0.2.0 — 2026-06-14

Project rename plus KEGG-reconstruction and CI improvements.

* **Project renamed `raven-python` → `raven-toolbox`.** The import package is now
  `raven_toolbox` (was `raven_python`) and the repository moved to
  `SysBioChalmers/raven-toolbox`. Update imports and any `raven-python` git/URL
  references accordingly. (#34)
* **De-novo KEGG query uses `hmmsearch` instead of `hmmscan`.**
  `get_kegg_model_from_sequences` now runs one `hmmsearch` over the concatenated KO
  library (`-Z` set to the profile count, so E-values and KO assignments are identical
  to the previous `hmmscan` path) — the faster, more parallel search direction.
  `ensure_kegg_hmm_library` no longer runs `hmmpress` (just gunzips); the published
  `.hmm.gz` artefact is unchanged. (#32)
* **Domain-mode `get_kegg_model_for_organism_from_artefacts` auto-resolves the
  taxonomy artefact** from the artefact directory, so `"prokaryotes"` /
  `"eukaryotes"` no longer require an explicit `taxonomy=` path. (#31)
* **Test data no longer ships real KEGG records.** The on-disk
  `tests/data/kegg_dump` is replaced by a session fixture that generates a fully
  fictional KEGG-format dump at runtime, so no KEGG-derived data is redistributed. (#33)
* **Removed the visualization stub and the `[visualization]` extra** — an
  unimplemented placeholder. (#30)
* **CI on Node 24** — `actions/checkout@v5`, `actions/setup-python@v6`. (#35)

## 0.1.0 — 2026-06-10

First release with **published, downloadable KEGG artefacts**, plus a cobra-aligned
hardening pass (no behaviour change on well-formed inputs). Highlights:

* **KEGG artefacts published (`kegg116`):** `ensure_kegg_data` /
  `ensure_kegg_hmm_library` fetch version-pinned, SHA256-verified files from the
  GitHub release. Every artefact is **gzip + version-prefixed**
  (`kegg116_<name>.gz`) so MATLAB and Windows read them with the built-in `gunzip`
  (no external tool) — `organism_gene_ko` moved from xz to gzip for this. The core
  model files (reference model + KO/reaction tables) ship as a single
  `kegg116_core.tar.gz` that `ensure_kegg_data` extracts on first use; the HMM
  libraries and `taxonomy` are separate assets. The **HMM
  libraries ship as one gzip concatenated flatfile per domain**
  (`kegg116_<domain>.hmm.gz`); the client decompresses and `hmmpress`-es once on
  first use, cutting the download ~10× versus the pressed index and letting the
  same artefact serve MATLAB RAVEN.
* **Taxonomy + phylogenetic distance:** publish `kegg116_taxonomy.gz` and add
  `reconstruction.kegg.phyl_dist` (with `PhylDist`), a faithful port of RAVEN's
  `getPhylDist` that regenerates the `keggPhylDist` distance matrix from the
  taxonomy file — so GECKO's organism-distance kcat selection needs no MATLAB
  `.mat`. `ensure_kegg_taxonomy` fetches the artefact.
* **Packaging:** `raven_toolbox.__version__` now derives from the installed package
  metadata (`importlib.metadata`) instead of a hard-coded literal that had drifted
  to `0.0.1`; the docs site reported the wrong version. Pinned `ruff==0.15.15` in
  both the `dev` extra and CI so the lint result is reproducible, and fixed two
  lint errors the unpinned ruff had started flagging.
* **Errors aligned to cobra:** solver/feasibility failures in `run_init`,
  `run_ftinit`, `fill_tasks` and `random_sampling` now raise
  `cobra.exceptions.OptimizationError` (already used elsewhere in the package)
  instead of a bare `RuntimeError`.
* **Consistency:** a single `utils.parse.subsystem_to_str` coerces a reaction
  `subsystem` to cobra's canonical `str` everywhere it is rendered/compared
  (`io.excel`, `comparison.compare`, `curation.batch`, `manipulation.add`) — fixes
  a crash on non-string subsystem items and the silent drop of multi-subsystem
  reactions. GPR score-aggregation (`AGGREGATORS` / `resolve_aggregators`) is now
  shared by `init.score` and `init.genes`. Maintainer-side KEGG-download progress
  uses a module logger instead of `print`.
* **Robustness:** path-traversal guard on bundled-ZIP extraction (`binaries.py`,
  matching the tarfile `filter="data"` precedent); `connect_blocked_reactions`
  rejects a non-positive `penalty`; `random_sampling` refuses a NaN-contaminated
  sample matrix; `ec_data` warns on an all-zero reaction↔enzyme coupling; optional
  `verify=` SHA256 re-check on `ensure_data_file` cache hits; reporter p-value
  guarded against non-finite z-scores. Regression tests added for each.

## 0.1.0a1 — 2026-05-30

First alpha release. Covers the functional scope of RAVEN built on cobrapy:
de-novo reconstruction (KEGG / homology), context-specific modeling (tINIT / ftINIT),
metabolic-task validation, connectivity gap-filling, HPA omics ingestion, sub-cellular
localisation, N-model comparison, reporter metabolites, FSEOF, flux sampling, and the
RAVEN-style I/O formats (YAML / SIF / Excel). Validated against MATLAB RAVEN on Human-GEM
(Jaccard 0.975–0.980).

* **Licensing:** released under the **MIT** license (previously GPL-3.0-or-later).
* **Docs:** Sphinx + MyST documentation site (sources under `docs/`).
* Not yet implemented: visualization (`visualization/`), metabolomics-based (f)tINIT scoring,
  and published binary / KEGG-artefact release bundles. See the README and
  [docs/todo.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/reference/todo.md).

The milestone sections below record the incremental development history leading to this release.

## Infrastructure

* **GitHub Actions CI** ([.github/workflows/ci.yml](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/.github/workflows/ci.yml)) —
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
`test_utils_balance.py`). [docs/known_issues.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/reference/known_issues.md) now
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

Closed all four "silent misbehaviour" items from [docs/known_issues.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/reference/known_issues.md):
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

Closed all six "latent edge-case bug" items from [docs/known_issues.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/reference/known_issues.md):
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

* **Sub-cellular localisation by MILP.** [`localization.predict_localization`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/localization/predict.py)
  + [`apply_localization`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/localization/predict.py). Deterministic (not simulated
  annealing); caller-passed `reactions_to_relocate` set with everything else pinned;
  incomplete-model tolerant (no silent reaction removal); `apply=False` returns a diff
  preview; multi-compartment by default with primary-free, extras-penalised scoring.
* **Predictor loaders.** [`load_wolfpsort`, `load_deeploc`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/localization/scores.py),
  with the `gene × compartment` DataFrame contract open for any predictor.
* **Compartment helpers** ([`manipulation/compartments.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/compartments.py)):
  `merge_compartments`, `copy_to_compartment` — useful standalone for model curation.
* **Real-data validation on yeast-GEM** ([docs/yeast_localization_benchmark.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/studies/yeast_localization_benchmark.md))
  — accuracy 0.72 → 0.39 on 298 GPR'd reactions as confident predictor mis-scoring rises
  from 0 % to 50 %; perfect on compartments with disjoint gene sets (c/g/lp/p/v/vm), and
  surfaces a `transport_cost` calibration insight for soft-probability score tables.

## Phase 5 — Data integration & analysis

* **Reporter metabolites, FSEOF, random sampling** ([`analysis/`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/analysis/)).
* **HPA omics ingestion** ([`omics.parse_hpa`, `parse_hpa_rna`, `hpa_gene_scores`, `rna_gene_scores`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/omics/hpa.py))
  — pandas-tidy DataFrames replace RAVEN's sparse-matrix layout; scoring adapters reuse the
  existing GPR walk.
* **N-model comparison** ([`comparison.compare_models`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/comparison/compare.py)).
* **Dynamic FBA** is **not ported** — established Python packages cover it (`dfba`,
  `reframed`, `mewpy`).

## Phase 4d — ftINIT

* **ftINIT pipeline** ([`init.ftinit`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/ftinit.py)) — staged MILP, linear merge,
  task-aware gap-filling, gene pruning.
* **Validated against MATLAB RAVEN on Human-GEM.** 5 Hart2015 cell-line models;
  Jaccard 0.973–0.977 (no-task) and 0.978–0.980 (task-constrained). See
  [docs/humangem_validation.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/studies/humangem_validation.md).
* **Parameter calibration & input-robustness study** ([docs/init_param_calibration.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/studies/init_param_calibration.md))
  — `mip_gap=0.01` is the genome-scale full-pipeline sweet spot (~37% faster than 0.001 at
  Jaccard 0.995); pipeline is robust to expression noise (Jaccard 0.92–0.95) but sensitive
  to sparsity (50–70% dropout → Jaccard 0.59–0.71); the task + gap-fill layer keeps the
  essential-task pass-rate at 67–69/69 across the gradient, whereas tINIT-without-it passes
  only 35/69 even on clean data.
* **Cross-solver portability** ([docs/init_solver_benchmark.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/studies/init_solver_benchmark.md))
  + [`tests/test_init_solvers.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/tests/test_init_solvers.py): Gurobi and GLPK pass at toy
  scale; only Gurobi is viable at genome scale today (HiGHS hits an upstream optlang
  `clone()` bug; GLPK ignores `configuration.timeout` on MIP).
* **Engineering wins surfaced by the genome-scale work:** `check_tasks` and
  `fill_tasks._feasible` rewritten in-place (~12× each); `optlang.symbolics.add` builds
  in the MILP construction (the O(n²) sympy `sum()` blow-up was the original genome-scale
  blocker); bounded gap-fill MILP; `rescaleModelForINIT` ported.

## Phase 4c — tINIT

* **INIT MILP and the tINIT pipeline** ([`init.run_init`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/init.py),
  [`init.get_init_model`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/build.py)). Clean optlang reformulation;
  RNA-seq scoring via `5·ln(level/ref)`-clamped.

## Phase 4b — Gap-filling

* **Connectivity gap-filling** ([`gapfilling.connect_blocked_reactions`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/gapfilling/fill.py))
  — MILP. Targeted (toward objective) mode delegates to `cobra.gapfill`.

## Phase 4a — Metabolic tasks

* **Task list parsing + `check_tasks`** ([`tasks/`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/tasks/)).

## Phase 3 — Reconstruction

* **Homology-based draft** from a template GEM + BLAST/DIAMOND wrappers
  ([`reconstruction/homology/`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/reconstruction/homology/)) — with structured
  improvements over RAVEN's `getModelFromHomology` (see IMPROVEMENTS H1–H6).
* **KEGG five-step pipeline** ([`reconstruction/kegg/`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/reconstruction/kegg/)):
  dump → parser → HMM library builder → species model → HMM-query draft.
* **MetaCyc reconstruction** **not ported** (and flagged for removal from MATLAB RAVEN —
  see IMPROVEMENTS R-MetaCyc).

## Phase 2 — I/O

* **YAML** aligned to cobra's `!!omap` writer + RAVEN-only fields preserved into `.notes`,
  plus geckopy `ec-*` for enzyme-constrained models
  ([`io/yaml.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/io/yaml.py)).
* **SIF**, **Excel export**, and **Standard-GEM `model/<fmt>/…` git layout**
  ([`io/`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/io/)). Excel import intentionally excluded.

## Phase 1 — Foundation

* **GPR / balance / validation / parsing helpers** ([`utils/`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/utils/)) —
  cobra-absent bits only; the rest are cheatsheeted.
* **Manipulation ergonomic layer** ([`manipulation/`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/)) —
  add/change/remove/transport/transfer/merge/simplify/variance + adopted transforms.
* **External-binary resolver** ([`binaries.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/binaries.py)) — version-pinned
  release-ZIP registry, SHA256-verified cache.

## Phase 0 — Scaffold

* Project structure, packaging, pytest skeleton, license alignment with MATLAB RAVEN
  (GPL-3.0-or-later).
