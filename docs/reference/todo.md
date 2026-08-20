# Open work

The item-level backlog. For the **plan** — phases, dependencies, exit criteria and the open
decisions — see [roadmap.md](roadmap.md); this page is the detail behind it.

Items are tagged **P0** (do next), **P1** (should land before 1.0), **P2** (nice to have /
watching brief), and by where the work lands:

* **→ raven-docs** — the dual-language documentation site
  ([edkerk/raven-docs](https://github.com/edkerk/raven-docs), MkDocs + Material on Read the
  Docs; RAVEN and raven-toolbox are git submodules, tracking `develop3` and `develop`).
  User-facing prose belongs there so it is online and covers both languages. Its
  `DESIGN.md` is the source of intent — read it before writing pages, and record decisions
  in its decisions log.
* **→ raven-toolbox** — code, tests, CI, and the internal reference documents that back the
  public pages.

**Page-granularity rule for raven-docs:** one topic per page, kept short. Prefer several
small linked pages over one long scrolling document; the site's left sidebar is the index.
The existing reconstruction protocol is the model to follow — one protocol is a *nav
section* of ~40–95-line chapter pages, not a single page.

Background references in this repo: [migration.md](migration.md) (function-by-function port
status), [matlab_raven_backports.md](matlab_raven_backports.md) (MATLAB ↔ Python differences
record), [improvements.md](improvements.md) (design decisions + proposed improvements),
[known_issues.md](known_issues.md) (review backlog — sections A–F are all closed).

---

## 1. MATLAB ↔ Python differences (→ raven-docs)

`docs/differences.md` on the site is a single page that mixes orientation, function mapping,
per-language feature lists and solver notes — and it is **currently wrong in ways that will
mislead users**:

* It lists Python functions that do not exist: `get_blast`, `get_diamond` (real names
  `run_blast` / `run_diamond`), `fill_gaps` (real names `connect_blocked_reactions`,
  `fill_gaps_fast_lp`, `fill_gaps_kumar_milp`), `import_model`, `export_model`,
  `import_excel_model`, `export_to_excel_format` (real name `export_to_excel`), `solve_lp`.
* It states that both toolboxes read and write RAVEN Excel "making models fully portable".
  raven-toolbox exports Excel but **deliberately has no Excel import**.
* It links twice to `matlab-vs-python.md`, which does not exist — and `mkdocs.yml` lists
  that same missing file as a **nav entry** under "MATLAB vs Python", so the build emits a
  broken nav reference.

The underlying cause is that the page was written by analogy (camelCase → snake_case)
instead of from the source, so the fix has to be structural, not just a correction pass.

* **P0 — Split `differences.md` into a small section**, one page per question, replacing the
  single page:
  * `differences/index.md` — orientation only: two independent implementations, large
    overlap, why the Python one is built on cobrapy. Short, with cards to the pages below.
  * `differences/mapping.md` — the full MATLAB ↔ Python function table the nav and the page
    body already link to as `matlab-vs-python.md`. Per `DESIGN.md` §4.3: RAVEN ↔
    raven-toolbox pairs auto-generated from the API data, cobrapy-alternative rows
    hand-curated.
  * `differences/python-only.md` — functionality with no MATLAB counterpart. From this
    repo's `__all__`s, 121 of 185 public names have no row in
    [migration.md](migration.md) today, including whole subsystems: confidence tracking,
    biomass helpers, growth conditions, batch curation, ΔG/SBO annotation, transport
    evidence, `assign_compartments`, `diff_models`, the KEGG artefact builders.
  * `differences/matlab-only.md` — `ravenCobraWrapper`, the `drawMap` family, MetaCyc
    reconstruction (flagged for removal upstream), dynamic FBA, ftINIT metabolomics scoring,
    Excel import, `printFluxes`; each with the reason it is absent and what to use instead.
  * `differences/behaviour.md` — **same function, different answer.** The page that does not
    exist anywhere today and is the most valuable one: differing defaults, arguments, return
    shapes, ordering/tie-breaking, and solver dependence. Rows already known from this repo:
    `check_tasks` (one model reused vs copy-per-task), `reporter_metabolites` (one-sided
    p-value, z-sorted vs RAVEN's two-tailed ordering), `fseof` (abs-slope classifier),
    `get_elemental_balance` (graded `unknown` class), `run_init` vs `ftinit` score-0
    semantics, `merge_models` / `add_reactions_from_model` (`name[comp]` matching),
    `convert_to_irreversible` / `expand_model` (geckopy-derived), `write_yaml_model`
    (`!!omap`, `metaData` first), `diff_models` (order-insensitive grRule logic), plus the
    two the protocol port will surface: `fillGaps` vs the three Python gap-fillers, and
    `contractModel` vs `remove_duplicate_reactions`.
  * `differences/parity.md` — what "the same result" means per function class: exact,
    set-level (alternate MILP optima), or statistical. Pairs with the test tiers in §2.
* **P0 — Generate the mapping, and validate names against the submodules.** A build-time
  generator (raven-docs `scripts/`, alongside `gen_api_pages.py`) that emits the pairs and
  **fails the build on any hand-written function name that no longer resolves** in either
  submodule. This is what prevents a repeat of the wrong names above.
* **P1 — Pin the baseline.** Record which RAVEN commit each differences page was verified
  against, plus a "last verified" date, so drift is visible.
* **P1 — Keep the source-of-truth split clear.** raven-docs carries the user-facing pages;
  this repo keeps [migration.md](migration.md) (port decisions),
  [matlab_raven_backports.md](matlab_raven_backports.md) (upstream port plans) and
  [improvements.md](improvements.md) (rationale) as the maintainer record the pages are
  derived from. Cross-link both ways; do not duplicate prose.

## 2. Cross-language equivalence tests (→ raven-toolbox)

What exists today: hand-transcribed RAVEN oracles (`tests/tinit_oracles.py`, from
`tinitTests.m`), the YAML round-trip parity gate (`tests/test_io_yaml_parity.py`), and the
validation studies (Human-GEM Jaccard 0.975–0.980, yeast, multi-organism) which are
*reported* in `docs/studies/` but never *asserted* in CI. Nothing regenerates MATLAB output
automatically, so a parity claim cannot fail a build.

* **P0 — Golden-fixture harness.** Committed small inputs plus a MATLAB driver that runs the
  RAVEN counterpart and writes JSON/CSV oracles into `tests/data/matlab/`; pytest loads
  those and asserts equality. Fixtures small enough to commit and to run without Gurobi.
  *Done when:* `pytest -m parity` passes on a clean checkout with no MATLAB installed, and a
  documented regeneration step reruns MATLAB to refresh the oracles.
* **P0 — Tier the parity contract**, because "identical output" is not achievable uniformly.
  The tiers are also what `differences/parity.md` publishes:
  1. **Exact** — YAML/Excel/SIF export, task-list parsing, GPR normalisation, elemental
     balance, `sort_identifiers`, `merge_models`, `convert_to_irreversible`, `expand_model`,
     KEGG table parsing, homology ortholog maps.
  2. **Set-level** — MILP outcomes ((f)tINIT extraction, gap-filling, localisation): assert
     Jaccard/containment bands against a recorded baseline rather than identity, since
     alternate optima are legitimate.
  3. **Statistical** — flux sampling and random sampling: distributional checks at a fixed
     seed.
* **P1 — The new tutorial is itself a parity fixture.** Written dual-language from the
  start (§3), it gives an end-to-end reconstruction with a MATLAB lane and a Python lane
  over identical inputs — the most realistic equivalence test available, and far better
  evidence than toy models. Design its fixtures so the nightly job can run both lanes and
  diff the resulting models. (The published hanpo-GEM protocol, with its known anticipated
  results, remains a useful *reference* benchmark even though the tutorial is not based on
  it.)
* **P1 — Determinism regression tests.** Recent fixes (#76, #83, `c239d2e`) made placement
  and gap-fill deterministic, but no CI test would catch a regression. Add repeated-run
  identity assertions for `assign_compartments`, `predict_localization`, and
  `run_init`/`ftinit` on toy models, plus row/column ordering of the built MILPs. The
  untracked `scripts/determinism_probe.py` / `master_determinism_probe.py` are the starting
  point.
* **P1 — Solver-dependent parity job.** Genome-scale (f)tINIT needs Gurobi, which free
  runners cannot install. A nightly / manually-triggered workflow on a licensed runner that
  runs the tier-2 checks and reports the Jaccard numbers, so the study documents stop being
  hand-refreshed.
* **P1 — Promote or delete the 14 untracked scripts** now sitting in `scripts/`
  (`cross_py_on_mat.py`, `full_pipeline_py.py`, `diff_drafts.py`, `export_draft*.py`,
  `export_scope.py`, `param_sweep.py`, `py_mps.py`, `roworder_test.py`, …). Several are
  exactly the cross-language drivers the harness above needs; the rest are spent scratch
  work. *Done when:* `git status` is clean and `scripts/README.md` documents what survived.
* **P2 — Assert the validation-study headline numbers** (Jaccard bands, task pass counts) as
  tolerance-checked tests in the nightly job, so a study document can never disagree with
  the code.

## 3. Protocols and tutorials (→ raven-docs)

### What is actually on the site

Reconciled against the checkout, not the repository listing:

* **One protocol exists, in ten chapter pages** — the hanpo-GEM homology reconstruction
  (Zorrilla & Kerkhoven 2022, *MMB* 2513), nav section "GEM reconstruction (*H. polymorpha*)":
  introduction → materials → template models → homology → biomass → lipid curation →
  gap-filling → save & simulate → manual curation → anticipated results. Its authoritative
  source is `code/reconstructionProtocol.m` in the hanpo-GEM submodule, transcluded whole
  into the last page.
* **Legacy tutorials 1–5** — MATLAB only, RAVEN 1 exercises, staying untouched.
* **"GEM extraction" and "GEM comparison" are listed as planned** on `protocol/index.md`,
  with nothing behind them yet.

### The finding that reorders this section

**There is no Python anywhere in the guides.** Across the ten protocol chapters and the five
legacy tutorials there are 72 MATLAB code blocks, **zero Python blocks and zero
MATLAB/Python tabs** — while `protocol/index.md` tells the reader protocols "are provided in
**both MATLAB and Python** with code shown in MATLAB/Python tabs", and `DESIGN.md` §4.4 and
§5 require linked, persistent tabs site-wide. The site currently promises a dual-language
protocol it does not have.

### What we are doing about it

**The existing protocol is not ported.** It stays as it is — a faithful rendering of a
published MATLAB protocol — alongside the untouched legacy tutorials. Instead we **define new
tutorials, dual-language from the first line**. The homology one may resemble the hanpo-GEM
protocol in structure and in several individual steps (homology draft → biomass →
gap-filling → simulation → curation is the natural arc of any reconstruction), but it is its
own exercise, with its own organism, its own inputs, and its own runnable scripts in both
languages.

**SLIME lipid curation is out of scope** while the parallel work on those methods proceeds.
The new tutorials do not depend on it and do not attempt a Python equivalent — which removes
what would otherwise have been the one hard blocker.

* **P0 — Correct the claim on `protocol/index.md` now.** Until Python material exists, the
  index must not state that protocols are provided in both languages. A one-line fix that
  stops the site over-promising.
* **P0 — Two new reconstruction tutorials, one per route** (roadmap decision **D1**), each a
  nav section of chapter pages with MATLAB/Python tabs throughout:
  * **Homology — a non-conventional yeast.** First: well-understood arc, small inputs, no new
    infrastructure. Stops short of lipid curation; biomass without the SLIME treatment.
  * **De novo from KEGG — a small prokaryote.** Second: the route with no coverage on the
    site today, but it must be explicit about artefact downloads and their size.
  Still to pick: the specific organisms and templates, and the published phenotype each
  tutorial validates against.
* **Feasibility, measured.** I checked every MATLAB call in the existing protocol's ten
  chapters against this package, as a proxy for what a Python lane of *any* reconstruction
  tutorial needs: **seven of the eight code-bearing chapters are already portable**, and the
  one that is not is the SLIME chapter now out of scope.

  | Chapter | MATLAB | Python | Status |
  |---|---|---|---|
  | Template models | `importModel`, `exportToExcelFormat`, `exportModel` | `cobra.io.read_sbml_model`, `export_to_excel`, `cobra.io.write_sbml_model` | ✅ portable |
  | Draft from homology | `getBlast`, `getModelFromHomology`, `contractModel`, `addRxnsGenesMets` | `run_blast`, `get_model_from_homology`, `remove_duplicate_reactions`, `add_reactions_from_model` | ✅ portable (`contractModel` ≠ `remove_duplicate_reactions` exactly — behaviour row) |
  | Biomass composition | `addRxnsGenesMets`, `changeRxns` | `add_reactions_from_model`, `change_reaction_equations` | ✅ portable |
  | Lipid curation (SLIME) | `addLipidReactions`, `addSLIMEreactions`, `scaleLipids` | — | ⏸️ **out of scope** — SLIMEr / hanpo-GEM helpers, not RAVEN core; parallel method work ongoing, so the new tutorial avoids this ground |
  | Gap-filling | `setParam`, `getExchangeRxns`, `removeReactions`, `fillGaps`, `solveLP`, `printFluxes`, `deleteUnusedGenes` | bounds/objective one-liners, `model.exchanges`, `remove_reactions`, `connect_blocked_reactions` / `fill_gaps_fast_lp` / `fill_gaps_kumar_milp`, `model.optimize`, — , `remove_genes` | ⚠️ portable, but `fillGaps` maps to three functions with different semantics, and `printFluxes` has no equivalent |
  | Save and simulate | `model.annotation.*`, `newCommit`, `newRelease`, `solveLP` | YAML `metaData`, `export_for_git`, `model.optimize` | ⚠️ portable; `newCommit`/`newRelease` are hanpo-GEM repo helpers over `exportForGit`, so the Python lane needs its own equivalents |
  | Manual curation | `changeGrRules`, `addRxns`, `addTransport`, `addRxnsGenesMets` | `change_gene_reaction_rules`, `add_reactions_from_equations`, `add_transport_reactions`, `add_reactions_from_model` | ✅ portable |
  | Gap-analysis tip | `canProduce`, `canConsume`, `checkProduction`, `getAllSubGraphs`, `haveFlux` | `analyse_topology` + `cobra.flux_analysis.find_blocked_reactions` | ✅ portable |

* **P1 — Then the two already-planned protocols**, each as its own nav section of chapter
  pages: **GEM extraction** (ftINIT — adapt the existing good material from the
  Human-GEM-guide rather than writing fresh, per the maintainer's instruction) and **GEM
  comparison** (`compare_models` / `diff_models`, reporter metabolites).
* **P1 — New protocols beyond those.** Revised against what the reconstruction chapters
  already cover — my earlier "draft to a model that grows" and "simulation and analysis"
  proposals are **dropped**, since the gap-filling and save-and-simulate chapters cover that
  ground:
  * **Compartmentalisation** — `predictLocalization` vs `assign_compartments`, evidence
    sources, reading the certification report. Largely Python-only.
  * **Confidence tracking** — running the facets, using the bands for curation triage.
    Python-only; currently only a study in this repo.
  * **ecModel handoff** — YAML `ec-*` round-trip with geckopy.
* **P1 — Decide how Python-only protocols present.** Several of the above have no MATLAB
  counterpart, which the MATLAB-default linked-tab convention does not cover. Options: a
  "Python only" badge with the MATLAB tab explaining the absence, or a separate Python-only
  grouping. Record the decision in the raven-docs decisions log.
* **P1 — Keep tutorial code executable.** The existing protocol transcludes a real `.m`
  script, so its MATLAB lane cannot rot silently. The new tutorial needs that on both sides:
  a runnable MATLAB script and a runnable Python script, version-controlled and transcluded
  into the pages, exercised in CI (roadmap decision **D4** fixes where they live). Without
  it the Python lane will drift exactly the way `differences.md` did.
* **P2 — Retire the thin capability guides in this repo** (`docs/guide/*.md`, 18–69 lines
  each) once the protocols cover the same ground, leaving the API reference and the studies
  here and the narrative on the site.

## 4. Maturity checklist (→ raven-toolbox, surfaced in raven-docs)

There is no single view of which functions are production-ready. The information exists
(tests, studies, `known_issues.md`) but is not aggregated, so users cannot tell a
five-times-validated path from a thin wrapper.

* **P0 — Per-function maturity table** with an explicit rubric — **stable** (unit-tested +
  validated on a real model + API frozen), **provisional** (tested, API may change),
  **experimental** (works, little validation). Maintained here; published as its own small
  page in the raven-docs differences section.
* **P0 — Coverage report in CI** with a floor, plus a first pass over the thin spots.
  Modules not obviously exercised by any test module today, to confirm or fix:
  `localization/transport_evidence.py`, `localization/triage.py`,
  `localization/substrate_ontology.py`, `localization/transporter_tables.py`,
  `manipulation/parameters.py`, `manipulation/boundary.py`, `biomass/*`,
  `conditions/apply.py`, `curation/batch.py`, `annotation/delta_g.py`, `annotation/sbo.py`,
  `reconstruction/kegg/assemble.py`, `reconstruction/homology/hits.py`.
* **P1 — Gaps the protocol port exposes.** `printFluxes` has no Python equivalent (a small,
  obviously useful addition); `contractModel`'s grRule-merging semantics differ from
  `remove_duplicate_reactions`; there is no Python counterpart to hanpo-GEM's
  `newCommit`/`newRelease` convenience layer over `export_for_git`. Decide port-or-document
  for each.
* **P1 — Docstring completeness for the generated API.** raven-docs renders this package's
  NumPy docstrings via `mkdocstrings`/griffe, so a thin docstring here is a thin page there.
  Audit the subsystems with no local API page — `confidence`, `biomass`, `conditions`,
  `curation`, `annotation`, `manifest`, `data`, `binaries` — and add them to
  `docs/reference/api/`.
* **P1 — State the known functional gaps in one place**: ftINIT metabolomics scoring
  (`NotImplementedError`), no Excel import, no MetaCyc reconstruction, no dynamic FBA,
  genome-scale (f)tINIT effectively requiring Gurobi. Feeds `differences/matlab-only.md`.
* **P1 — Implement or formally drop the 💡 proposals** in `improvements.md` (A4 compartment
  inference from structured metabolite ids, Y4 a first-class home for
  `deltaG`/`confidence_score`, R4 the `remove_metabolites` wrapper review, G7
  follow-through).
* **P2 — Stability policy for 0.x → 1.0**: what counts as public API, the deprecation
  window, and how `notes` / `annotation` key names are versioned.

## 5. Housekeeping

* **P1 — Land or close the open PRs**:
  [#75](https://github.com/SysBioChalmers/raven-toolbox/pull/75) (assignment ablation
  benchmark) and [#85](https://github.com/SysBioChalmers/raven-toolbox/pull/85) (opt-in
  deterministic ftINIT extraction, still draft).
* **P1 — Post-0.3.0 doc debt**: restore the documentation that PR #53 removed, and land the
  unmerged empirical parameter-default results.
* **P2 — MATLAB-side back-port queue**: `assignCompartments.m` (plan already written in
  `matlab_raven_backports.md`), plus the FS4 / B2 / G1 / G5 items.

## Upstream blockers (not raven-toolbox work, but worth tracking)

* `optlang.hybrid_interface.Configuration.clone()` bug — blocks HiGHS at any scale (CI
  catches it in `tests/test_init_solvers.py`).
* GLPK's MIP solve ignores `configuration.timeout` at genome scale — blocks GLPK on large
  MILPs.
* Both documented in [init_solver_benchmark.md](../studies/init_solver_benchmark.md) with
  concrete fix suggestions.
