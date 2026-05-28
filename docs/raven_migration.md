# RAVEN → ravengem migration reference

A function-by-function map from the MATLAB RAVEN Toolbox to ravengem (and cobrapy where
appropriate). For each RAVEN function we record one of four outcomes:

* ✅ **ported** — there's a direct ravengem replacement; see the link.
* 🗒️ **cheatsheet** — cobrapy already covers it; a one-liner or short idiom does the job
  (recorded under each row).
* ⛔ **not ported** — explicit decision not to bring it across, with rationale.
* 🆕 **new in ravengem** — functionality ravengem adds that has no RAVEN counterpart.

For ravengem's deliberate improvements over RAVEN (and which of them are candidates to
upstream into MATLAB RAVEN), see [IMPROVEMENTS.md](../IMPROVEMENTS.md).

## Design principle

The in-memory object is always a [`cobra.Model`](https://cobrapy.readthedocs.io). There is
no parallel RAVEN struct, no `ravenCobraWrapper`-style adapter. RAVEN fields that cobra
doesn't model natively (`rxnMiriams`, `metDeltaG`, `rxnConfidenceScores`, …) live in
cobra's `annotation` / `notes` dictionaries. ravengem's YAML I/O follows the cobra YAML
standard plus the geckopy enzyme-constrained extension, so ecModels round-trip.

---

## Foundation: utilities, manipulation, I/O

| RAVEN | ravengem | Notes |
|---|---|---|
| `addRxns` (equations) | ✅ [`manipulation.add_reactions_from_equations`](../src/ravengem/manipulation/add.py) | The keystone: equation-string → reactions, matching mets by id, name, or `name[comp]`; strict / auto policies for new mets and genes. |
| `addTransport` | ✅ [`manipulation.add_transport_reactions`](../src/ravengem/manipulation/transport.py) | Transport across compartments, matching mets by name, sequential `tr_NNNN` ids. cobra has no transport primitive. |
| `addRxnsGenesMets` | ✅ [`manipulation.add_reactions_from_model`](../src/ravengem/manipulation/transfer.py) | Copy reactions from a source model, matching mets by `name[comp]` (vs cobra's strict-by-id merge). |
| `mergeModels` | ✅ [`manipulation.merge_models`](../src/ravengem/manipulation/merge.py) | N-model merge, unify mets by `name[comp]` (or id), keep all reactions (id collisions renamed). cobra's `merge` is pairwise / strict-by-id. |
| `changeRxns`, `changeGrRules` | ✅ [`manipulation.change_reaction_equations`](../src/ravengem/manipulation/change.py), `change_gene_reaction_rules` | Stoichiometry change in place; batch GPR set/append (`(old) or (new)`). |
| `setParam('var', …)` | ✅ [`manipulation.set_variance_bounds`](../src/ravengem/manipulation/parameters.py) | ±% band around measured values. **Other modes (lb/ub/eq/obj/unc)**: 🗒️ cobra one-liners (`reaction.bounds`, `model.objective`, `Configuration().bounds`). |
| `setExchangeBounds` | ⛔ not ported | cobra's `model.medium = {ex_id: uptake}` covers it. |
| `simplifyModel` (gap modes) | ✅ [`manipulation.remove_dead_end_reactions`](../src/ravengem/manipulation/simplify.py), `remove_duplicate_reactions`, `constrain_reversible_reactions`, `group_linear_reactions` | Cobra-covered modes (no-flux→`find_blocked_reactions`, zero-interval, unconstrained) are 🗒️ cheatsheeted. `group_linear` is lossy (drops genes), per RAVEN. |
| `removeMets`, `removeGenes` | ✅ [`manipulation.remove_metabolites`](../src/ravengem/manipulation/remove.py), `remove_genes` | Delegate to cobra; add `by_name` cross-compartment deletion (mets) and a `blocked_reactions` policy (genes). `removeReactions` itself ⛔ — cobra's `remove_reactions` covers it. |
| `convertToIrrev` | ✅ [`manipulation.convert_to_irreversible`](../src/ravengem/manipulation/irreversible.py) | Splits reversible non-exchange reactions into forward + `_REV`. Adopted from geckopy. |
| `expandModel` | ✅ [`manipulation.expand_model`](../src/ravengem/manipulation/expand.py) | Splits OR-GPR (isozyme) reactions into one per AND-clause. Adopted from geckopy. |
| `mergeCompartments` | ✅ [`manipulation.merge_compartments`](../src/ravengem/manipulation/compartments.py) | Collapse a multi-compartment model into one; deduplicate identical reactions; optionally drop one-met collapses (`drop_single_metabolite_reactions`). |
| `copyToComps` | ✅ [`manipulation.copy_to_compartment`](../src/ravengem/manipulation/compartments.py) | Duplicate reactions into a target compartment (idempotent; `delete_original=True` makes it a move). |
| `mapCompartments` | ⛔ not ported | Overlaps with `comparison.compare_models` on the reaction-id intersection. |
| `getElementalBalance` | ✅ [`utils.get_elemental_balance`](../src/ravengem/utils/balance.py) | Graded `balanced` / `unbalanced` / `unknown` — `unknown` catches a missing formula that cobra's `check_mass_balance` silently miscounts. |
| `checkModelStruct` (curation subset) | ✅ [`utils.check_model`](../src/ravengem/utils/validate.py) | Structured curation report. RAVEN's struct/type checks are moot in cobra. |
| `is_dnf` / GPR check (from `standardizeGrRules`) | ✅ [`utils.is_dnf`](../src/ravengem/utils/gpr.py), `find_non_dnf_grrules` | Lint-only half; cobra auto-normalises GPRs on assignment, so the rewriting half isn't ported. |
| `getIndexes` (`metcomps` sliver) | ✅ [`utils.parse_name_comp`](../src/ravengem/utils/parse.py) | The only `getIndexes` bit cobra doesn't already cover. |
| `editMiriam`, `extractMiriam` | ⛔ not ported | cobra's `.annotation` is already a `{namespace: id(s)}` dict — read/write it directly. |
| `getRxnsInComp`, `getMetsInComp` | ⛔ not ported | One-liners over cobra's `reaction.compartments` / `metabolite.compartment`. |
| `constructEquations` | ⛔ not ported | `reaction.build_reaction_string(use_metabolite_names=...)` already does both id and name equations. |
| `sortIdentifiers` | ✅ [`utils.sort_identifiers`](../src/ravengem/utils/sort.py) | Model-wide alphabetical sort; also via `sort_ids=` on `write_yaml_model`. |

### I/O

| RAVEN | ravengem | Notes |
|---|---|---|
| `readYAMLmodel`, `writeYAMLmodel` | ✅ [`io.read_yaml_model`, `write_yaml_model`](../src/ravengem/io/yaml.py) | Aligned to cobra's `!!omap` writer (RAVEN `fa281a1`). Adds the RAVEN-only top-level per-entry keys (inchis/deltaG/metFrom/notes, confidence_score/references/rxnFrom/deltaG, protein) into `.notes`, plus `version`/`metaData`/GECKO `ec-*`. cobra-readable output verified. |
| `exportModelToSIF` | ✅ [`io.export_model_to_sif`](../src/ravengem/io/sif.py) | Cytoscape SIF (`rc`/`rr`/`cc` graphs). |
| `exportToExcelFormat` (export only) | ✅ [`io.export_to_excel`](../src/ravengem/io/excel.py) | RAVEN 5-sheet xlsx (RXNS / METS / COMPS / GENES / MODEL). Excel **import** is intentionally excluded. |
| `exportForGit` | ✅ [`io.export_for_git`](../src/ravengem/io/git.py) | Standard-GEM repo layout (`model/<fmt>/…`). |
| `importYAML/SBML/Mat/Excel` | 🗒️ cobra's standard readers | `cobra.io.read_sbml_model` / `load_json_model` / etc.; Excel import not ported. |

## Reconstruction

| RAVEN | ravengem | Notes |
|---|---|---|
| `getModelFromHomology` + `getBlast`/`getDiamond` | ✅ [`reconstruction.homology.get_model_from_homology`](../src/ravengem/reconstruction/homology/homology.py), `run_blast`, `run_diamond`, `blast_from_table` | Core homology reconstruction with structured improvements (bidirectional / best-hits-only, AST GPR rewrite, complex policy, bitscore best-hits, DataFrame ortholog map). |
| KEGG download → species model (5 steps) | ✅ [`reconstruction.kegg.*`](../src/ravengem/reconstruction/kegg/) | All five steps: `download.fetch_keggdb`, `parse.read_kegg_table` + reference model, `hmm.build_libraries`, `organism.build_kegg_model_for_organism` (no-FASTA), `query.assign_kos` + `run_hmmscan`. |
| `getPhylDist` | ⛔ not ported | Fixed prok90/euk90 libraries make the distance matrix moot. |
| `getMetaCycModelForOrganism` | ⛔ not ported (and **flagged for removal from MATLAB RAVEN**) | BLAST-to-single-representatives is low-precision at every cutoff. See [IMPROVEMENTS.md](../IMPROVEMENTS.md) under `R-MetaCyc`. |

## Tasks, gap-filling, INIT, ftINIT

| RAVEN | ravengem | Notes |
|---|---|---|
| `parseTaskList`, `checkTasks` | ✅ [`tasks.parse_task_list`](../src/ravengem/tasks/tasklist.py), `check_tasks` ([tasks/check.py](../src/ravengem/tasks/check.py)) | `check_tasks` reuses one model across the task list (no per-task model copy) — at genome scale ~12× faster than the copy-per-task implementation it replaced. |
| `fillGaps` (connectivity mode) | ✅ [`gapfilling.connect_blocked_reactions`](../src/ravengem/gapfilling/fill.py) | MILP (min penalty-weighted template reactions s.t. blocked reactions carry flux). Targeted mode → 🗒️ `cobra.gapfill`. |
| `runINIT`, `scoreComplexModel`, `getINITModel` | ✅ [`init.run_init`](../src/ravengem/init/init.py), [`init.score_reactions_from_genes`, `gene_scores_from_expression`](../src/ravengem/init/score.py), [`init.get_init_model`](../src/ravengem/init/build.py) | INIT MILP rewritten in optlang (no sparse-matrix construction); RNA-seq scoring is `5·ln(level/ref)`-clamped. |
| `ftINIT`, `prepINITModel`, `ftINITInternalAlg`, `getINITSteps` | ✅ [`init.ftinit`](../src/ravengem/init/ftinit.py), [`init.prep_init_model`](../src/ravengem/init/prep.py), [`init.run_ftinit`](../src/ravengem/init/ftinit.py), [`init.get_init_steps`](../src/ravengem/init/steps.py) | Staged MILP + linear merge + scaling (`rescaleModelForINIT`). Validated against RAVEN on Human-GEM (Jaccard 0.975–0.980; see [humangem_validation.md](humangem_validation.md)). Metabolomics-based scoring is the one piece not yet implemented (raises `NotImplementedError`). |
| `ftINITFillGaps`, `ftINITFillGapsMILP`, `ftINITFillGapsForAllTasks` | ✅ [`init.fill_tasks`](../src/ravengem/init/taskfill.py) | Task-aware gap-filling within ftINIT; in-place `_feasible` check + bounded fill MILP (`mip_gap`, `time_limit`). |
| `mergeLinear` | ✅ [`init.merge_linear`](../src/ravengem/init/merge.py) | Linear merge of unit-stoichiometry chains; bookkeeping (`group_ids`, `reversed_rxns`) to map back to the reference model. |
| `removeLowScoreGenes` | ✅ [`init.remove_low_score_genes`](../src/ravengem/init/genes.py) | Final gene-prune step of ftINIT. |
| `fitTasks` | ⛔ not ported | Niche; tasks are normally consumed via `checkTasks` / ftINIT's task layer. |

## Omics, analysis, comparison

| RAVEN | ravengem | Notes |
|---|---|---|
| `parseHPA`, `parseHPArna`, `scoreModel` | ✅ [`omics.parse_hpa`, `parse_hpa_rna`, `hpa_gene_scores`, `rna_gene_scores`](../src/ravengem/omics/hpa.py) | Pandas-tidy DataFrames; scoring adapters reuse `score_reactions_from_genes` (single source of truth for the GPR walk). |
| `reporterMetabolites` | ✅ [`analysis.reporter_metabolites`](../src/ravengem/analysis/reporter.py) | Exact closed-form background replaces RAVEN's Monte-Carlo (RM1 in IMPROVEMENTS). |
| `FSEOF` | ✅ [`analysis.fseof`](../src/ravengem/analysis/fseof.py) | Regression slope + correlation, amplify/knockdown/knockout classes, gene aggregation (FS1–FS4 in IMPROVEMENTS). |
| `randomSampling` | ✅ [`analysis.sample`](../src/ravengem/analysis/sampling.py) | Wraps cobra's flux sampling. |
| `compareMultipleModels` | ✅ [`comparison.compare_models`](../src/ravengem/comparison/compare.py) | Tidy DataFrames (reactions / mets / genes / subsystems presence + pairwise Jaccard + optional `check_tasks` pass/fail). Plotting and tSNE/MDS are 🗒️ one-liners in seaborn/scikit-learn; intentionally not in the function. |
| `runDynamicFBA` | ⛔ not ported | Established Python implementations exist: [`dfba`](https://pypi.org/project/dfba/) (Pinheiro et al.; CVODES-backed), [`reframed`](https://pypi.org/project/reframed/) (Machado lab), [`mewpy`](https://pypi.org/project/mewpy/) (Cunha lab). Cobrapy itself has none, but re-porting would duplicate maintained prior art. |

## Localisation (Phase 7)

| RAVEN | ravengem | Notes |
|---|---|---|
| `predictLocalization` | ✅ [`localization.predict_localization`](../src/ravengem/localization/predict.py) | Deterministic MILP (not simulated annealing). Caller-passed `reactions_to_relocate` set (everything else pinned). Multi-compartment by default: primary "free", extras pay `multi_compartment_penalty`. Tolerates incomplete models (no silent reaction removal). `apply=False` returns a `LocalizationProposal` diff. Real-data validation against curated yeast-GEM in [yeast_localization_benchmark.md](yeast_localization_benchmark.md). |
| `getWoLFScores`, `parseScores('wolf')` | ✅ [`localization.load_wolfpsort`](../src/ravengem/localization/scores.py) | Parses WoLF PSORT summary output (RAVEN-compatible); row-normalised. Does not shell out to the WoLF PSORT binary — run that separately and feed in the output. |
| `parseScores('deeploc')` | ✅ [`localization.load_deeploc`](../src/ravengem/localization/scores.py) | DeepLoc 2 per-protein CSV (Protein_ID / Localizations / Signals + one column per compartment). |

## Things deliberately not ported

* **`ravenCobraWrapper` / RAVEN struct adapter** — cobra is the canonical object; no parallel struct.
* **`checkModelStruct` struct/type checks** — moot in cobra.
* **`runDynamicFBA`** — see Omics/analysis row.
* **`getMetaCycModelForOrganism`** — see Reconstruction row; flagged for upstream removal.
* **`getPhylDist`** — fixed prok90/euk90 libraries make it moot.
* **`mapCompartments`** — overlaps with `compare_models`.
* **`editMiriam`, `extractMiriam`, `getRxnsInComp`, `getMetsInComp`, `constructEquations`, `getIndexes`** (most), **`setExchangeBounds`**, **most `setParam` modes**, **`getBlastFromExcel` Excel branch** — cobra one-liners; recorded above.

## "New in ravengem" entry points

These are not direct RAVEN ports:

* [`comparison.compare_models`](../src/ravengem/comparison/compare.py) — already in RAVEN, but the ravengem version returns tidy DataFrames suitable for downstream analysis (RAVEN's version embeds heatmap/tSNE plotting).
* [`reconstruction.homology.run_blast` / `run_diamond` / `blast_from_table`](../src/ravengem/reconstruction/homology/blast.py) — generic subprocess wrappers; a DataFrame is the canonical "hits" object (no `blastStructure` struct).
* [`binaries.resolve_binary`, `ensure_binary`](../src/ravengem/binaries.py) — version-pinned release-ZIP registry for external tools (BLAST/DIAMOND/HMMER), SHA256-verified.
* [`localization.LocalizationProposal`](../src/ravengem/localization/predict.py) — the diff-preview mode (`apply=False`).
