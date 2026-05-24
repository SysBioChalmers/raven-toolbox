# ravengem — Refactoring Progress

Living status tracker for the RAVEN (MATLAB) → ravengem (Python/cobrapy) port.
**PLAN.md** is the design spec (what to port, why, and how it maps to cobrapy);
this file tracks **how far along** the port is. Update it whenever code lands.

> This is **not** a one-to-one transcription. When porting a function, also judge whether it can be
> made smarter/faster or whether RAVEN is missing something that fits — and log those in
> **[IMPROVEMENTS.md](IMPROVEMENTS.md)** (candidates to also back-port to MATLAB RAVEN).

_Last updated: 2026-05-24_

---

## Status at a glance

| Phase | Theme | Status |
|---|---|---|
| 0 | Scaffold & decisions | ✅ done |
| 1 | Foundation (`utils/`, `manipulation/`) | 🟢 functions done; package pip-installable (own `.venv`); CI remains |
| 2 | I/O (`io/`) | 🟢 done in scope — YAML/SIF/Excel-export/exportForGit (+sort_ids); Excel import excluded |
| 3a | Reconstruction — homology (`reconstruction/homology/`) | 🟢 implemented (get_model_from_homology + BLAST/DIAMOND wrappers); binary ZIPs/CI pending |
| 3b | Reconstruction — KEGG (`reconstruction/kegg/`) — 5-step pipeline: download → parse dump → build HMMs → model-for-species → model-by-HMM-query | 🟢 all 5 steps done (3b.1 download, 3b.2 dump parser, 3b.3 HMM libraries, 3b.4 species model, 3b.5 HMM query). `getPhylDist` distance-matrix deliberately not ported (fixed prok90/euk90 libs make it moot). |
| 3c | Reconstruction — MetaCyc | ❌ **dropped** (2026-05-24) — BLAST-to-single-representatives is low-precision at every cutoff; also to be removed from MATLAB RAVEN. See IMPROVEMENTS R-MetaCyc. |
| 4a | Metabolic tasks (`tasks/` — `parseTaskList`, `checkTasks`) — the task file | 🟢 done — `parse_task_list` + `check_tasks`; `fitTasks` + essential-rxn output deferred to 4c (tINIT consumer) |
| 4b | Gap-filling (`gapfilling/`) | 🟢 done — `connect_blocked_reactions` (connectivity, MILP via cobra/optlang); targeted mode → `cobra.gapfill` (cheatsheet) |
| 4c | tINIT (`init/` — original INIT MILP `getINITModel`/`runINIT` + scoring) | 🟡 in progress — `run_init` (the INIT MILP) done; `getINITModel` (expression scoring + task integration) pending |
| 4d | ftINIT (`init/` — fast staged INIT) — **⚠️ critical review of MATLAB code; most complex port** | ⬜ not started |
| 5 | Data integration & analysis (`omics/`, `analysis/`, `comparison/`) | ⬜ not started |
| 6 | Visualization (`plotting/`) | ⬜ not started |
| 7 | Localization (`localization/`) — `predictLocalization` + pluggable predictors (WoLF PSORT, DeepLoc, …); self-contained | ⬜ not started |

Legend: ✅ done · 🟡 in progress · ⬜ not started

---

## Ported files

Functions that exist as working, tested Python in ravengem.

| ravengem location | function(s) | RAVEN origin | tests | notes |
|---|---|---|---|---|
| `manipulation/irreversible.py` | `convert_to_irreversible` | `convertToIrrev.m` | ✅ `tests/test_manipulation_irreversible.py` | Adopted from geckopy `pipeline/preprocess.py`. Splits reversible non-exchange reactions into forward + `_REV`. cobra's old `convert_to_irreversible` was removed, so this is a real port. |
| `manipulation/expand.py` | `expand_model`, `_gpr_to_dnf`, `_node_to_dnf` | `expandModel.m` | ✅ `tests/test_manipulation_expand.py` | Adopted from geckopy `pipeline/expand.py`. Splits isozyme (OR-GPR) reactions into one reaction per AND-clause (`_EXP_N`), using cobra's GPR AST. |
| `utils/gpr.py` | `is_dnf`, `find_non_dnf_grrules`, `GPRIssue` | `standardizeGrRules.m` (`findPotentialErrors` half) | ✅ `tests/test_utils_gpr.py` | Lint half only — flags GPRs not in disjunctive normal form via cobra's GPR AST, structured `GPRIssue` output. The normalization half is **not** ported (cobra auto-normalizes GPRs on assignment). |
| `utils/parse.py` | `parse_name_comp` | `getIndexes.m` (`metcomps` sliver) | ✅ (via `test_manipulation_add.py`) | Parse a `name[comp]` token → `(name, compartment)`. The only cobra-absent bit of `getIndexes`. |
| `manipulation/add.py` | `add_reactions_from_equations` | `addRxns.m` | ✅ `tests/test_manipulation_add.py` | Keystone. Adds reactions from equation strings; matches mets by id, name, or `name[comp]`; assigns compartment to new mets; strict/auto policies for new mets & genes; duplicate-ID guard. Equation parsing/arrows/gene-creation delegated to cobra. |
| `manipulation/change.py` | `change_reaction_equations`, `change_gene_reaction_rules` | `changeRxns.m`, `changeGrRules.m` | ✅ `tests/test_manipulation_change.py`, `tests/test_change_grrules.py` | Stoichiometry change in place (reuses the `add` parser); and batch GPR set/append (`(old) or (new)`), with gene creation + normalization free from cobra. |
| `manipulation/parameters.py` | `set_variance_bounds` | `setParam.m` (`'var'` mode) | ✅ `tests/test_parameters.py` | The only `setParam` mode cobra has no idiom for: a ±% band around measured values (sign-aware). Other modes (`lb`/`ub`/`eq`/`obj`/`unc`) are cobra one-liners (cheatsheet). |
| `utils/balance.py` | `get_elemental_balance`, `ElementalBalance` | `getElementalBalance.m` | ✅ `tests/test_utils_balance.py` | Graded `balanced`/`unbalanced`/`unknown` per reaction — `unknown` catches a missing formula that cobra's `check_mass_balance` silently miscounts. |
| `utils/validate.py` | `check_model`, `ModelIssue` | `checkModelStruct.m` (curation subset) | ✅ `tests/test_utils_validate.py` | Structured curation report: orphan mets/genes, empty reactions, duplicate name+compartment, empty names, objective sanity. RAVEN's struct/type/duplicate-ID checks are moot in cobra. |
| `manipulation/transport.py` | `add_transport_reactions` | `addTransport.m` | ✅ `tests/test_manipulation_transport.py` | Transport reactions across compartments, matching mets by name, sequential `tr_0001` IDs, optional target-metabolite creation. cobra has no transport primitive. |
| `manipulation/transfer.py` | `add_reactions_from_model` | `addRxnsGenesMets.m` | ✅ `tests/test_manipulation_transfer.py` | Copy reactions from a source model, matching mets by `name[comp]` (not id), adding only new mets/genes. cobra's `merge` is strict-by-id. |
| `manipulation/merge.py` | `merge_models` | `mergeModels.m` | ✅ `tests/test_manipulation_merge.py` | Merge N models, unify mets by `name[comp]` (or id), keep all reactions (id collisions renamed), merge genes, provenance in `notes['origin']`. cobra's `merge` is pairwise/strict-by-id. |
| `manipulation/simplify.py` | `remove_dead_end_reactions`, `remove_duplicate_reactions`, `constrain_reversible_reactions`, `group_linear_reactions` | `simplifyModel.m` (gap modes) | ✅ `tests/test_manipulation_simplify.py` | The cobra-absent reduction modes; cobra-covered modes (no-flux→`find_blocked_reactions`, zero-interval, unconstrained) cheatsheeted. `group_linear` is lossy (drops genes), per RAVEN. |
| `io/sif.py` | `export_model_to_sif` | `exportModelToSIF.m` | ✅ `tests/test_io_sif.py` | Cytoscape SIF export (`rc`/`rr`/`cc` graphs). cobra has no network export. |
| `io/excel.py` | `export_to_excel` | `exportToExcelFormat.m` (export only) | ✅ `tests/test_io_excel.py` | RAVEN 5-sheet xlsx (RXNS/METS/COMPS/GENES/MODEL); RAVEN fields pulled from cobra annotation/notes. `openpyxl` (lazy, `[excel]` extra). Excel **import** intentionally excluded. cobra has no Excel I/O. |
| `io/git.py` | `export_for_git` | `exportForGit.m` | ✅ `tests/test_io_git.py` | Standard-GEM repo layout (`model/<fmt>/`) for yml/xml/mat/xlsx/txt + `dependencies.txt`; sorts a copy first. Orchestrates the other writers; cobra has no repo-layout writer. |
| `utils/sort.py` | `sort_identifiers` | `sortIdentifiers.m` | ✅ `tests/test_utils_sort.py` | Model-wide alphabetical `DictList.sort`; also `sort_ids=` on `write_yaml_model`. cobra has per-list sort only. |
| `io/yaml.py` | `read_yaml_model`, `write_yaml_model` | `readYAMLmodel.m` / `writeYAMLmodel.m` (RAVEN `fa281a1`) | ✅ `tests/test_io_yaml.py` | Aligned to RAVEN's cobra-native `!!omap` writer (`fa281a1`). cobra owns standard fields + the `annotation` block (smiles/ec-code/MIRIAM); this adds the RAVEN-only top-level per-entry keys (inchis/deltaG/metFrom/notes; confidence_score/references/rxnFrom/deltaG; protein) → `.notes`, plus `version`/`metaData`/GECKO `ec-*`. Output verified cobra-readable; legacy id-in-metaData supported. |
| `manipulation/remove.py` | `remove_metabolites`, `remove_genes` | `removeMets.m` / `removeGenes.m` | ✅ `tests/test_manipulation_remove.py` | Delegate to cobra; add the gaps: `by_name` cross-compartment deletion (mets — flagged as a deletion candidate if unused), and a `blocked_reactions` remove/constrain/keep policy for gene knockouts (genes). `removeReactions` **not** ported (coupled orphan cleanup = cobra's `remove_reactions`). |
| `reconstruction/homology/hits.py` | `make_ortholog_hits`, `HIT_COLUMNS` | `makeFakeBlastStructure.m` | ✅ `tests/test_reconstruction_homology.py` | Bidirectional hits DataFrame (the `blastStructure` replacement); the no-BLAST seeding/testing path. |
| `reconstruction/homology/homology.py` | `get_model_from_homology`, `HomologyResult` | `getModelFromHomology.m` | ✅ `tests/test_reconstruction_homology.py` | Core homology reconstruction with logic improvements H1–H6 (bidirectional/best_hits_only, AST GPR rewrite, complex_policy, bitscore best-hits, DataFrame ortholog map, provenance). No BLAST needed. |
| `reconstruction/homology/blast.py` | `run_blast`, `run_diamond`, `blast_from_table` | `getBlast.m`/`getDiamond.m`/`getBlastFromExcel.m` | ✅ `tests/test_reconstruction_blast.py` | Subprocess wrappers → hits DataFrame; `blast_from_table` loads a CSV (no Excel). `run_blast` verified against installed BLAST+. |
| `binaries.py` | `resolve_binary`, `ensure_binary` | `software/` provisioning | ✅ `tests/test_binaries.py` | Generic binary resolver (arg→env→PATH→bundled ZIP) + version-pinned release-ZIP registry (SHA256-verified cache). Registry empty until ZIPs published. |

**Test status:** 320 tests passing (incl. smoke) under cobra 0.31.1, run via ravengem's own `.venv` (`pip install -e '.[dev,excel]'`).

---

## Subpackage scaffold

All subpackages exist as importable stubs (purpose docstring only) unless noted above.

| subpackage | purpose | port status |
|---|---|---|
| `utils/` | GPR hygiene + balance + validation + parse helpers (`is_dnf`/`find_non_dnf_grrules` ✅, `get_elemental_balance` ✅, `check_model` ✅, `parse_name_comp` ✅) — **no** struct adapter; `getRxnsInComp`/`getMetsInComp`, MIRIAM/ID-prefix **not** ported (cobra covers) | 🟢 foundation done |
| `manipulation/` | model construction, editing & structural transforms (ergonomic layer, see PLAN §1b) | ✅ done — add/change/remove/transport/transfer/merge/simplify/variance + 2 adopted transforms |
| `io/` | RAVEN YAML/SIF/Excel-export/git-layout (Excel import excluded) | 🟢 done in scope |
| `reconstruction/homology/` | homology-based draft from a template GEM + BLAST/DIAMOND (3a) | ⬜ stub |
| `reconstruction/kegg/download.py` | Download + arrange KEGG FTP dump, pure stdlib (`fetch_keggdb.sh`, step 3b.1) | 🟢 done |
| `reconstruction/kegg/parse.py` | Parse KEGG dump → gene-free reference model + gzipped-TSV tables (`getRxnsFromKEGG`/`getMetsFromKEGG`/`getGenesFromKEGG`/`getModelFromKEGG`, step 3b.2) | 🟢 done |
| `reconstruction/kegg/organism.py` | Build draft model for a KEGG species from artefacts (`getKEGGModelForOrganism` no-FASTA branch, step 3b.4) | 🟢 done |
| `reconstruction/kegg/hmm.py` + `taxonomy.py` | Per-KO multi-FASTA + CD-HIT/MAFFT/hmmbuild → prok90/euk90 pressed HMM libraries (`constructMultiFasta` + train stages, step 3b.3) | 🟢 done |
| `reconstruction/kegg/query.py` + `assemble.py` | HMM-query de-novo path: hmmscan → assign_kos (cutoff + score ratios) → shared model assembler (`getKEGGModelForOrganism` FASTA branch, step 3b.5) | 🟢 done |
| `data.py` | `ensure_data` — fetch/verify/cache published artefacts (KEGG reference model + tables + HMM libs) under `~/.cache/ravengem/data/`, mirroring `ensure_binary`; auto-fetch wired into the `…_from_artefacts` entry points | 🟢 done (registry empty until published) |
| `gapfilling/fill.py` | `connect_blocked_reactions` | `fillGaps.m` (connectivity mode) | ✅ `tests/test_gapfilling.py` | MILP (min penalty-weighted template reactions s.t. blocked reactions carry flux) via cobra/optlang. No cobra equivalent. Templates matched by `name[comp]`. Targeted mode (fill toward objective) → `cobra.gapfill` (cheatsheet, §1). |
| `tasks/tasklist.py` + `check.py` | `parse_task_list` + `Task`; `check_tasks` + `TaskResult` | `parseTaskList.m` / `checkTasks.m` | ✅ `tests/test_tasks.py` | Task-list parser (TSV/xlsx; multi-row tasks, `;`-split, defaults, ALLMETS/ALLMETSIN) + feasibility checker that imposes inputs/outputs via relaxed metabolite mass-balance bounds (RAVEN's `b`), adds task equations, changes bounds, closes boundaries. `fitTasks`/essential-rxns deferred to 4c. |
| `init/init.py` | `run_init` + `InitResult` | `runINIT.m` | ✅ `tests/test_init.py` | The INIT MILP, clean optlang reformulation (split reversibles; `eps·x ≤ v ≤ ub·x` include-indicators; `x_fwd+x_rev≤1` for `no_rev_loops`; per-met reward sinks for `prod_weight`; LP-relaxation `present_mets` test). `getINITModel` scoring wrapper still pending (4c). |
| `init/` (rest) | `getINITModel` (4c scoring) + ftINIT (4d — critical review) | ⬜ pending |
| `omics/` | HPA omics → reaction scores | ⬜ stub |
| `localization/` | subcellular localization (Phase 7) — `predictLocalization` + pluggable predictors (WoLF PSORT, DeepLoc, …) | ⬜ stub |
| `analysis/` | reporterMetabolites, FSEOF, dFBA, … | ⬜ stub |
| `comparison/` | multi-model comparison | ⬜ stub |
| `plotting/` | pathway maps / omics overlay | ⬜ stub |

---

## Key decisions (summary — full rationale in [PLAN.md](PLAN.md))

1. **Name `ravengem`** — PyPI dist == import name. `ravenpy` was taken (a hydrology package).
2. **`cobra.Model` is the canonical object** — no `ravenCobraWrapper` adapter, no parallel RAVEN
   struct. RAVEN-only fields live in cobra `.annotation`/`.notes`.
3. **Do not re-port what cobrapy covers** (~70 RAVEN functions → cobra calls; documented as a
   migration cheatsheet, not code).
4. **YAML I/O** = cobrapy standard + geckopy ec extension keys (`ec-rxns`, `ec-enzymes`,
   `gecko_light`, `metaData`); plus read the legacy RAVEN/MATLAB `!!omap` dialect geckopy rejects.
5. **tINIT/tasks**: functional equivalence, not bit-exact.
6. **KEGG/MetaCyc data**: configurable — live REST (+disk cache) or reused RAVEN dumps.
7. **geckopy relocation strategy**: adopt copy now, converge later — ravengem holds the canonical
   copy; geckopy switches to `from ravengem.manipulation import ...` once ravengem is on PyPI.
   Until then, keep the two copies in sync.
8. **License**: GPL-3.0-or-later (derivative of GPLv3 RAVEN).

9. **Ergonomic layer (PLAN §1b)** — re-examined the RAVEN "manipulation" functions first dismissed
   as cobra-redundant. Many earn a port because they batch, chain steps, match by name, or
   auto-create dependencies cobra leaves to you. PORT keystones: `addRxns` (+ family), `setParam`,
   `removeReactions`/`removeMets`/`removeGenes`, `simplifyModel`, `mergeModels`,
   `sortModel` (deterministic core). These are a thin layer **on top of** cobra primitives, not a
   parallel data model. `standardizeGrRules` shrank to its lint half only — cobra auto-normalizes
   GPR syntax, so only `find_non_dnf_grrules`/`is_dnf` were ported (✅).
10. **`getIndexes` NOT ported** — cobra's `DictList` (`get_by_any`/`get_by_id`/`query`/`index`)
    already covers mixed id/object/index/name lookup more idiomatically; RAVEN needed a central
    resolver only because of its struct-of-parallel-arrays design. Kept just the `name[comp]`
    composite resolver as a small helper. Two findings (hash lookup, `[1 1 1]` mask/index bug)
    logged in IMPROVEMENTS.md as MATLAB-RAVEN-only back-port candidates.

### Open questions
- geckopy → ravengem dependency direction once ravengem is published (currently duplicated).

---

## Change log

Keyed to commits on `main`.

| commit | summary |
|---|---|
| `5b39f54` | Initial ravenpy scaffold + RAVEN→cobrapy port plan |
| `3a50e79` | Fix `.gitignore` shadowing the `reconstruction/metacyc` package |
| `1e24c00` | Record author decisions: name, repo home, INIT fidelity, DB sourcing |
| `000f031` | Drop RAVEN-COBRA adapter; scope YAML I/O to cobra standard + geckopy ec keys |
| `cb88c02` | Rename package to `ravengem` (PyPI dist == import name) |
| `ba71a90` | Plan to relocate RAVEN-derived code from geckopy; add `manipulation/` |
| `6798cb3` | Adopt `convert_to_irreversible` and `expand_model` from geckopy |
| `e6ce70d` | Reclassify ergonomic RAVEN functions as port targets; add PROGRESS.md |
| `30b1582` | Add IMPROVEMENTS.md; getIndexes improvement proposals |
| `62b43d1` | Demote getIndexes (cobra `DictList` covers it) |
| `4c65c8a` | Port GPR lint half of standardizeGrRules (`is_dnf`/`find_non_dnf_grrules`) |
| `2224eed` | Port addRxns as `add_reactions_from_equations`; add `parse_name_comp` |
| `5779f47` | Port changeRxns as `change_reaction_equations` |
| `5a4e292` | Port YAML I/O as `read_yaml_model`/`write_yaml_model` |
| `f869968` | Route RAVEN-only YAML fields to annotation/notes by meaning |
| `9169c25` | Drop `remove_reactions`; flag `remove_metabolites` by_name |
| `77f39cc` | Split reconstruction plan into homology/KEGG/MetaCyc tracks |
| `5972fed` | Remove compartment selectors (too thin over cobra) |
| `da80d9d` | Port changeGrRules (compartment selectors removed in `5972fed`) |
| `9dacc75` | Port setParam + getElementalBalance (3 simplest skipped) |
| `d0e63b6` | Trim set_parameters → set_variance_bounds |
| `fee1e6b` | Port addTransport as `add_transport_reactions` |
| `23d3ceb` | Skip setExchangeBounds (cobra `model.medium` covers it) |
| `6b557b2` | Port addRxnsGenesMets as `add_reactions_from_model` |
| `5d367b4` | Port checkModelStruct curation subset as `check_model` |
| `e6020c2` | Port mergeModels as `merge_models` |
| `cd2eea9` | Port simplifyModel gap modes (`manipulation/simplify.py`) |
| `a9c90cc` | Realign YAML I/O to RAVEN fa281a1 (cobra-native !!omap) |
| `53d94df` | Port exportModelToSIF as `export_model_to_sif` |
| `8879b57` | Plan exportForGit + sortIdentifiers |
| `5d4fef2` | Port exportToExcelFormat as `export_to_excel` (export only) |
| `df315ab` | Port sortIdentifiers + exportForGit (incl. xlsx) |
| `1b0e4a6` | Port KEGG dump parser (3b.2) → reference model + gzipped-TSV tables |
| `17f7eb3` | Port KEGG download/arrange (3b.1) as pure-stdlib tooling |
| `a61beb7` | Port KEGG model-for-species (3b.4) |
| `9d6f0af` | Port KEGG HMM-library construction (3b.3) |
| `217a367` | Port KEGG HMM-query path (3b.5) + domain mode; Phase 3b complete |
| `33d13c1` | Stream organism_gene_ko (real-data memory fix); 3b validated end-to-end on real KEGG dump |
| `16d1ef4` | Reference model as gzipped YAML; implement ensure_data artefact fetch/cache |
| `44ea6b5` | Add maintainer scripts: build KEGG artefacts + emit registry snippets |
| `b8d2764` | Calibrate MAFFT budget from measured memory curve; remove max_sequences cap |
| `4b9cae0` | Drop MetaCyc reconstruction (3c) — low-precision; split gap-filling into Phase 4b |
| `92deaa8` | Port gap-filling (Phase 4b): fill_gaps (connectivity) + gapfill_to_objective (targeted) |
| `c063daf` | Simplify gap-filling: drop targeted wrapper; rename to connect_blocked_reactions |
| `b66da06` | Restructure roadmap: split Phase 4 (4a/4c/4d); number localization as Phase 7 |
| `fdae0f5` | Port metabolic tasks (Phase 4a): parse_task_list + check_tasks |

---

## Next up

**Done:** Phases 1 (foundation), 2 (I/O), 3a (homology), 3b (KEGG), 4b (gap-filling).
**Dropped:** 3c (MetaCyc). Candidate next steps:

1. **Phase 4a** — metabolic `tasks/` (`parseTaskList`, `checkTasks`/`fitTasks`). The task
   file; foundation for the INIT phases.
2. **Phase 4c** — tINIT (original INIT MILP) — depends on 4a; needs a MIP solver.
3. **Phase 4d** — ftINIT (fast staged INIT). **⚠️ Critical, non-transcriptive review of the
   MATLAB code required — most complex port.** Depends on 4a, 4c.
4. **Phase 7** — Localization (`predictLocalization` + pluggable predictors: WoLF PSORT,
   DeepLoc, …). Self-contained (only needs Phase 1).
5. **Phase 5 / 6** — omics/analysis/comparison; visualization.
6. **CI** — GitHub Actions running the suite in ravengem's own venv.
