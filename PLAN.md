# ravengem — Port Plan

A plan for porting the [RAVEN Toolbox 2](https://github.com/SysBioChalmers/RAVEN)
(MATLAB, v2.11.1, ~271 `.m` files) to a Python package built on **cobrapy**.

The guiding principle is **do not re-port what cobrapy already provides.** RAVEN and cobrapy
overlap heavily on simulation and model-manipulation primitives. The value of `ravengem` is in
RAVEN's *reconstruction* and *context-specific modeling* features, which cobrapy lacks. This
document (1) maps every RAVEN functional area to its cobrapy equivalent or "port" verdict, and
(2) lays out a phased roadmap.

---

## 0. Foundational design decisions

| Decision | Choice | Rationale |
|---|---|---|
| In-memory model object | `cobra.Model` (not a re-implemented RAVEN struct) | Reuse cobrapy's data model, solver interface, SBML I/O, ecosystem interop. |
| RAVEN-only fields (`rxnMiriams`, `metDeltaG`, `rxnConfidenceScores`, `compMiriams`, `geneShortNames`, `eccodes`, `inchis`, …) | Stored in `cobra` `.annotation` / `.notes`; `subSystems` → SBML groups | cobrapy already round-trips these through SBML; no parallel struct needed. |
| RAVEN↔COBRA conversion | **None.** ravengem functions consume and produce `cobra.Model` directly — there is no `ravenCobraWrapper` port and no parallel RAVEN struct | Adhere to the cobrapy model format throughout; full COBRA-ecosystem interop for free, nothing to keep in sync. |
| Solver | cobrapy's optlang (GLPK/Gurobi/CPLEX/HiGHS) | Replaces RAVEN's `solveLP`/`solveQP`/`optimizeProb` MILP/LP layer. MILP for INIT/gap-filling needs a MIP-capable solver (Gurobi/CPLEX/SCIP). |
| Package name | `ravengem` (PyPI distribution **and** import name) | "raven" keeps discoverability; "gem" = genome-scale metabolic model. The obvious `ravenpy` is taken on PyPI by an unrelated hydrology package, so it was avoided; `ravengem` is free and dist==import. |
| Repo home | `SysBioChalmers/ravengem` | Alongside the MATLAB RAVEN. |
| tINIT/ftINIT & tasks fidelity | **Functional equivalence** (not bit-exact) | Same algorithm/intent; minor numerical differences from solver behavior are acceptable. |
| KEGG/MetaCyc data source | **Configurable**: live REST (+disk cache) *or* reuse RAVEN's pre-built dumps | Live for currency; dumps for reproducibility with MATLAB results. |
| License | GPL-3.0-or-later | Derivative of GPLv3 RAVEN. |

---

## 1. Already in cobrapy — DO NOT PORT (use/wrap instead)

These RAVEN functions have direct cobrapy equivalents. `ravengem` should provide, at most, thin
convenience wrappers or documentation mapping the old names.

> **Caveat (see §1b):** "cobra has an equivalent" is not the same as "cobra makes it as
> easy." A first pass put many model-*manipulation* functions here; a closer read of the MATLAB
> source shows several do batch input, multi-step orchestration, name-based matching, or
> auto-creation of dependencies that cobra forces you to write by hand. Those have been **moved
> to §1b as genuine port targets.** The table below is now restricted to functions that really do
> collapse to a one-liner with no ergonomic loss (simulation, analysis, SBML I/O, trivial lookups).

| RAVEN function(s) | cobrapy equivalent |
|---|---|
| `importModel`, `exportModel` (SBML) | `cobra.io.read_sbml_model`, `write_sbml_model` |
| `solveLP`, `optimizeProb`, `checkSolution` | `model.optimize()`, `model.slim_optimize()` |
| `solveQP`, `qMOMA` | `cobra.flux_analysis.moma` / `room` |
| `getAllowedBounds` (FVA) | `cobra.flux_analysis.flux_variability_analysis` |
| `getMinNrFluxes` (parsimonious) | `cobra.flux_analysis.pfba` |
| `findGeneDeletions` | `single_gene_deletion`, `double_gene_deletion`, `single/double_reaction_deletion` |
| `getEssentialRxns` | `cobra.flux_analysis.find_essential_reactions` / `find_essential_genes` |
| `randomSampling`, `analyzeSampling` | `cobra.sampling` (OptGP, ACHR) |
| `runProductionEnvelope`, `runPhenotypePhasePlane` | `cobra.flux_analysis.production_envelope` |
| `getExchangeRxns`, `getTransportRxns` | `model.exchanges`, `model.boundary`, custom filters |
| `buildEquation`, `parseRxnEqu` | `reaction.reaction`, `reaction.build_reaction_from_string` |
| `constructS` | `cobra.util.create_stoichiometric_matrix` |
| `printFluxes`, `printModel`, `printModelStats` | `model.summary()`, `reaction.summary()`, `metabolite.summary()` |
| `getGenesFromGrRules` | `cobra.core.gene.GPR` (parse/eval/AST) |
| `getIndexes` | `DictList.get_by_any` (mixed id/object/index→objects), `get_by_id`, `query` (name/regex), `index` (position). See §1b note — only the `name[comp]` resolver sliver is kept. |
| `parseFormulas` | `metabolite.formula` / `elements` |
| `deleteUnusedGenes` | `cobra.manipulation.prune_unused_metabolites`, `prune_unused_reactions` |
| `sortIdentifiers`, `permuteModel` | native Python indexing / cobra DictList |

**Verdict:** these collapse into cobrapy calls with no ergonomic loss. Capture them as a "migration
cheatsheet" in the docs rather than as code. (Functions reclassified out of this list — `addRxns`,
`removeReactions`, `setParam`, `getIndexes`, `mergeModels`, `simplifyModel`, etc. — are in §1b.)

---

## 1b. cobra has it, but RAVEN's version is worth porting (ergonomic layer)

Re-examined the MATLAB source of every "manipulation" function first parked in §1. The verdict
rubric: a function earns a **PORT** if it (1) batches over many objects in one call, (2) chains
multiple steps, (3) does matching / validation / auto-creation / cascading cleanup cobra leaves to
you, or (4) lets you work with human-readable equation strings. **WRAP** = minor convenience, often
just there to support a PORT. These form a thin RAVEN-style ergonomic layer over `cobra.Model`;
each is still implemented *on top of* cobra primitives, not as a parallel data model.

### Targets for `manipulation/` — model construction & editing

| RAVEN | verdict | Why (beyond the cobra one-liner) |
|---|---|---|
| `addRxns` | **PORT** (keystone) | Adds a batch of reactions from **equation strings** *or* mets+coeffs; auto-creates missing metabolites (`allowNewMets`, 3 matching modes: by id, by `name+comp`, or `name[comp]` syntax) and missing genes from grRules (`allowNewGenes`); validates GPR genes. cobra needs hand-built objects + coeff dicts + pre-created genes. |
| `addRxnsGenesMets` | **PORT** | Copies a batch of reactions from a source model into a draft, matching mets by `name[comp]` (not id), skipping/​reporting duplicates, auto-adding only genuinely new mets/genes with annotation carried over. The post-homology reaction-transfer workflow; cobra makes you hand-write the merge+dedup. |
| `addTransport` | **PORT** | Batch-creates transport reactions from one compartment to many, matching mets **by name** across comps, auto-naming (`tr_0001`), and auto-creating the target-compartment metabolite when missing. cobra has **no** transport primitive at all. |
| `changeRxns` | **PORT** (cheap once `addRxns` exists) | Replace a batch of reactions' stoichiometry via **equation strings**, preserving all other fields and original ordering. cobra has no string-equation edit path. |
| `changeGrRules` | **PORT** | Batch-set grRules with an **append** mode (`(old) or (new)`), auto-adding new genes, then re-standardizing and rebuilding `rxnGeneMat`. cobra's `gene_reaction_rule=` is per-reaction, no append, no consistent gene auto-creation. |
| `setParam` | **PORT** | One call sets `lb`/`ub`/`eq`/`obj`/`var`(±% band)/`unc` over a list of reactions (ids/index/mask), broadcasts a scalar, silently skips missing reactions, resets objective on `obj`, validates `lb≤ub`. cobra scatters these across attributes + manual loops. (Drop RAVEN-only `rev`.) |
| `setExchangeBounds` | **PORT** | Real media-definition logic: finds exchanges, maps mets by name/id/index, auto-detects import direction, refuses inconsistent `closeOthers`, optionally restricts to the extracellular compartment, can close all other imports. Richer than `model.medium`. |
| `removeReactions` | **PORT** | Three **separable** cascade flags (`removeUnusedMets`/`Genes`/`Comps`) vs cobra's single coupled `remove_orphans`; accepts ids / mask / index interchangeably. |
| `removeMets` | **PORT** | Delete mets **by name across all compartments** at once (`isNames`), then cascade to orphaned reactions/genes/comps and remap `metComps`. cobra's `remove_metabolites` has only `destructive`. |
| `removeGenes` | **PORT** | Flux-aware: rewrites GPRs dropping the gene, and a `removeBlockedRxns` toggle to either delete reactions that can no longer carry flux **or** constrain them to 0. (Reimplement GPR eval via cobra's GPR AST, not MATLAB `eval`.) |
| `simplifyModel` | **PORT** (stage by mode) | Orchestrator with 8 reduction modes + reserved-reaction protection + a deletion audit log: delete unconstrained / duplicate (`contractModel`) / zero-interval / inaccessible (dead-end) / no-flux (FVA) reactions, `groupLinear` enzyme chains, `constrainReversible`. Modes 1/3/4 are pure-graph (easy); 5 needs FVA; 6 is complex+lossy. No bundled cobra equivalent. |
| `mergeModels` | **PORT** | Merge **N** models at once, matching mets by `name+comp` (or id), auto-renaming id conflicts, reconciling compartments, tracking `rxnFrom`/`metFrom`/`geneFrom` provenance. cobra's `merge` is pairwise and strict-by-id. |
| `sortModel` (core) | **PORT** core / **SKIP** optimizer | Port the **deterministic** canonical ordering (mets by `name[comp]`; reversible reactions flipped to lexicographic-first reactant) — exactly what makes YAML diffable. **Skip** the stochastic `sortReactionOrder` annealer. |
| `addMets` | **WRAP** | Batch add with dedupe + `copyInfo` (copy formula/charge/InChI/MIRIAM from same-named met in another comp). Mostly exists to back `addRxns`. |
| `addExchangeRxns` | **WRAP** | Batch over a met list + RAVEN naming (`EXC_OUT_<id>`); `model.add_boundary` already covers the rest. |
| `constructEquations` | **WRAP** | Inverse of `addRxns`: batch reactions → readable equation strings (met id/name/formula, sorting). Backs `addRxnsGenesMets`; cobra has per-reaction `.reaction`. |

### Targets for `utils/` — lookup, GPR hygiene, balance

| RAVEN | verdict | Why |
|---|---|---|
| `getIndexes` | **DO NOT PORT** (keep `name[comp]` sliver only) | RAVEN needs a central index resolver because it's a struct of parallel arrays. cobra is object-oriented and already covers mixed lookup more idiomatically: `DictList.get_by_any` (mixed id/object/**index** → objects), `get_by_id` (O(1)), `query` (name/substring/regex), `index` (position), comprehensions for filtering. A 1-based-index port would be redundant and un-Pythonic. **Keep only** the `name[comp]` composite resolver (`parse_name_comp`) as a small helper for `addRxns`/`addTransport`/`mergeModels` — that's the one bit cobra lacks. |
| `standardizeGrRules` | **PORT lint only** ✅ | Two halves: (1) syntax normalization — **not ported**, cobra auto-normalizes every GPR on assignment (case, whitespace, redundant brackets); (2) the `findPotentialErrors` lint for non-DNF rules — **ported** as `find_non_dnf_grrules` + `is_dnf` ([utils/gpr.py](src/ravengem/utils/gpr.py)), reworked onto cobra's GPR AST and returning structured `GPRIssue`s instead of printing. |
| `getElementalBalance` | **WRAP** | Batch graded balance table (status: balanced / unbalanced / missing-info / parse-error) with InChI fallback — more informative than per-reaction `check_mass_balance`. |
| `getRxnsInComp`, `getMetsInComp` | **WRAP** | "Objects in compartment" accessors cobra lacks as first-class calls; the only non-trivial bit is `getRxnsInComp(include_partial=False)` (fully-contained vs touching). |

---

## 2. RAVEN-UNIQUE — the real port targets

Organized by the `src/ravengem/` subpackage that will host them.

### 2.1 `utils/` — model helpers (NO struct adapter)  *(Phase 1, foundational)*
The `ravenCobraWrapper` adapter is **explicitly not ported** — ravengem works on `cobra.Model`
directly. Only the genuinely useful, RAVEN-flavored helpers survive, as small functions on cobra
objects:
| RAVEN | Notes |
|---|---|
| `checkModelStruct` | Validate a `cobra.Model` against RAVEN reconstruction expectations (beyond cobra's own `model.validate` / SBML validation). |
| `editMiriam`, `extractMiriam` | Convenience get/set for MIRIAM-style entries inside cobra `.annotation` dicts. |
| `addIdentifierPrefix`, `removeIdentifierPrefix` | `R_`/`M_`/`G_` prefix handling for interop with COBRA-Toolbox-style IDs (only where cobra's SBML layer doesn't already cover it). |
| `parse_name_comp` (from `getIndexes` `metcomps`) | **PORT** (small helper) — resolve the `name[comp]` composite (metabolite by name + compartment) cobra can't. The rest of `getIndexes` is **not ported** — `DictList.get_by_any`/`get_by_id`/`query`/`index` already cover mixed id/object/index/name lookup more idiomatically (see §1b note). |
| `standardizeGrRules` | **PORT lint only** ✅ — normalization is cobra-automatic; non-DNF lint ported as `find_non_dnf_grrules`/`is_dnf`. See §1b. |
| `getElementalBalance` | **WRAP** — batch graded mass-balance table with InChI fallback; see §1b. |
| `getRxnsInComp`, `getMetsInComp` | **WRAP** — "objects in compartment" accessors (`include_partial` containment logic); see §1b. |
| ~~`ravenCobraWrapper`, `standardizeModelFieldOrder`, `cobraNamespaces.csv`, `COBRA_structure_fields.csv`~~ | **Dropped** — no parallel struct to convert or order. |

### 2.1b `manipulation/` — model construction, editing & structural transforms  *(Phase 1)*
The home for the RAVEN ergonomic layer (§1b) that *mutates* a `cobra.Model`, plus structural
transforms cobra lacks. Two transforms are **already ported in geckopy** and relocated here (see §7).

**Structural transforms:**
| RAVEN | Notes |
|---|---|
| `convertToIrrev` | Split reversible non-exchange reactions into forward + `_REV` pair. **Already ported** as `convert_to_irreversible` (geckopy `pipeline/preprocess.py`); cobra's old `convert_to_irreversible` was removed, so this is a real port, not a wrapper. |
| `expandModel` | Split isozyme (OR-GPR) reactions into one reaction per AND-clause (`_EXP_N`). **Already ported** as `expand_model` + `_gpr_to_dnf`/`_node_to_dnf` (geckopy `pipeline/expand.py`), using cobra's GPR AST instead of RAVEN string manipulation. |
| `simplifyModel` | **PORT, stage by mode** — 8 reduction modes + reserved-rxn protection + deletion log; see §1b. |
| `mergeModels` | **PORT** — N-way merge with name+comp matching, conflict rename, provenance; see §1b. |
| `sortModel` | **PORT** deterministic canonical ordering (skip stochastic optimizer); see §1b. |
| `contractModel`, `mergeCompartments`, `copyToComps` | Other RAVEN structural transforms to port as needed. |

**Construction & editing (ergonomic layer, §1b):** `addRxns` (keystone — equation-string batch add
with met/gene auto-creation), `addRxnsGenesMets`, `addTransport`, `changeRxns`, `changeGrRules`,
`setParam`, `setExchangeBounds`, `removeReactions`, `removeMets`, `removeGenes`; plus `addMets`,
`addExchangeRxns`, `constructEquations` as supporting **WRAP**s. See §1b for per-function rationale
and verdicts. **Build order:** `parse_name_comp` (utils) → `addMets`/
`constructEquations` → `addRxns` → everything that depends on it (`changeRxns`, `addRxnsGenesMets`,
`addTransport`).

### 2.2 `io/` — RAVEN-specific formats  *(Phase 1–2)*
| RAVEN | Notes |
|---|---|
| `readYAMLmodel`, `writeYAMLmodel` | **Follow the cobrapy YAML standard** (`cobra.io.dict.model_to_dict` / `model_from_dict` for the metabolites/reactions/genes/compartments portion), and **additionally support geckopy's enzyme-constrained extension keys** (`ec-rxns`, `ec-enzymes`, `gecko_light`, `metaData`) so ecModels round-trip. Empty/NaN fields omitted; loader fills them back. **Also read the legacy RAVEN/MATLAB YAML dialect** (`!!omap` tags) that geckopy's loader rejects — this is a ravengem-unique capability. `ruamel.yaml`. High priority — heavily used. |
| `importExcelModel`, `exportToExcelFormat`, `SBMLFromExcel` | RAVEN Excel model format. `openpyxl` (optional dep). |
| `exportToTabDelimited`, `exportForGit` | Plain-text / git-friendly model dumps. |
| `exportModelToSIF` | Cytoscape SIF export (visualization). |
| `getToolboxVersion`, `getMD5Hash` | Provenance helpers. |

### 2.3 `reconstruction/` — de novo reconstruction  *(Phase 3, flagship)*
The single biggest reason RAVEN exists and cobrapy does not cover this at all.

**`reconstruction/kegg/`**
| RAVEN | Notes |
|---|---|
| `getKEGGModelForOrganism` | Top-level: build a draft GEM for an organism from KEGG. Orchestrates the below. |
| `getModelFromKEGG`, `getRxnsFromKEGG`, `getMetsFromKEGG`, `getGenesFromKEGG` | Parse KEGG flat-files / REST into a model. |
| `getPhylDist` | Phylogenetic-distance weighting of KEGG orthologs. |
| `constructMultiFasta` | Build per-KO FASTA for homology search. |

**`reconstruction/metacyc/`**
| RAVEN | Notes |
|---|---|
| `getMetaCycModelForOrganism`, `getModelFromMetaCyc` | MetaCyc-based draft reconstruction. |
| `getRxnsFromMetaCyc`, `getMetsFromMetaCyc`, `getEnzymesFromMetaCyc` | MetaCyc flat-file parsers. |
| `linkMetaCycKEGGRxns`, `combineMetaCycKEGGModels`, `addSpontaneousRxns` | Cross-DB reconciliation. |

**`reconstruction/homology/`**
| RAVEN | Notes |
|---|---|
| `getModelFromHomology` | Transfer reactions from a template GEM via bidirectional best hits. |
| `getBlast`, `getDiamond`, `getBlastFromExcel` | Wrap external BLAST+/DIAMOND executables (subprocess). |
| `makeFakeBlastStructure`, `parseScores` | Homology-score plumbing. |

### 2.4 `init/` — context-specific models (tINIT / ftINIT)  *(Phase 4, flagship)*
RAVEN-unique MILP algorithm; no cobrapy equivalent. Needs a MIP solver.
| RAVEN | Notes |
|---|---|
| `ftINIT`, `getINITModel`, `runINIT` | Top-level extraction of a context model from omics scores. |
| `prepINITModel`, `getINITSteps`, `ftINITInternalAlg`, `INITStepDesc` | ftINIT staged algorithm. |
| `ftINITFillGaps`, `ftINITFillGapsMILP`, `ftINITFillGapsForAllTasks` | Task-aware gap-filling within INIT. |
| `scoreComplexModel`, `getExprForRxnScore`, `groupRxnScores`, `removeLowScoreGenes` | Reaction scoring from gene expression. |
| `mergeLinear`, `rescaleModelForINIT`, `reverseRxns` | MILP preprocessing. |

### 2.5 `tasks/` — metabolic task validation  *(Phase 4, flagship)*
No cobrapy equivalent.
| RAVEN | Notes |
|---|---|
| `checkTasks`, `fitTasks` | Run a task list (required/forbidden production) against a model. |
| `parseTaskList` | Parse the Excel/text task definition format. |
| `checkProduction`, `getExpressionStructure` | Production checks underpinning tasks. |

### 2.6 `gapfilling/` — RAVEN gap-filling  *(Phase 4)*
RAVEN's template-based MILP gap-filling differs from cobrapy's `gapfilling.GapFiller`.
| RAVEN | Notes |
|---|---|
| `fillGaps` | MILP fill from a reference/template model. |
| `gapReport` | Connectivity / dead-end / blocked report. |
| `checkProduction`, `canProduce`, `canConsume`, `makeSomething`, `consumeSomething` | Production/consumption diagnostics. |
| `removeBadRxns` | Remove reactions enabling erroneous production. |
| `fitParameters` | Parameter fitting. |

### 2.7 `omics/` & `localization/` — data integration  *(Phase 5)*
| RAVEN | Notes |
|---|---|
| `parseHPA`, `parseHPArna`, `scoreModel` | Human Protein Atlas → reaction scores (feeds INIT). |
| `predictLocalization`, `getWoLFScores` | WoLF PSORT-based subcellular localization → compartmentalize model. |
| `mapCompartments`, `mergeCompartments`, `copyToComps`, `getMetsInComp`, `getRxnsInComp` | Compartment manipulation (some overlaps cobra; port the RAVEN-specific logic). |

### 2.8 `analysis/` — RAVEN-specific analyses  *(Phase 5)*
Not in cobrapy core (some exist in cameo/straindesign — evaluate reuse before porting).
| RAVEN | Notes |
|---|---|
| `reporterMetabolites` | Patil & Nielsen reporter-metabolite algorithm. **Port** (cobrapy lacks it). |
| `FSEOF` | Flux Scanning with Enforced Objective Function. **Port** (or wrap straindesign). |
| `runRobustnessAnalysis`, `runSimpleOptKnock` | Evaluate vs `cobra`/`straindesign` before porting. |
| `runDynamicFBA` | cobrapy has no polished dFBA — **port** a lightweight version. |
| `getAllSubGraphs` | Graph connectivity (networkx). |
| `checkRxn`, `getFluxZ`, `followChanged`, `followFluxes` | Small diagnostics — port as needed. |

### 2.9 `comparison/` — multi-model comparison  *(Phase 5)*
| RAVEN | Notes |
|---|---|
| `compareMultipleModels`, `compareRxnsGenesMetsComps` | Structural comparison / Jaccard across models. |
| `compareModels`, `getModelFromHomology`-based diffs | Reporting. |

### 2.10 `plotting/` & `pathway/` — visualization  *(Phase 6, lowest priority)*
| RAVEN | Notes |
|---|---|
| `drawMap`, `drawPathway`, `colorPathway`, `colorSubsystem` | KEGG/pathway maps with flux & expression overlay. |
| `markPathwayWithFluxes`, `markPathwayWithExpression`, `setOmicDataToRxns` | Omics overlay. |
| others (`getColorCodes`, `plotLabels`, …) | Port last; consider Escher/matplotlib instead of MATLAB plotting. |

---

## 3. Explicitly OUT of scope (initially)
- `legacy/` (11 files) — deprecated RAVEN 1 compatibility.
- `software/`, `installation/`, `INIT/startup`, MATLAB-path/Java setup (`addJavaPaths`, `loadWorkbook`) — MATLAB-runtime plumbing.
- `solver/` — replaced wholesale by cobrapy/optlang.
- `external/updateDocumentation`, `getWSLpath`, `parallelPoolRAVEN` — MATLAB infra.
- `testing/` MATLAB tests — replaced by a fresh pytest suite (but reused as oracle data).

---

## 4. Phased roadmap

| Phase | Theme | Deliverables | Depends on |
|---|---|---|---|
| **1** | Foundation | `utils/` helpers (`is_dnf`/`find_non_dnf_grrules` ✅, `parse_name_comp`, `checkModelStruct`, MIRIAM/annotation + ID-prefix — **no** struct adapter; `getIndexes` and grRule *normalization* **not** ported, cobra covers them) and the `manipulation/` ergonomic layer (§1b: `addRxns` & co., `setParam`, `removeReactions` & co., `simplifyModel`, `mergeModels`, `sortModel`), packaging, CI, pytest skeleton. Migration cheatsheet doc. | — |
| **2** | I/O | `readYAMLmodel`/`writeYAMLmodel`, Excel import/export, tab-delimited & SIF export. | 1 |
| **3** | Reconstruction | homology (BLAST/DIAMOND) → KEGG → MetaCyc reconstruction. | 1, 2 |
| **4** | Context-specific & tasks | metabolic `tasks/`, `gapfilling/`, then tINIT/ftINIT. (Tasks first — INIT depends on them.) | 1, 2, MIP solver |
| **5** | Data integration & analysis | HPA/omics, localization, `reporterMetabolites`, FSEOF, dFBA, model comparison. | 1–4 |
| **6** | Visualization | pathway maps / omics overlay (consider Escher). | 1–2 |

**Suggested order rationale:** each phase produces something usable on its own. Reconstruction
(Phase 3) is RAVEN's headline feature and only needs the foundation + I/O. tINIT (Phase 4)
depends on the task framework, so tasks are built first within the same phase.

---

## 5. Cross-cutting concerns
- **Validation strategy:** for every ported function, run the MATLAB original on a small fixture
  (e.g. *S. cerevisiae* / *S. coelicolor* draft) and store its output as a reference oracle in
  `tests/data/` to assert numerical/structural parity.
- **MIP solver:** tasks/INIT/gap-filling need MILP — document Gurobi/CPLEX/SCIP setup; GLPK MILP
  is too slow for genome scale.
- **External tools:** reconstruction wraps BLAST+/DIAMOND/HMMER binaries via `subprocess`; KEGG/
  MetaCyc access via REST/flat-files with on-disk caching (respect KEGG API terms; MetaCyc license).
- **Reference data:** the `keggModelFiles`/MetaCyc dumps RAVEN ships should be downloaded/cached,
  not vendored, to keep the repo small.

---

## 6. Resolved decisions
1. **Package name:** `ravengem` — both PyPI distribution and import name (`pip install ravengem`,
   `import ravengem`). The obvious `ravenpy` is taken on PyPI by an unrelated hydrology package. ✅
2. **Repo home:** `SysBioChalmers/ravengem`. ✅
3. **tINIT/ftINIT & tasks:** functional equivalence is sufficient (not bit-exact). ✅
4. **KEGG/MetaCyc data:** support both live REST (with disk cache) and reused RAVEN dumps,
   selectable via configuration. ✅
5. **No RAVEN⇄COBRA adapter:** `ravenCobraWrapper` is not ported; ravengem adheres to the
   `cobra.Model` format directly. ✅
6. **YAML I/O:** follow the cobrapy YAML standard plus geckopy's ec extension keys
   (`ec-rxns`, `ec-enzymes`, `gecko_light`, `metaData`). ✅

### Still open
- **geckopy → ravengem dependency direction** (see §7).

---

## 7. Code relocated from geckopy
During the geckopy refactor, some RAVEN functionality was ported to Python as generic
`cobra.Model` operations. ravengem is their canonical home.

| function(s) | RAVEN origin | ravengem location | status |
|---|---|---|---|
| `convert_to_irreversible` | `convertToIrrev.m` | `manipulation/irreversible.py` | ✅ adopted + tests passing |
| `expand_model`, `_gpr_to_dnf`, `_node_to_dnf` | `expandModel.m` | `manipulation/expand.py` | ✅ adopted + tests passing |

**Stays in geckopy** (enzyme-constrained or GECKO-only, despite RAVEN mentions in comments):
`ec_fseof` (= `ecFSEOF.m`, operates on `EcModel`/protein usage), `remove_pseudoreaction_gprs`,
`invert_backwards_only_reactions`, and `setParam`/legacy-YAML MATLAB-COMPAT notes.

**Resolved — strategy: adopt copy now, converge later.** ✅ ravengem holds the canonical copy
(provenance comments preserved; geckopy's pytest cases adopted as the initial suite — 34 tests
pass under cobra 0.31). geckopy keeps its own copy untouched for now and will switch to importing
from `ravengem.manipulation` once ravengem is published on PyPI. Until then there is brief, tracked
duplication; keep the two implementations in sync if either changes.
