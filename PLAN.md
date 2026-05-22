# ravenpy — Port Plan

A plan for porting the [RAVEN Toolbox 2](https://github.com/SysBioChalmers/RAVEN)
(MATLAB, v2.11.1, ~271 `.m` files) to a Python package built on **cobrapy**.

The guiding principle is **do not re-port what cobrapy already provides.** RAVEN and cobrapy
overlap heavily on simulation and model-manipulation primitives. The value of `ravenpy` is in
RAVEN's *reconstruction* and *context-specific modeling* features, which cobrapy lacks. This
document (1) maps every RAVEN functional area to its cobrapy equivalent or "port" verdict, and
(2) lays out a phased roadmap.

---

## 0. Foundational design decisions

| Decision | Choice | Rationale |
|---|---|---|
| In-memory model object | `cobra.Model` (not a re-implemented RAVEN struct) | Reuse cobrapy's data model, solver interface, SBML I/O, ecosystem interop. |
| RAVEN-only fields (`rxnMiriams`, `metDeltaG`, `rxnConfidenceScores`, `compMiriams`, `geneShortNames`, `eccodes`, `inchis`, …) | Stored in `cobra` `.annotation` / `.notes`; `subSystems` → SBML groups | cobrapy already round-trips these through SBML; no parallel struct needed. |
| RAVEN↔COBRA conversion | **None.** ravenpy functions consume and produce `cobra.Model` directly — there is no `ravenCobraWrapper` port and no parallel RAVEN struct | Adhere to the cobrapy model format throughout; full COBRA-ecosystem interop for free, nothing to keep in sync. |
| Solver | cobrapy's optlang (GLPK/Gurobi/CPLEX/HiGHS) | Replaces RAVEN's `solveLP`/`solveQP`/`optimizeProb` MILP/LP layer. MILP for INIT/gap-filling needs a MIP-capable solver (Gurobi/CPLEX/SCIP). |
| Naming on PyPI | `raven-toolbox` (import as `ravenpy`) | ⚠️ `ravenpy` on PyPI is taken by an unrelated hydrology package. Distribution name mirrors the "RAVEN Toolbox" branding; import name stays `ravenpy`. |
| Repo home | `SysBioChalmers/ravenpy` | Alongside the MATLAB RAVEN. |
| tINIT/ftINIT & tasks fidelity | **Functional equivalence** (not bit-exact) | Same algorithm/intent; minor numerical differences from solver behavior are acceptable. |
| KEGG/MetaCyc data source | **Configurable**: live REST (+disk cache) *or* reuse RAVEN's pre-built dumps | Live for currency; dumps for reproducibility with MATLAB results. |
| License | GPL-3.0-or-later | Derivative of GPLv3 RAVEN. |

---

## 1. Already in cobrapy — DO NOT PORT (use/wrap instead)

These RAVEN functions have direct cobrapy equivalents. `ravenpy` should provide, at most, thin
convenience wrappers or documentation mapping the old names.

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
| `addRxns`, `addMets`, `addExchangeRxns`, `addTransport` | `Model.add_reactions`, `add_metabolites`, `add_boundary` |
| `removeReactions`, `removeMets`, `removeGenes` | `Model.remove_reactions`, `remove_metabolites`, `cobra.manipulation.remove_genes` |
| `changeRxns`, `setParam`, `setExchangeBounds` | direct attribute assignment (`rxn.bounds`, `model.medium`, `model.objective`) |
| `getExchangeRxns`, `getTransportRxns` | `model.exchanges`, `model.boundary`, custom filters |
| `constructEquations`, `buildEquation`, `parseRxnEqu` | `reaction.build_reaction_string`, `reaction.reaction`, `reaction.build_reaction_from_string` |
| `convertToIrrev` | manual split, or `cobra.util.array` helpers |
| `constructS` | `cobra.util.create_stoichiometric_matrix` |
| `getElementalBalance` | `reaction.check_mass_balance`, `metabolite.elements` |
| `deleteUnusedGenes`, parts of `simplifyModel` | `cobra.manipulation.prune_unused_metabolites`, `prune_unused_reactions` |
| `mergeModels` (basic) | `model.merge` |
| `printFluxes`, `printModel`, `printModelStats` | `model.summary()`, `reaction.summary()`, `metabolite.summary()` |
| `standardizeGrRules`, `changeGrRules`, `getGenesFromGrRules` | `cobra.core.gene.GPR` (parse/eval/AST) |
| `parseFormulas` | `metabolite.formula` / `elements` |
| `getIndexes`, `sortModel`, `sortIdentifiers`, `permuteModel` | native Python indexing / cobra DictList |

**Verdict:** ~70 RAVEN functions collapse into cobrapy calls. Capture these as a "migration
cheatsheet" in the docs rather than as code.

---

## 2. RAVEN-UNIQUE — the real port targets

Organized by the `src/ravenpy/` subpackage that will host them.

### 2.1 `utils/` — model helpers (NO struct adapter)  *(Phase 1, foundational)*
The `ravenCobraWrapper` adapter is **explicitly not ported** — ravenpy works on `cobra.Model`
directly. Only the genuinely useful, RAVEN-flavored helpers survive, as small functions on cobra
objects:
| RAVEN | Notes |
|---|---|
| `checkModelStruct` | Validate a `cobra.Model` against RAVEN reconstruction expectations (beyond cobra's own `model.validate` / SBML validation). |
| `editMiriam`, `extractMiriam` | Convenience get/set for MIRIAM-style entries inside cobra `.annotation` dicts. |
| `addIdentifierPrefix`, `removeIdentifierPrefix` | `R_`/`M_`/`G_` prefix handling for interop with COBRA-Toolbox-style IDs (only where cobra's SBML layer doesn't already cover it). |
| ~~`ravenCobraWrapper`, `standardizeModelFieldOrder`, `cobraNamespaces.csv`, `COBRA_structure_fields.csv`~~ | **Dropped** — no parallel struct to convert or order. |

### 2.2 `io/` — RAVEN-specific formats  *(Phase 1–2)*
| RAVEN | Notes |
|---|---|
| `readYAMLmodel`, `writeYAMLmodel` | **Follow the cobrapy YAML standard** (`cobra.io.dict.model_to_dict` / `model_from_dict` for the metabolites/reactions/genes/compartments portion), and **additionally support geckopy's enzyme-constrained extension keys** (`ec-rxns`, `ec-enzymes`, `gecko_light`, `metaData`) so ecModels round-trip. Empty/NaN fields omitted; loader fills them back. `ruamel.yaml`. High priority — heavily used. |
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
| **1** | Foundation | `utils/` model helpers (`checkModelStruct` validation, MIRIAM/annotation + ID-prefix helpers — **no** struct adapter), packaging, CI, pytest skeleton. Migration cheatsheet doc. | — |
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
1. **PyPI name:** `raven-toolbox` (import as `ravenpy`). ✅
2. **Repo home:** `SysBioChalmers/ravenpy`. ✅
3. **tINIT/ftINIT & tasks:** functional equivalence is sufficient (not bit-exact). ✅
4. **KEGG/MetaCyc data:** support both live REST (with disk cache) and reused RAVEN dumps,
   selectable via configuration. ✅

5. **No RAVEN⇄COBRA adapter:** `ravenCobraWrapper` is not ported; ravenpy adheres to the
   `cobra.Model` format directly. ✅
6. **YAML I/O:** follow the cobrapy YAML standard plus geckopy's ec extension keys
   (`ec-rxns`, `ec-enzymes`, `gecko_light`, `metaData`). ✅

### Still open
- Final package / PyPI name (see brainstorm below — to be decided).
