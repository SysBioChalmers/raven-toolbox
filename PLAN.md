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
| `buildEquation`, `parseRxnEqu`, `constructEquations` | `reaction.reaction`, `reaction.build_reaction_string(use_metabolite_names=...)`, `reaction.build_reaction_from_string` |
| `addMets` | `model.add_metabolites([cobra.Metabolite(...), ...])` |
| `addExchangeRxns` | `model.add_boundary(met, type="exchange" / "demand" / "sink")` |
| `setParam` (`lb`/`ub`/`eq`/`obj`/`unc`) | `rxn.bounds = (lo, hi)` / `rxn.bounds = (v, v)` / `model.objective = {rxn: c}` / `rxn.bounds = cobra.Configuration().bounds`; loop for batch. (`var` band → `set_variance_bounds`.) |
| `setExchangeBounds`, `getExchangeRxns` | `model.medium = {ex_id: uptake}` (sets uptake, closes others, handles direction); `model.exchanges` / `model.sinks` / `model.demands` |
| `simplifyModel` `deleteMinMax` | `model.remove_reactions(cobra.flux_analysis.find_blocked_reactions(model))` |
| `simplifyModel` `deleteZeroInterval` | `model.remove_reactions([r for r in model.reactions if r.bounds == (0, 0)])` then `cobra.manipulation.prune_unused_metabolites` |
| `editMiriam`, `extractMiriam` | `met.annotation["kegg.compound"] = "C00031"` / `met.annotation.get(...)` (annotation is a `{namespace: id(s)}` dict) |
| `addIdentifierPrefix`, `removeIdentifierPrefix` | handled by `cobra.io` SBML read/write |
| `constructS` | `cobra.util.create_stoichiometric_matrix` |
| `printFluxes`, `printModel`, `printModelStats` | `model.summary()`, `reaction.summary()`, `metabolite.summary()` |
| `getGenesFromGrRules` | `cobra.core.gene.GPR` (parse/eval/AST) |
| `getMetsInComp` | `[m for m in model.metabolites if m.compartment == c]` |
| `getRxnsInComp` | `[r for r in model.reactions if r.compartments == {c}]` (fully in `c`); use `c in r.compartments` to include transport |
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
| `addRxns` | **PORT (keystone)** ✅ | Done as `add_reactions_from_equations` ([manipulation/add.py](src/ravengem/manipulation/add.py)). cobra's `build_reaction_from_string` already parses equations/arrows/coeffs and auto-creates ID-matched mets, and setting `gene_reaction_rule` auto-creates genes — so the port keeps only what cobra lacks: **name-based matching** (`mets_by="name"`) and **`name[comp]` syntax**, **compartment** assignment for new mets (cobra leaves them `None`), **strict** policies (error vs auto-create mets/genes), and a duplicate-ID guard. RAVEN's `eqnType` 1/2/3 integer became the readable `mets_by` keyword (+ auto `name[comp]`). |
| `addRxnsGenesMets` | **PORT** ✅ | Done as `add_reactions_from_model` ([manipulation/transfer.py](src/ravengem/manipulation/transfer.py)). Copies reactions from a source model into a draft, matching mets by `name[comp]` (not id), skipping duplicates, adding only genuinely new mets (copying id/formula/charge/annotation) and genes. cobra's `merge` is strict-by-id. |
| `addTransport` | **PORT** ✅ | Done as `add_transport_reactions` ([manipulation/transport.py](src/ravengem/manipulation/transport.py)). Creates transport reactions from one compartment to many, matching mets **by name** across comps, sequential `tr_0001` IDs, optionally creating the target metabolite (copying formula/charge/annotation). cobra has **no** transport primitive at all. |
| `changeRxns` | **PORT** ✅ | Done as `change_reaction_equations` ([manipulation/change.py](src/ravengem/manipulation/change.py)). Replaces stoichiometry from equation strings, reusing the `add` parser (id/name/`name[comp]`). cobra edits the same `Reaction` object in place, so RAVEN's remove→re-add→re-sort dance is unnecessary — other fields and order are preserved automatically. Bounds left unchanged, per RAVEN. |
| `changeGrRules` | **PORT** ✅ | Done as `change_gene_reaction_rules` ([manipulation/change.py](src/ravengem/manipulation/change.py)). Batch-set with an **append** mode (`(old) or (new)`); gene auto-creation and normalization come free from cobra's `gene_reaction_rule=` setter. |
| `setParam` | **MOSTLY DO NOT PORT** | All modes but one are cobra one-liners (`lb`/`ub`/`eq` → `reaction.bounds`; `obj` → `model.objective`; `unc` → `Configuration().bounds`; batch → a loop) — see §1 cheatsheet. The only mode without a clean cobra idiom, the `var` ±% band, is kept as the focused `set_variance_bounds` ([manipulation/parameters.py](src/ravengem/manipulation/parameters.py)) ✅. A 6-mode `set_parameters` was built then trimmed (review: not Pythonic, mostly re-wrapped cobra). |
| `setExchangeBounds` | **DO NOT PORT** | cobra's `model.medium = {ex_id: uptake}` already does the substance — sets uptake bounds, **closes other exchanges' uptake** (`closeOthers`), handles import direction **per reaction** (RAVEN needs uniform direction), and `model.exchanges` excludes sinks/demands (`mediaOnly`). Gaps are thin: keying by metabolite *name* (vs exchange-rxn id) and setting the secretion bound (a one-line `rxn.upper_bound=`). See §1 cheatsheet. |
| `removeReactions` | **DO NOT PORT** | Decided not to separate orphan-metabolite from orphan-gene cleanup; coupled, it is exactly `cobra.Model.remove_reactions(remove_orphans=...)`. Use cobra directly. |
| `removeMets` | **PORT (thin)** ⚠️ | Done ([manipulation/remove.py](src/ravengem/manipulation/remove.py)). Delegates to cobra; the **only** add is `by_name` (delete a metabolite across all compartments). That need is likely rare — flagged as a **deletion candidate** if unused. |
| `removeGenes` | **PORT** ✅ | Done. cobra's `remove_genes` already rewrites GPRs via AST with correct AND/OR semantics (RAVEN did this with `eval`); the port adds the `blocked_reactions` policy — `"remove"` / `"constrain"` (bounds→0, RAVEN default) / `"keep"` — for reactions left with no enzyme. |
| `simplifyModel` | **PORT (gap modes)** ✅ | Done as focused functions in [manipulation/simplify.py](src/ravengem/manipulation/simplify.py): `remove_dead_end_reactions` (deleteInaccessible), `remove_duplicate_reactions` (deleteDuplicates), `constrain_reversible_reactions` (FVA), `group_linear_reactions` (lossy). The cobra-covered modes are **not** ported: `deleteMinMax`→`find_blocked_reactions`, `deleteZeroInterval`→filter+prune, `deleteUnconstrained`→moot (§1 cheatsheet). |
| `mergeModels` | **PORT** ✅ | Done as `merge_models` ([manipulation/merge.py](src/ravengem/manipulation/merge.py)). Merges **N** models, unifying mets by `name[comp]` (or id), keeping **all** reactions (RAVEN does not dedup; id collisions renamed `id_<source>`), merging genes, tracking provenance in `notes['origin']`. cobra's `merge` is pairwise and strict-by-id. |
| `sortModel` (core) | **PORT** core / **SKIP** optimizer | Port the **deterministic** canonical ordering (mets by `name[comp]`; reversible reactions flipped to lexicographic-first reactant) — exactly what makes YAML diffable. **Skip** the stochastic `sortReactionOrder` annealer. |
| `addMets` | **DO NOT PORT** | `model.add_metabolites([Metabolite(...), ...])` covers batch add; `copyInfo` is niche. `addRxns` already auto-creates new mets. |
| `addExchangeRxns` | **DO NOT PORT** | `model.add_boundary(met, type="exchange"/"demand"/"sink")` covers it; only RAVEN's `EXC_OUT_*` naming differs. |
| `constructEquations` | **DO NOT PORT** | `reaction.build_reaction_string(use_metabolite_names=...)` gives id or name equations; formula-rendering is niche. See §1 cheatsheet. |

### Targets for `utils/` — lookup, GPR hygiene, balance

| RAVEN | verdict | Why |
|---|---|---|
| `getIndexes` | **DO NOT PORT** (keep `name[comp]` sliver only) | RAVEN needs a central index resolver because it's a struct of parallel arrays. cobra is object-oriented and already covers mixed lookup more idiomatically: `DictList.get_by_any` (mixed id/object/**index** → objects), `get_by_id` (O(1)), `query` (name/substring/regex), `index` (position), comprehensions for filtering. A 1-based-index port would be redundant and un-Pythonic. **Keep only** the `name[comp]` composite resolver (`parse_name_comp`) as a small helper for `addRxns`/`addTransport`/`mergeModels` — that's the one bit cobra lacks. |
| `standardizeGrRules` | **PORT lint only** ✅ | Two halves: (1) syntax normalization — **not ported**, cobra auto-normalizes every GPR on assignment (case, whitespace, redundant brackets); (2) the `findPotentialErrors` lint for non-DNF rules — **ported** as `find_non_dnf_grrules` + `is_dnf` ([utils/gpr.py](src/ravengem/utils/gpr.py)), reworked onto cobra's GPR AST and returning structured `GPRIssue`s instead of printing. |
| `getElementalBalance` | **PORT** ✅ | Done as `get_elemental_balance` ([utils/balance.py](src/ravengem/utils/balance.py)). Graded `balanced`/`unbalanced`/`unknown` status — distinguishing a missing formula, which cobra's `check_mass_balance` silently miscounts. |
| `getRxnsInComp`, `getMetsInComp` | **DO NOT PORT** | One-liners over cobra's `reaction.compartments` / `metabolite.compartment`; not worth a wrapper (see §1 cheatsheet). |

---

## 2. RAVEN-UNIQUE — the real port targets

Organized by the `src/ravengem/` subpackage that will host them.

### 2.1 `utils/` — model helpers (NO struct adapter)  *(Phase 1, foundational)*
The `ravenCobraWrapper` adapter is **explicitly not ported** — ravengem works on `cobra.Model`
directly. Only the genuinely useful, RAVEN-flavored helpers survive, as small functions on cobra
objects:
| RAVEN | Notes |
|---|---|
| `checkModelStruct` | **PORT (curation subset)** ✅ — `check_model` ([utils/validate.py](src/ravengem/utils/validate.py)). RAVEN's struct/field-type/duplicate-ID/`lb>ub`/`rev` checks are moot in cobra (object model enforces them); the surviving curation bundle (orphan mets/genes, empty reactions, duplicate name+compartment, empty names, objective sanity) is returned as structured `ModelIssue`s. |
| ~~`editMiriam`, `extractMiriam`~~ | **DO NOT PORT** — cobra's `.annotation` is already a `{namespace: id(s)}` dict; read/write it directly (cheatsheet). |
| ~~`addIdentifierPrefix`, `removeIdentifierPrefix`~~ | **DO NOT PORT** — cobra's SBML layer handles `R_`/`M_`/`G_` prefixing on read/write. |
| `parse_name_comp` (from `getIndexes` `metcomps`) | **PORT** (small helper) ✅ — resolve the `name[comp]` composite (metabolite by name + compartment) cobra can't. Done in [utils/parse.py](src/ravengem/utils/parse.py). The rest of `getIndexes` is **not ported** — `DictList.get_by_any`/`get_by_id`/`query`/`index` already cover mixed id/object/index/name lookup more idiomatically (see §1b note). |
| `standardizeGrRules` | **PORT lint only** ✅ — normalization is cobra-automatic; non-DNF lint ported as `find_non_dnf_grrules`/`is_dnf`. See §1b. |
| `getElementalBalance` | **PORT** ✅ — graded balanced/unbalanced/unknown (catches cobra's silent missing-formula); see §1b. Lives in `utils/balance.py`. |
| `getRxnsInComp`, `getMetsInComp` | **DO NOT PORT** | Too thin over cobra (one-liners via `reaction.compartments` / `metabolite.compartment`); see §1 cheatsheet. Reconsider only if a downstream consumer (e.g. localization) genuinely needs the `include_partial` distinction in several places. |
| ~~`ravenCobraWrapper`, `standardizeModelFieldOrder`, `cobraNamespaces.csv`, `COBRA_structure_fields.csv`~~ | **Dropped** — no parallel struct to convert or order. |

### 2.1b `manipulation/` — model construction, editing & structural transforms  *(Phase 1)*
The home for the RAVEN ergonomic layer (§1b) that *mutates* a `cobra.Model`, plus structural
transforms cobra lacks. Two transforms are **already ported in geckopy** and relocated here (see §7).

**Structural transforms:**
| RAVEN | Notes |
|---|---|
| `convertToIrrev` | Split reversible non-exchange reactions into forward + `_REV` pair. **Already ported** as `convert_to_irreversible` (geckopy `pipeline/preprocess.py`); cobra's old `convert_to_irreversible` was removed, so this is a real port, not a wrapper. |
| `expandModel` | Split isozyme (OR-GPR) reactions into one reaction per AND-clause (`_EXP_N`). **Already ported** as `expand_model` + `_gpr_to_dnf`/`_node_to_dnf` (geckopy `pipeline/expand.py`), using cobra's GPR AST instead of RAVEN string manipulation. |
| `simplifyModel` | **PORT (gap modes)** ✅ — dead-end / duplicate / constrain-reversible / group-linear ported ([manipulation/simplify.py](src/ravengem/manipulation/simplify.py)); cobra-covered modes cheatsheeted. See §1b. |
| `mergeModels` | **PORT** ✅ — N-way merge with name+comp matching, conflict rename, provenance ([manipulation/merge.py](src/ravengem/manipulation/merge.py)); see §1b. |
| `sortModel` | **PORT** deterministic canonical ordering (skip stochastic optimizer); see §1b. |
| `contractModel`, `mergeCompartments`, `copyToComps` | Other RAVEN structural transforms to port as needed. |

**Construction & editing (ergonomic layer, §1b):** `addRxns` (keystone — equation-string batch add
with met/gene auto-creation) ✅, `changeRxns` ✅, `changeGrRules` ✅, `setParam` (→ `set_variance_bounds`
only) ✅, `removeMets`/`removeGenes` ✅, `addTransport` ✅, `addRxnsGenesMets` ✅ — the editing layer
is essentially complete. `setExchangeBounds` is **not ported** (cobra's `model.medium` covers it; §1 cheatsheet).
`addMets`, `addExchangeRxns`, `constructEquations` are **not ported** — cobra covers them (§1
cheatsheet). See §1b for per-function rationale and verdicts.

### 2.2 `io/` — RAVEN-specific formats  *(Phase 1–2)*
| RAVEN | Notes |
|---|---|
| `readYAMLmodel`, `writeYAMLmodel` | **DONE** ✅ as `read_yaml_model`/`write_yaml_model` ([io/yaml.py](src/ravengem/io/yaml.py)). **Aligned to RAVEN commit `fa281a1`** (`feat/geckopy-compat-yaml`, in `/mnt/c/Work/GitHub/RAVEN`), whose writer emits **cobra's native `!!omap`**. Because the format *is* cobra's, cobra reads/writes the standard content **and the whole `annotation` block** (which holds `smiles`, `ec-code`, MIRIAM). This port adds only what cobra drops: (1) RAVEN-only **top-level per-entry keys** — `inchis`/`deltaG`/`metFrom`/`notes`(metNotes) on mets, `confidence_score`/`references`/`rxnFrom`/`deltaG`/`notes` on rxns, `protein` on genes — stashed in `.notes` on read, lifted to top-level on write; (2) model-level `version`, `metaData` provenance, and GECKO `ec-*`/`gecko_light` sections; (3) legacy files with id/name in `metaData`. Note vs the earlier port: `smiles` is **annotation** (cobra-owned), `inchis` is **top-level**, and the top-level `notes` **string** is handled (no longer crashes). |
| `importExcelModel`, `SBMLFromExcel` | **DO NOT PORT** (per decision) — Excel *import* is not enabled in ravengem. |
| `exportToExcelFormat` | Optional Excel *export* (`openpyxl`); export only, lower priority. |
| `exportToTabDelimited`, `exportForGit` | Plain-text / git-friendly model dumps. |
| `exportModelToSIF` | **DONE** ✅ as `export_model_to_sif` ([io/sif.py](src/ravengem/io/sif.py)) — Cytoscape SIF (`rc`/`rr`/`cc` graphs). cobra has no network export. |
| `getToolboxVersion`, `getMD5Hash` | Provenance helpers. |

### 2.3 `reconstruction/` — de novo reconstruction  *(Phase 3, flagship)*
RAVEN's single biggest reason to exist; cobrapy does not cover it at all. There are
**three independent reconstruction approaches**, each a self-contained track with its
own top-level entry point, inputs, and external dependencies. They are planned and
built separately (suggested order **homology → KEGG → MetaCyc**); each can ship on its
own. The shared piece they all rely on is the model-construction layer (§2.1b
`add_reactions_from_equations` etc.) and I/O (§2.2).

#### 2.3a Homology-based — `reconstruction/homology/`  *(Phase 3a)*
Transfer reactions from an existing **template GEM** of a related organism, via
bidirectional best BLAST/DIAMOND hits between protein sequences. The most
self-contained approach (no online database; inputs are a model + FASTA), so built
first.
- **Entry point:** `getModelFromHomology` — given a template `cobra.Model` and a
  homology mapping, transfer reactions and remap genes to the target organism.
- **Inputs:** a template `cobra.Model` + query/target protein FASTA.
- **External deps:** BLAST+ / DIAMOND executables (via `subprocess`).

| RAVEN | Notes |
|---|---|
| `getModelFromHomology` | Transfer reactions from a template GEM via bidirectional best hits. |
| `getBlast`, `getDiamond`, `getBlastFromExcel` | Wrap external BLAST+/DIAMOND executables (subprocess). |
| `makeFakeBlastStructure`, `parseScores` | Homology-score plumbing. |

#### 2.3b KEGG-based — `reconstruction/kegg/`  *(Phase 3b)*
Build a draft GEM from **KEGG** orthology assignments for the organism. Heavier external
data than homology; independent of MetaCyc.
- **Entry point:** `getKEGGModelForOrganism` — orchestrates KO assignment → reaction/
  metabolite/gene retrieval → draft `cobra.Model`.
- **Inputs:** organism KEGG id, or a proteome FASTA (for de novo KO assignment).
- **External deps:** KEGG REST (with on-disk cache) **or** RAVEN's pre-built KEGG dumps
  (configurable, per §0); HMMER for KO assignment.

| RAVEN | Notes |
|---|---|
| `getKEGGModelForOrganism` | Top-level: build a draft GEM for an organism from KEGG. Orchestrates the below. |
| `getModelFromKEGG`, `getRxnsFromKEGG`, `getMetsFromKEGG`, `getGenesFromKEGG` | Parse KEGG flat-files / REST into a model. |
| `getPhylDist` | Phylogenetic-distance weighting of KEGG orthologs. |
| `constructMultiFasta` | Build per-KO FASTA for homology search. |

#### 2.3c MetaCyc-based — `reconstruction/metacyc/`  *(Phase 3c)*
Build a draft from **MetaCyc** reactions/pathways, and optionally reconcile it with a
KEGG draft (so this track can build on 3b but does not require it).
- **Entry point:** `getMetaCycModelForOrganism` — draft `cobra.Model` from MetaCyc;
  `combineMetaCycKEGGModels` merges a MetaCyc and a KEGG draft.
- **Inputs:** MetaCyc flat-file dumps (license-restricted); optionally a KEGG draft from 3b.
- **External deps:** MetaCyc data dumps; cross-DB linking to KEGG (`linkMetaCycKEGGRxns`).

| RAVEN | Notes |
|---|---|
| `getMetaCycModelForOrganism`, `getModelFromMetaCyc` | MetaCyc-based draft reconstruction. |
| `getRxnsFromMetaCyc`, `getMetsFromMetaCyc`, `getEnzymesFromMetaCyc` | MetaCyc flat-file parsers. |
| `linkMetaCycKEGGRxns`, `combineMetaCycKEGGModels`, `addSpontaneousRxns` | Cross-DB reconciliation (with the KEGG track, 3b). |

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
| **1** | Foundation | `utils/` helpers (`is_dnf`/`find_non_dnf_grrules` ✅, `get_elemental_balance` ✅, `check_model` ✅, `parse_name_comp` ✅ — **no** struct adapter; `getIndexes`, grRule *normalization*, MIRIAM/ID-prefix helpers **not** ported, cobra covers them) and the `manipulation/` ergonomic layer (§1b: largely done — add/change/remove/transport/transfer/variance ✅), packaging, CI, pytest skeleton. Migration cheatsheet doc. | — |
| **2** | I/O | `readYAMLmodel`/`writeYAMLmodel`, Excel import/export, tab-delimited & SIF export. | 1 |
| **3a** | Reconstruction — homology | `getModelFromHomology` + BLAST/DIAMOND wrappers (§2.3a). Self-contained; build first. | 1, 2 |
| **3b** | Reconstruction — KEGG | `getKEGGModelForOrganism` + KEGG retrieval/KO assignment (§2.3b). | 1, 2 |
| **3c** | Reconstruction — MetaCyc | `getMetaCycModelForOrganism` + MetaCyc parsers + KEGG reconciliation (§2.3c). | 1, 2, (3b for combine) |
| **4** | Context-specific & tasks | metabolic `tasks/`, `gapfilling/`, then tINIT/ftINIT. (Tasks first — INIT depends on them.) | 1, 2, MIP solver |
| **5** | Data integration & analysis | HPA/omics, localization, `reporterMetabolites`, FSEOF, dFBA, model comparison. | 1–4 |
| **6** | Visualization | pathway maps / omics overlay (consider Escher). | 1–2 |

**Suggested order rationale:** each phase produces something usable on its own. Reconstruction
(Phase 3) is RAVEN's headline feature and only needs the foundation + I/O; it splits into three
**independent** tracks (3a homology, 3b KEGG, 3c MetaCyc), built in that order — homology is most
self-contained, and MetaCyc can optionally reconcile against a KEGG draft. tINIT (Phase 4)
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
