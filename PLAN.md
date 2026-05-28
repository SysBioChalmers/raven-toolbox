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
| `analyzeSampling` (post-hoc stats on a sample set) | pandas on the sample matrix (mean/std/quantiles) |
| `runProductionEnvelope`, `runPhenotypePhasePlane` | `cobra.flux_analysis.production_envelope` |
| `getExchangeRxns`, `getTransportRxns` | `model.exchanges`, `model.boundary`, custom filters |
| `buildEquation`, `parseRxnEqu`, `constructEquations` | `reaction.reaction`, `reaction.build_reaction_string(use_metabolite_names=...)`, `reaction.build_reaction_from_string` |
| `addMets` | `model.add_metabolites([cobra.Metabolite(...), ...])` |
| `addExchangeRxns` | `model.add_boundary(met, type="exchange" / "demand" / "sink")` |
| `setParam` (`lb`/`ub`/`eq`/`obj`/`unc`) | `rxn.bounds = (lo, hi)` / `rxn.bounds = (v, v)` / `model.objective = {rxn: c}` / `rxn.bounds = cobra.Configuration().bounds`; loop for batch. (`var` band → `set_variance_bounds`.) |
| `setExchangeBounds`, `getExchangeRxns` | `model.medium = {ex_id: uptake}` (sets uptake, closes others, handles direction); `model.exchanges` / `model.sinks` / `model.demands` |
| `fillGaps` (`useModelConstraints=true`, fill toward the objective) | `cobra.flux_analysis.gapfill(model, universal, lower_bound=…, penalties=…)` — **align the template's metabolite ids to the draft first** (cobra matches by id; RAVEN matched by name). The *other* `fillGaps` mode (connectivity — unblock blocked reactions) has no cobra equivalent → `gapfilling.connect_blocked_reactions` (§4b). |
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
| `permuteModel` | native Python indexing / cobra DictList (`sortIdentifiers` → planned `sort_identifiers`, §2.2) |

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
| `exportToExcelFormat` | **DONE** ✅ as `export_to_excel` ([io/excel.py](src/ravengem/io/excel.py)) — export only (RAVEN 5-sheet xlsx: RXNS/METS/COMPS/GENES/MODEL), `openpyxl` (lazy, `[excel]` extra). cobra has no Excel I/O. |
| `exportToTabDelimited` | **DO NOT PORT** — served the excluded Excel/tab importer; a generic table is a pandas one-liner (cheatsheet). |
| `exportForGit` | **DONE** ✅ as `export_for_git` ([io/git.py](src/ravengem/io/git.py)) — standard-GEM repo layout (`model/<fmt>/<prefix>.<fmt>`) for `yml`/`xml`/`mat`/`xlsx`/`txt` via `write_yaml_model` / cobra `write_sbml_model` / `save_matlab_model` / `export_to_excel` / single-file reaction table, plus `dependencies.txt` (python/cobra/ravengem). Sorts a copy via `sort_identifiers` first (caller's model untouched). |
| `sortIdentifiers` | **DONE** ✅ as `sort_identifiers` ([utils/sort.py](src/ravengem/utils/sort.py)) — model-wide alphabetical `DictList.sort`; plus `sort_ids=` on `write_yaml_model` (sorts the output, not the model). |
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
| `getModelFromHomology` | → `get_model_from_homology` — the core algorithm; pure Python on cobra models + a hits table. |
| `getBlast`, `getDiamond` | → `run_blast` / `run_diamond` — subprocess wrappers (BLAST+ / DIAMOND). |
| `getBlastFromExcel` | → `blast_from_table` — load hits from a CSV/DataFrame (no Excel). |
| `makeFakeBlastStructure`, `parseScores` | → `make_ortholog_hits` + the tabular parser inside `run_blast`. |

**Design (Phase 3a):**

*Data structure — `blastStructure` → a tidy DataFrame.* RAVEN's struct array of directional
hit sets becomes one `pandas.DataFrame` of bidirectional hits with columns
`from_id, to_id, from_gene, to_gene, evalue, identity, align_len, bitscore, ppos`
(one row per hit). This makes the strictness/best-hit filtering plain groupby logic. A thin
`BlastResult` wrapper (or just the DataFrame with a documented schema) is the currency between
the BLAST wrappers and `get_model_from_homology`.

*Functions (`reconstruction/homology/`):*
- **`make_ortholog_hits(ortholog_pairs, source_model_id, target_id)`** (makeFakeBlastStructure) —
  build a hits DataFrame from a predefined ortholog list, sentinel metrics (e=0, ident=100,
  len=1000) so all pass filters. Pure Python — the testing/seeding entry point.
- **`get_model_from_homology(...)`** (getModelFromHomology) — the core: filter hits, build the
  ortholog map, rewrite GPRs, transfer reactions, merge templates by name[comp]. **Detailed design +
  proposed logic improvements in [docs/plan_get_model_from_homology.md](docs/plan_get_model_from_homology.md)**
  (clearer `bidirectional`/`best_hits_only` params replacing the overloaded `strictness`; **AST-based**
  GPR rewriting instead of regex; explicit `complex_policy` for unmapped AND-subunits; bitscore-based
  best hits; DataFrame ortholog map; provenance). Pure Python on `cobra.Model`; reuses §1b
  reaction-transfer/merge. **Fully testable via `make_ortholog_hits`** (no BLAST).
- **`run_blast(organism_id, fasta, model_ids, ref_fastas, *, evalue=1e-5)`** (getBlast) — bidirectional
  BLAST+ via `subprocess` (`makeblastdb` + `blastp -outfmt 6 'qseqid sseqid evalue length pident
  bitscore ppos'`), parse tabular output → hits DataFrame. Detect binaries via `shutil.which`,
  clear error if missing.
- **`run_diamond(...)`** (getDiamond) — same contract via DIAMOND (`diamond makedb`/`blastp`).
- **`ensure_binary("diamond"/"blastp", version=...)`** (shared `binaries.py`) — registry-driven
  download of the version-pinned ZIP from ravengem releases into a cache; the auto-resolve fallback
  used by all wrappers (generic across tools, not DIAMOND-specific).
- **`blast_from_table(path_or_df, ...)`** (getBlastFromExcel) — accept a precomputed homology
  table (CSV/DataFrame), validate columns → hits DataFrame.

*External tools & binary access (decided): a generic, version-pinned binary provisioner.*
BLAST+/DIAMOND are compiled binaries — **no official pip wheels** (the PyPI names `diamond`/`blast`
are *unrelated* packages; never depend on them). ravengem ships its own **version-pinned binaries as
ZIPs attached to ravengem GitHub releases**, served through one shared module (not homology-specific —
reused by KEGG 3b's HMMER, etc.). Lives in a top-level **`ravengem/binaries.py`** (or
`ravengem/external/`):
- **`ensure_binary(tool, version=None) -> Path`** — consult a **registry** mapping
  `(tool, version, os, arch) → {url, sha256}` (urls = ravengem release assets); download the pinned
  ZIP into a `platformdirs` cache, **verify SHA256**, `chmod +x`, return the path. Actionable error
  (conda/manual options) if the `(os, arch)` combo isn't hosted.
- **Resolution order in every wrapper:** explicit `binary=` arg → env var
  (`RAVENGEM_DIAMOND` / `RAVENGEM_BLASTP` / ...) → `shutil.which` (system/conda/apt/brew) →
  `ensure_binary` (fetch pinned ZIP) → error. So a pre-installed binary always wins; the bundle is
  the zero-setup fallback.
- Pins exact tool versions → **reproducible reconstruction**; uniform across diamond/blast+/future
  tools. bioconda (`conda install -c bioconda diamond blast`) documented as *one* option, not primary.

Maintainer workflow for updating versions and building minimal ZIPs is documented in
[docs/maintaining_binaries.md](docs/maintaining_binaries.md) (registry format, asset naming, the
exact executables to ship — BLAST+ = only `blastp`+`makeblastdb` — stripping/compression, licensing).

*Follow-ups this implies (tracked, not blocking 3a core):* a CI build/release pipeline to produce the
per-OS/arch ZIPs (coverage = what we build: Linux x86-64 first; macOS/ARM/Windows as built); license
compliance for redistribution (BLAST+ = US-gov public domain ✓; DIAMOND = GPLv3, ship licence/notice);
published SHA256s. `platformdirs` becomes a (small) dependency for the cache dir.

*Test strategy:* the **core (`get_model_from_homology` + `make_ortholog_hits` + tabular parsing) is
tested without any external tool**; `run_blast`/`run_diamond` execution gets a `skipif`-guarded test
(only when the binary is present), with command-construction + output-parsing unit-tested against a
captured tabular fixture.

*Build order:* `make_ortholog_hits` + hits schema → `get_model_from_homology` → `run_blast`/
`run_diamond` → `blast_from_table`.

*Improvements to log:* struct-array→DataFrame (filterable); `getBlastFromExcel`→CSV (no Excel,
consistent with the Excel-import exclusion); core testable without BLAST installed.

#### 2.3b KEGG-based — `reconstruction/kegg/`  *(Phase 3b)*
Build a draft GEM from **KEGG** orthology (KO) assignments. cobra covers none of this. The track
is a **pipeline of five sub-steps** — the first three build a *shared, reusable* KEGG reference +
HMM library (done once per KEGG version), the last two build a model for a *specific organism* in
two modes. RAVEN's `getKEGGModelForOrganism` documents exactly these stages via its `dataDir`
sub-folders (`keggdb` / `fasta` / `aligned` / `hmms`).

| Step | ravengem (proposed) | RAVEN | What it does |
|---|---|---|---|
| **3b.1 Obtain KEGG** (maintainer, paid FTP) ✅ | `download_kegg_dump(dest)` / `fetch_kegg_files` + `extract_kegg_dump` (`reconstruction/kegg/download.py`) | bulk FTP dump (ports `fetch_keggdb.sh`) | **Implemented.** A maintainer with a **paid KEGG FTP subscription** pulls the bulk dump (reaction/compound/glycan/ko archives + euk/prok proteomes + taxonomy) **once per KEGG release**; arranged into the flat layout 3b.2/3b.3 consume. **Pure stdlib** (`urllib`/`tarfile`/`gzip`/`netrc`) — no wget/tar/gunzip/Cygwin; credentials via `~/.netrc` (see [docs/maintaining_kegg_data.md](docs/maintaining_kegg_data.md)). End users never need KEGG access — they fetch the published ravengem artefacts (see Data access). |
| **3b.2 Parse dump → reference model + tables** ✅ | `parse_kegg_dump(keggdb_dir, out)` (`reconstruction/kegg/parse.py`) | `getModelFromKEGG` + `getRxnsFromKEGG`/`getMetsFromKEGG`/`getGenesFromKEGG` | **Implemented.** Parse into a **gene-free** reference GEM (reactions + metabolites only) **and** minimal **gzipped-TSV** tables (`ko_reaction`, `ko_names`, `organism_gene_ko`, `rxn_flags`). `phyl_dist` deferred to `getPhylDist` (3b.5). Maintainer-side; published as artefacts. |
| **3b.3 Construct HMMs** (maintainer, build our own) ✅ | `build_ko_fastas` → `build_ko_hmm` → `build_hmm_library(domain=…)` (`reconstruction/kegg/hmm.py` + `taxonomy.py`) | `constructMultiFasta` + align + `hmmbuild`/`hmmpress` | **Implemented.** We build the HMMs ourselves (no KOfam, no third-party pre-built), from **all KEGG organisms** using the paid-FTP bulk sequences from 3b.1 — split by domain into **prok90 / euk90** libraries (taxonomy-driven). Per KO: gather member-gene sequences → multi-FASTA → CD-HIT (~90 %) → MAFFT align → `hmmbuild`; then concatenate + `hmmpress` into one pressed library (K7). Binaries via the `binaries.py` registry (HMMER/MAFFT/CD-HIT). Maintainer-side, once per KEGG release; output published as a ravengem artefact. |
| **3b.4 Model for a KEGG species** ✅ | `get_kegg_model_for_organism(organism_id, …)` / `…_from_artefacts` (`reconstruction/kegg/organism.py`) | `getKEGGModelForOrganism(organismID)` (no-FASTA branch) | **Implemented.** For an organism already in KEGG: read its gene↔KO from `organism_gene_ko`, map KO→reaction via `ko_reaction`, OR-join the organism's genes into each reaction's GPR (genes added here), keep reactions with genes (+ spontaneous when `keep_spontaneous`), and apply the `keep*` quality filters from `rxn_flags`. No homology search. Domain mode (`eukaryotes`/`prokaryotes`) deferred to 3b.5 (needs `getPhylDist`). |
| **3b.5 Model by HMM sequence query** ✅ | `get_kegg_model_from_sequences(fasta, reference_model, ko_reaction, library, …)` (`reconstruction/kegg/query.py`) | `getKEGGModelForOrganism(fastaFile)` (FASTA branch) | **Implemented.** One `hmmscan` of a proteome FASTA against the pressed library (3b.3) → `assign_kos` applies the E-value cut-off + the two score-ratio filters (`min_score_ratio_ko`/`_g`) → shared assembler builds the draft. The de-novo path for organisms not in KEGG. **Phyl-dist subsampling dropped** (our fixed prok90/euk90 libraries make per-organism distance weighting moot; domain choice = library choice). `getPhylDist` distance-matrix therefore not ported; domain mode (3b.4) uses the taxonomy classification directly. |

- **Inputs:** KEGG dump (3b.1); then either a KEGG `organism_id` (3b.4) or a proteome FASTA (3b.5).
- **Data access (decided):** a **maintainer** with a **paid KEGG FTP** account builds, **once per
  KEGG release**, ravengem's own version-pinned artefacts: (i) the parsed KEGG **reference model** +
  KO map (3b.2) and (ii) the **prok90/euk90 HMM library** built from all KEGG (3b.3). These are
  published via the data/release registry (same mechanism as binaries). **End users just download
  the pinned artefacts** — no KEGG account, no FTP, no REST, no per-user harvest. The REST API is
  *not* used for bulk building (dropped — too slow/rate-limited for all of KEGG).
- **Licensing (resolved ✅):** a KEGG **redistribution licence has been obtained**, so ravengem may
  publish the KEGG-derived artefacts (reference model, SQLite tables, HMMs) as version-pinned
  downloads. (Maintainer still uses the paid FTP subscription to build them once per release.)
- **External tools:** **HMMER** (`hmmbuild`/`hmmpress`/`hmmsearch`), an aligner (**MAFFT**), and
  optionally **CD-HIT** (identity dereplication) — all via the shared `binaries.py` `ensure_binary`
  registry (add `hmmer`/`mafft`/`cd-hit` bundles; same pattern as BLAST/DIAMOND in 3a).
  `getPhylDist` → `phylogenetic_distance` helper for 3b.5 score weighting.
- **Storage & distribution (decided):** **not** RAVEN's `.mat` structs.
  - **Reference GEM = gene-free** — only reactions + metabolites (the chemistry), **no genes/GPRs**
    (organism genes number in the millions → would dwarf the model). Stored as **gzipped
    RAVEN/cobra YAML** (`reference_model.yml.gz`) — RAVEN-native, MATLAB-readable, and gzipped to
    match the tables (the YAML I/O is gzip-aware on a `.gz` suffix). Per-organism GPRs are built at
    runtime (3b.4/3b.5) from the KO↔reaction + organism gene↔KO tables.
  - **Relational tables = minimal, stored as gzipped TSV** — store *only* what 3b.4/3b.5/HMM-build
    consume: `ko_reaction` (KO↔reaction), `organism_gene_ko` (the large one; for 3b.4), `phyl_dist`
    (for 3b.5 weighting), KO names, and reaction-quality flags (spontaneous/incomplete/general/
    undefined-stoich, for the `keep*` filters). Nothing the functions don't use. **Format = gzipped
    TSV** (`.tsv.gz`), partitioned per organism for the large `organism_gene_ko` table. This
    supersedes the earlier SQLite/Parquet idea: gzipped TSV is the **dependency-free cross-language
    format** — pandas reads/writes it with **no extra package** (`read_csv`/`to_csv`,
    `compression="gzip"` built in) and MATLAB reads it natively (`readtable`, no toolbox). SQLite
    would need MATLAB's Database Toolbox; Parquet would need `pyarrow`/`fastparquet` on the Python
    side. See [docs/kegg_data_format.md](docs/kegg_data_format.md) for the rationale and the future
    options (Parquet/SQLite) we may revisit if a table grows large enough to justify the dependency.
  - **Distribution ✅ (`ensure_data`):** the gene-free reference GEM, the gzipped-TSV tables, and the
    prok90/euk90 **HMM library** are **separate, version-pinned downloads** — fetched on first use,
    SHA256-verified, cached under `~/.cache/ravengem/data/kegg-<version>/` (platformdirs), via the
    **`ensure_data`** registry in `data.py` mirroring `binaries.py`'s `ensure_binary`
    (`ensure_kegg_data` for the core set, `ensure_kegg_hmm_library` for a domain library). **Not
    bundled in the pip wheel** (size). The `…_from_artefacts` entry points fetch automatically when
    no local dir is given. Registry empty until the artefacts are published (same as the binary
    registry). Redistribution is licensed (above).

**Improvements to log:** split the overloaded `getKEGGModelForOrganism` (organism-vs-FASTA modes,
~15 params) into the two clear entry points **3b.4 / 3b.5**; ship version-pinned ravengem KEGG
artefacts (reference model + HMMs) so end users need no KEGG access; HMMER/MAFFT/CD-HIT via the
same version-pinned binary registry as 3a.

*Two audiences:* **3b.1–3b.3 are maintainer build-time tools** (paid FTP, run once per KEGG
release) producing published artefacts; **3b.4/3b.5 are the end-user runtime API** consuming those
artefacts. *Build order:* 3b.2 (parse — testable against a tiny dump fixture) → 3b.4 (annotation
mode, no tools) → 3b.3 (HMM build) → 3b.5 (HMM query). The model-building core (3b.2/3b.4) is
testable against small fixtures with no external tools; HMMER steps sit behind `skipif`.

#### 2.3c MetaCyc-based — **DROPPED (not ported)**
**Decision (2026-05-24): MetaCyc reconstruction will not be ported**, and is slated
for **removal from MATLAB RAVEN** too. See IMPROVEMENTS.md (R-MetaCyc) for the
empirical justification and the list of MATLAB functions/data to remove.

Rationale in brief: `getMetaCycModelForOrganism` calls genes by BLASTing the query
proteome against MetaCyc's **single representative sequence per enzyme** (~11.6k
seqs), with no profile to discriminate family members. A leave-organism-out
precision/recall test on real KEGG 118 / MetaCyc data showed this is unreliable at
**every** cutoff: at RAVEN's default (bitscore 100, ppos 45) only **36 % reaction-level
precision** (~64 % of assignments wrong) and **59 % even at EC-family level**;
tightening to bitscore 300 reaches only ~44 %/65 % precision while recall halves.
No cutoff makes it usable, and real proteomes (with non-enzyme decoys) would fare
worse. Accurate gene-calling already comes from KEGG HMMs (3b) and homology (3a);
MetaCyc's database value (extra reactions/pathways) does not justify a separate,
low-precision, data-heavy track. Its `addSpontaneousRxns`/reconciliation ideas can
be revisited as small standalone helpers if a concrete need arises.

### 2.4 `tasks/` — metabolic task validation  *(Phase 4a — done ✅)*
The foundation the INIT phases build on. No cobrapy equivalent.
| RAVEN | Notes |
|---|---|
| `parseTaskList` | ✅ `parse_task_list` + `Task` ([tasks/tasklist.py](src/ravengem/tasks/tasklist.py)) — TSV/xlsx, multi-row tasks, `;`-split, defaults, ALLMETS/ALLMETSIN. |
| `checkTasks` | ✅ `check_tasks` + `TaskResult` ([tasks/check.py](src/ravengem/tasks/check.py)) — inputs/outputs via relaxed metabolite mass-balance bounds (RAVEN's `b`), task equations, bound changes, closed boundaries; feasibility verdict (handles `should_fail`). |
| `fitTasks`, essential-reaction output (`getEssential`) | Deferred to **4c** — only tINIT needs them. |
| `checkProduction`, `getExpressionStructure` | Production checks underpinning tasks (port if needed). |

### 2.5 `init/` — tINIT (original INIT MILP)  *(Phase 4c — done ✅)*
RAVEN-unique MILP; no cobrapy equivalent. Needs a MIP solver. Depends on tasks (4a).
| RAVEN | Notes |
|---|---|
| `runINIT` | ✅ `run_init` ([init/init.py](src/ravengem/init/init.py)) — clean optlang reformulation of the INIT MILP. ⚠️ RAVEN's scale-dependent magic numbers (`eps`, `prod_weight`) exposed for tuning. |
| `scoreComplexModel` | ✅ `score_reactions_from_genes` + `gene_scores_from_expression` ([init/score.py](src/ravengem/init/score.py)) — GPR scoring + the common RNA-seq `5·ln(level/ref)` gene scoring. |
| `getINITModel` | ✅ `get_init_model` ([init/build.py](src/ravengem/init/build.py)) — scores → dead-end removal → `run_init`. **Deferred:** HPA/single-cell data ingestion → Phase 5 omics; automatic task-essential discovery + task gap-filling → 4d (shared with ftINIT). RNA-seq is the common input; single-cell/HPA are alternative upstream sources. |
| `removeLowScoreGenes`, `mergeLinear`, `rescaleModelForINIT`, `reverseRxns` | RAVEN MILP preprocessing/cleanup — port if profiling shows they matter (clean reformulation makes some moot). |

### 2.6 `init/` — ftINIT (fast staged INIT)  *(Phase 4d — CRITICAL REVIEW)*
> ⚠️ **ftINIT needs a lot of special attention.** It is the most complex algorithm in
> RAVEN — a multi-step MILP with task-aware gap-filling. **Review the MATLAB code very
> critically before porting**: do not transcribe blindly. Understand each step, question
> the formulation, check for bugs/edge cases, and validate against RAVEN outputs on real
> models. Likely the largest single port in the project.

| RAVEN | Notes |
|---|---|
| `ftINIT` | Top-level staged context-model extraction (the newer, faster INIT). |
| `prepINITModel`, `getINITSteps`, `ftINITInternalAlg`, `INITStepDesc` | The staged algorithm + step descriptors. |
| `ftINITFillGaps`, `ftINITFillGapsMILP`, `ftINITFillGapsForAllTasks` | Task-aware gap-filling within INIT. |

### 2.7 `gapfilling/` — connectivity gap-filling  *(Phase 4b — done ✅)*
Implemented as `connect_blocked_reactions` ([gapfilling/fill.py](src/ravengem/gapfilling/fill.py)):
MILP to add the fewest template reactions so blocked draft reactions carry flux. RAVEN's
objective-feasibility mode → `cobra.flux_analysis.gapfill` (§1 cheatsheet). Remaining
production/consumption diagnostics (`canProduce`/`canConsume`/`makeSomething`/`gapReport`)
are small follow-ups if needed.

### 2.8 `omics/` + `analysis/` + `comparison/` — data integration & analysis  *(Phase 5)*
| RAVEN | Notes |
|---|---|
| `reporterMetabolites` | ✅ `reporter_metabolites` ([analysis/reporter.py](src/ravengem/analysis/reporter.py)) — exact closed-form background replaces RAVEN's Monte-Carlo (RM1). |
| `parseHPA`, `parseHPArna`, `scoreModel` | ✅ `parse_hpa` / `parse_hpa_rna` / `hpa_gene_scores` / `rna_gene_scores` ([omics/hpa.py](src/ravengem/omics/hpa.py)) — pandas-tidy DataFrames replace RAVEN's sparse-matrix + cell-array layout; scoring adapters reuse the existing `score_reactions_from_genes` so the GPR walk has one source of truth. |
| `FSEOF` | ✅ `fseof` ([analysis/fseof.py](src/ravengem/analysis/fseof.py)) — redesigned output: regression slope+correlation, amplify/knockdown/knockout classes, gene aggregation (FS1–FS4). |
| `compareMultipleModels` | ✅ `compare_models` ([comparison/compare.py](src/ravengem/comparison/compare.py)) — returns `ModelComparison` of tidy DataFrames (reactions/metabolites/genes/subsystems presence + Jaccard similarity + optional `check_tasks` pass/fail). Plotting and tSNE/MDS deliberately not ported — one-liners in seaborn / scikit-learn on the returned DataFrames; keeping plotting out of the core function keeps it useful in pipelines. |
| `runDynamicFBA` (dynamic FBA) | **DO NOT PORT** — established Python implementations cover this: [`dfba`](https://pypi.org/project/dfba/) (Pinheiro et al., 2021; CVODES-backed), [`reframed`](https://pypi.org/project/reframed/) (Machado lab), [`mewpy`](https://pypi.org/project/mewpy/) (Cunha lab). Cobrapy itself has no dFBA. Re-porting `runDynamicFBA` would duplicate well-maintained prior art with no obvious value-add; users should reach for one of these. Same call as MetaCyc / gap-fill targeted mode (use established Python tooling, document the migration). |

### 2.9 `localization/` — subcellular localization  *(Phase 7, self-contained)*
A self-contained track (depends only on Phase 1, can be done anytime): predict subcellular
localization of gene products and compartmentalize a single-compartment model accordingly. No
cobra equivalent.
- **Entry point:** `predict_localization` (`predictLocalization`) — assign reactions to compartments
  from per-gene localization scores, iteratively moving reactions to minimise cross-membrane transport.
- **Pluggable predictors (not just WoLF PSORT):** the algorithm consumes a generic **gene→compartment
  score table**, so any localization tool can feed it. Provide loaders for **WoLF PSORT**
  (`getWoLFScores`) and modern predictors such as **DeepLoc** (and leave the table format open for
  others, e.g. TargetP/LocTree). The compartmentalization algorithm itself is predictor-agnostic.

| RAVEN | Notes |
|---|---|
| `predictLocalization` | ✅ `predict_localization` ([localization/predict.py](src/ravengem/localization/predict.py)) — **deterministic MILP** (not simulated annealing), **caller-passed relocate set** (`reactions_to_relocate=[…]`; everything else pinned), **incomplete-model-tolerant** (no silent reaction removal), `apply=False` returns a `LocalizationProposal` diff, **multi-compartment-by-default scoring** (primary free, extras pay `multi_compartment_penalty` *plus* their lower predictor score is the implicit penalty — no hard cutoff). Existing compartmentalisation respected by default. See [docs/localization_design.md](docs/localization_design.md). |
| `getWoLFScores` | ✅ `load_wolfpsort` ([localization/scores.py](src/ravengem/localization/scores.py)) — parses WoLF PSORT summary output (RAVEN-compatible); row-normalised to max=1.0. Does *not* call the WoLF PSORT binary (RAVEN's `getWoLFScores` shells out to Perl); run that separately and feed in the output. |
| *(new)* DeepLoc loader | ✅ `load_deeploc` ([localization/scores.py](src/ravengem/localization/scores.py)) — parses DeepLoc 2 per-protein CSV (Protein_ID, Localizations, Signals, then one column per compartment). |
| `mergeCompartments` | ✅ `merge_compartments` ([manipulation/compartments.py](src/ravengem/manipulation/compartments.py)) — collapses all metabolite compartments into one, deduplicates resulting identical reactions, optionally drops one-metabolite collapses (RAVEN's `deleteRxnsWithOneMet`). Useful **independently of `predict_localization`** for flattening before an analysis that doesn't care about topology. |
| `copyToComps` | ✅ `copy_to_compartment` ([manipulation/compartments.py](src/ravengem/manipulation/compartments.py)) — duplicates a set of reactions into a target compartment (mirror a pathway into mitochondria/peroxisomes). Idempotent; optional `delete_original=True` turns the copy into a move. |
| `mapCompartments` | **NOT PORTED** — its "transfer compartment assignments from a curated reference to a draft" use case overlaps with `compare_models` on the reaction id intersection; add a small adapter only if a real workflow needs it. |

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
| **3b** | Reconstruction — KEGG | 5-step pipeline (§2.3b): download KEGG → parse dump to reference model → build HMMs → model-for-species (annotation) → model-by-HMM-query (FASTA). | 1, 2 |
| ~~**3c**~~ | ~~Reconstruction — MetaCyc~~ | **DROPPED** (2026-05-24) — BLAST-to-single-representatives is low-precision at every cutoff (§2.3c, IMPROVEMENTS R-MetaCyc); also to be removed from MATLAB RAVEN. | — |
| **4a** | Metabolic tasks (the task file) | `tasks/` — `parseTaskList`, `checkTasks`/`fitTasks`. Foundation for the INIT phases. | 1, 2 |
| **4b** | Gap-filling | `gapfilling/` — `connect_blocked_reactions` ✅ (done). | 1, 2 |
| **4c** | tINIT | `init/` — original INIT MILP (`getINITModel`/`runINIT`) + reaction scoring. | 1, 2, 4a, MIP solver |
| **4d** | ftINIT | `init/` — fast staged INIT (`ftINIT` + task-aware gap-filling). **⚠️ critical review of the MATLAB code required; most complex port.** | 1, 2, 4a, 4c, MIP solver |
| **5** | Data integration & analysis | HPA/omics scoring, `reporterMetabolites`, FSEOF, dFBA, model comparison. | 1–4 |
| **6** | Visualization | pathway maps / omics overlay (consider Escher). | 1–2 |
| **7** | Localization | `localization/` — `predictLocalization` + pluggable predictors (WoLF PSORT, DeepLoc, …). Self-contained. | 1 |

**Suggested order rationale:** each phase produces something usable on its own. Reconstruction
(Phase 3) is RAVEN's headline feature and only needs the foundation + I/O. It has two
**independent** tracks — 3a homology and 3b KEGG (both done). The planned 3c MetaCyc track was
**dropped** (§2.3c): its homology gene-calling is low-precision and its database value doesn't
justify the track. The INIT phases (4c tINIT, 4d ftINIT) depend on the task framework (4a), so
tasks are built first; ftINIT (4d) is split out as its own phase because it is the most complex
algorithm in RAVEN and needs a critical, non-transcriptive port. Localization (7) is self-contained
(only needs Phase 1) and can be slotted in anytime.

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
