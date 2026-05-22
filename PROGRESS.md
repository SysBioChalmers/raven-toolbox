# ravengem — Refactoring Progress

Living status tracker for the RAVEN (MATLAB) → ravengem (Python/cobrapy) port.
**PLAN.md** is the design spec (what to port, why, and how it maps to cobrapy);
this file tracks **how far along** the port is. Update it whenever code lands.

> This is **not** a one-to-one transcription. When porting a function, also judge whether it can be
> made smarter/faster or whether RAVEN is missing something that fits — and log those in
> **[IMPROVEMENTS.md](IMPROVEMENTS.md)** (candidates to also back-port to MATLAB RAVEN).

_Last updated: 2026-05-22_

---

## Status at a glance

| Phase | Theme | Status |
|---|---|---|
| 0 | Scaffold & decisions | ✅ done |
| 1 | Foundation (`utils/`, `manipulation/`) | 🟡 in progress |
| 2 | I/O (`io/`) | ⬜ not started |
| 3a | Reconstruction — homology (`reconstruction/homology/`) | ⬜ not started |
| 3b | Reconstruction — KEGG (`reconstruction/kegg/`) | ⬜ not started |
| 3c | Reconstruction — MetaCyc (`reconstruction/metacyc/`) | ⬜ not started |
| 4 | Context-specific & tasks (`tasks/`, `gapfilling/`, `init/`) | ⬜ not started |
| 5 | Data integration & analysis (`omics/`, `localization/`, `analysis/`, `comparison/`) | ⬜ not started |
| 6 | Visualization (`plotting/`) | ⬜ not started |

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
| `utils/compartments.py` | `get_metabolites_in_compartment`, `get_reactions_in_compartment` | `getMetsInComp.m`, `getRxnsInComp.m` | ✅ `tests/test_utils_compartments.py` | "Objects in compartment" selectors over cobra's `reaction.compartments`, with the `include_partial` (fully-contained vs touching) distinction. Return objects, not masks. |
| `io/yaml.py` | `read_yaml_model`, `write_yaml_model` | `readYAMLmodel.m` / `writeYAMLmodel.m` | ✅ `tests/test_io_yaml.py` | Wraps cobra YAML (which already reads the `!!omap` RAVEN/Human-GEM format) and adds what cobra drops: model identity/provenance from `metaData`, RAVEN-only per-entry fields routed by meaning (chemical IDs `smiles`/`inchis` → `annotation`; `deltaG`/`confidence_score`/`*From`/`protein` → `notes`), and verbatim preservation of foreign sections (GECKO ec). Verified on a real yeast-GEM.yml. |
| `manipulation/remove.py` | `remove_metabolites`, `remove_genes` | `removeMets.m` / `removeGenes.m` | ✅ `tests/test_manipulation_remove.py` | Delegate to cobra; add the gaps: `by_name` cross-compartment deletion (mets — flagged as a deletion candidate if unused), and a `blocked_reactions` remove/constrain/keep policy for gene knockouts (genes). `removeReactions` **not** ported (coupled orphan cleanup = cobra's `remove_reactions`). |

**Test status:** 107 tests passing (incl. smoke) under cobra 0.31.1, run via geckopy's `.venv`.

---

## Subpackage scaffold

All subpackages exist as importable stubs (purpose docstring only) unless noted above.

| subpackage | purpose | port status |
|---|---|---|
| `utils/` | GPR hygiene + compartment selectors + model helpers (`is_dnf`/`find_non_dnf_grrules` ✅, `get_*_in_compartment` ✅, `parse_name_comp`, `checkModelStruct`, MIRIAM/annotation, ID-prefix) — **no** struct adapter | 🟡 GPR lint + comp selectors ported |
| `manipulation/` | model construction, editing & structural transforms (ergonomic layer, see PLAN §1b) | 🟡 `add_reactions_from_equations` + 2 transforms ported |
| `io/` | RAVEN YAML/Excel/SIF formats | 🟡 YAML read/write ported |
| `reconstruction/homology/` | homology-based draft from a template GEM + BLAST/DIAMOND (3a) | ⬜ stub |
| `reconstruction/kegg/` | KEGG-based draft (orthology/KO assignment) (3b) | ⬜ stub |
| `reconstruction/metacyc/` | MetaCyc-based draft + KEGG reconciliation (3c) | ⬜ stub |
| `init/` | tINIT/ftINIT context extraction | ⬜ stub |
| `tasks/` | metabolic task validation | ⬜ stub |
| `gapfilling/` | template-based MILP gap-filling | ⬜ stub |
| `omics/` | HPA omics → reaction scores | ⬜ stub |
| `localization/` | WoLF PSORT subcellular localization | ⬜ stub |
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

---

## Next up

Candidate next steps, in rough priority order:

1. **More light functions** — `setParam` (batch bounds/objective + var-band), `addTransport`,
   `setExchangeBounds`, `sortModel` (deterministic core), `getElementalBalance`.
2. **`utils/` foundation** — `checkModelStruct` validation + MIRIAM/annotation + ID-prefix helpers.
3. **`mergeModels`** — N-way merge with name+comp matching, conflict rename, provenance (heavier).
4. **`simplifyModel`** — stage by mode (pure-graph modes first; FVA/groupLinear later).
5. **`io/` Excel / tab-delimited / SIF** exporters.
