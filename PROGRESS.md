# ravengem — Refactoring Progress

Living status tracker for the RAVEN (MATLAB) → ravengem (Python/cobrapy) port.
**PLAN.md** is the design spec (what to port, why, and how it maps to cobrapy);
this file tracks **how far along** the port is. Update it whenever code lands.

_Last updated: 2026-05-22_

---

## Status at a glance

| Phase | Theme | Status |
|---|---|---|
| 0 | Scaffold & decisions | ✅ done |
| 1 | Foundation (`utils/`, `manipulation/`) | 🟡 in progress |
| 2 | I/O (`io/`) | ⬜ not started |
| 3 | Reconstruction (`reconstruction/`) | ⬜ not started |
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

**Test status:** 34 tests passing (incl. smoke) under cobra 0.31.1, run via geckopy's `.venv`.

---

## Subpackage scaffold

All subpackages exist as importable stubs (purpose docstring only) unless noted above.

| subpackage | purpose | port status |
|---|---|---|
| `utils/` | lookup + GPR hygiene + model helpers (`getIndexes`, `standardizeGrRules`, `checkModelStruct`, MIRIAM/annotation, ID-prefix) — **no** struct adapter | ⬜ stub |
| `manipulation/` | model construction, editing & structural transforms (ergonomic layer, see PLAN §1b) | 🟡 2 functions ported |
| `io/` | RAVEN YAML/Excel/SIF formats | ⬜ stub |
| `reconstruction/{kegg,metacyc,homology}/` | de novo reconstruction (flagship) | ⬜ stub |
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
   auto-create dependencies cobra leaves to you. PORT keystones: `getIndexes`, `addRxns` (+ family),
   `setParam`, `removeReactions`/`removeMets`/`removeGenes`, `simplifyModel`, `mergeModels`,
   `standardizeGrRules`, `sortModel` (deterministic core). These are a thin layer **on top of**
   cobra primitives, not a parallel data model.

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

---

## Next up

Candidate next steps (Phase 1–2), in rough priority order:

1. **`io/` YAML reader/writer** — high-use, self-contained, clear oracle (round-trip a
   yeast-GEM/Human-GEM file), gives immediate ecModel interop with geckopy.
2. **`utils/` foundation** — `checkModelStruct` validation + MIRIAM/annotation + ID-prefix helpers.
3. **More `manipulation/` transforms** — `simplifyModel`, `contractModel`, `mergeCompartments`,
   `copyToComps` (the parts cobra's `prune_*` helpers don't cover).
