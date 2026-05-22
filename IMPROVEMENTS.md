# Proposed improvements over RAVEN

This project is **not** a one-to-one MATLAB→Python transcription. Where a RAVEN function can be
made smarter/faster, or where a logical gap in RAVEN's feature set is worth filling, we record the
change here — with enough detail that it can **also be back-ported to MATLAB RAVEN** later.

Each entry states: what RAVEN does today, the proposed improvement, the rationale, and whether it
is a candidate to upstream into MATLAB RAVEN.

Categories:
- **EFFICIENCY** — same behavior, faster/smarter implementation.
- **ERGONOMICS** — same job, less friction / fewer foot-guns / clearer contract.
- **NEW** — functionality RAVEN lacks but that fits naturally alongside what it already has.

Status legend: 💡 proposed · 🔨 implemented in ravengem · ⬆️ upstreamed to MATLAB RAVEN

---

## addRxns

RAVEN `core/addRxns.m` — add reactions from equation strings (or mets+coeffs), auto-creating
metabolites/genes. Ported as `add_reactions_from_equations`
([manipulation/add.py](src/ravengem/manipulation/add.py)).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| A1 | ERGONOMICS | ravengem 🔨 + MATLAB RAVEN 💡 | 🔨 | **Readable matching mode instead of the `eqnType` integer.** RAVEN takes `eqnType=1/2/3` (match by id / by name in a given compartment / `name[comp]`), which is opaque at call sites. ravengem uses `mets_by="id"\|"name"` and auto-detects `name[comp]` per token. **MATLAB back-port:** accept a string keyword. |
| A2 | ERGONOMICS (bug-class) | ravengem 🔨 | 🔨 | **Error on duplicate reaction IDs explicitly.** RAVEN errors; cobra's `add_reactions` *silently ignores* a duplicate. ravengem keeps RAVEN's stricter behaviour (raise) rather than cobra's silent drop. |
| A3 | EFFICIENCY (reuse) | ravengem 🔨 | 🔨 | **Delegate equation/arrow/coefficient parsing and gene/met creation to cobra** (`build_reaction_from_string` semantics, GPR auto-creation) instead of re-implementing RAVEN's `constructS`/`addGenesRaven`. Only the genuinely cobra-absent pieces (name matching, compartment for new mets, strict policies) are hand-written. |
| A4 | NEW | both 💡 | 💡 | **Infer compartment from a structured metabolite ID** (e.g. `atp_c` → `c`) as an alternative to requiring `compartment`. Not yet implemented; would reduce boilerplate for SBML-style IDs. Revisit alongside `addMets`. |

## changeGrRules

Ported as `change_gene_reaction_rules` ([manipulation/change.py](src/ravengem/manipulation/change.py)).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| G8 | EFFICIENCY (reuse) | ravengem 🔨 | 🔨 | **changeGrRules: delegate gene creation + normalization to cobra.** RAVEN calls `getGenesFromGrRules` + `addGenesRaven` + `standardizeGrRules` + rebuilds `rxnGeneMat`; cobra does all of that automatically on `gene_reaction_rule =`. The port keeps only the batch loop and the append (`(old) or (new)`) option. |

## setParam / getElementalBalance

Ported as `set_parameters` ([manipulation/parameters.py](src/ravengem/manipulation/parameters.py))
and `get_elemental_balance` ([utils/balance.py](src/ravengem/utils/balance.py)).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| ~~P1~~ | ERGONOMICS | — | ↩ revised | A 6-mode keyword `set_parameters` was built then **trimmed** (review: not Pythonic — it re-wrapped cobra one-liners for `lb`/`ub`/`eq`/`obj`/`unc`). Only the `var` ±% band, which cobra has no idiom for, is kept as `set_variance_bounds`; the rest are documented as cobra idioms in the §1 cheatsheet. |
| B1 | ERGONOMICS (correctness) | ravengem 🔨 | 🔨 | **getElementalBalance: report `unknown` for missing formulas.** cobra's `check_mass_balance` silently treats a metabolite with no formula as contributing nothing, so the reaction can read as (un)balanced on incomplete data. ravengem flags those as `unknown` rather than guessing — preserving RAVEN's distinction (its `-1` status). |

## getRxnsInComp / getMetsInComp — not ported

Briefly ported, then **removed** (user review): too thin over cobra (`metabolite.compartment` /
`reaction.compartments` one-liners). Mapped in the §1 migration cheatsheet instead. Reconsider only
if a downstream consumer needs the `include_partial` (fully-contained vs touching) distinction in
several places — and ask before re-adding (see process note: argue pros/cons for marginal WRAPs).

Ported as `remove_metabolites` / `remove_genes`
([manipulation/remove.py](src/ravengem/manipulation/remove.py)). `removeReactions` was **not**
ported: with orphan cleanup kept coupled (decision: don't separate metabolites from genes), it is
identical to `cobra.Model.remove_reactions(remove_orphans=...)`.

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| ~~R1~~ | ERGONOMICS | — | ❌ rejected | Separable orphan-met vs orphan-gene cleanup was considered, then **dropped** by decision — keep them coupled like cobra. (Removed the `remove_reactions` wrapper entirely as a result.) |
| R2 | EFFICIENCY (reuse) | ravengem 🔨 | 🔨 | **GPR rewriting delegated to cobra's AST**, not RAVEN's `eval` of a `&&`/`\|\|`-substituted rule string. cobra's `remove_genes` already gives correct AND/OR semantics (removing one gene of `A and B` empties the rule; of `A or B` keeps the other). **MATLAB back-port:** replace `canRxnCarryFlux`'s `eval` with a parsed boolean tree (safer, no eval). |
| R3 | ERGONOMICS | ravengem 🔨 | 🔨 | **`blocked_reactions` policy as a clear keyword** (`remove`/`constrain`/`keep`) instead of RAVEN's `removeBlockedRxns` boolean — and `keep` (rewrite GPR, leave bounds) is a third option RAVEN lacks. |
| R4 | (review) | ravengem ⚠️ | 💡 | **`remove_metabolites` is a deletion candidate.** Its only value over cobra is `by_name` cross-compartment deletion, likely rarely used; revisit and possibly drop the wrapper. |

## readYAMLmodel / writeYAMLmodel

RAVEN `io/readYAMLmodel.m` + `writeYAMLmodel.m` (+ private legacy parser). Ported as
`read_yaml_model`/`write_yaml_model` ([io/yaml.py](src/ravengem/io/yaml.py)).

**Lens correction (no separate legacy parser).** RAVEN ships a 462-line `parseYAMLLegacy.m` for the
`!!omap` dialect, and geckopy refuses it ("re-save from MATLAB"). But `!!omap` is **cobra's own YAML
format**: `cobra.io.load_yaml_model` reads a real yeast-GEM.yml (4102 rxns) directly. So the
ravengem-unique capability the PLAN imagined (a legacy reader) is unnecessary; the real cobra-absent
value is preserving `metaData` identity and RAVEN-only per-entry fields, which is what was built.

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| Y1 | EFFICIENCY (scope) | ravengem 🔨 | 🔨 | **Drop the bespoke legacy YAML parser; delegate to cobra's `!!omap` loader.** ~460 lines of RAVEN parsing not reimplemented. ravengem reads old RAVEN/Human-GEM YAML in pure Python with no MATLAB needed (geckopy can't — it tells users to re-save from MATLAB). |
| Y2 | ERGONOMICS (data loss) | ravengem 🔨 | 🔨 | **Don't silently drop model identity/provenance or RAVEN-only fields.** A plain `cobra.io.load_yaml_model` of a RAVEN file yields `model.id is None` and discards `smiles`/`deltaG`/`confidence_score`/etc. ravengem preserves them. **Routed by meaning** (not blindly to notes): chemical-structure identifiers `smiles`/`inchis` → cobra `annotation` (the MIRIAM-style store other tools read); genuinely non-standard data (`deltaG`, `confidence_score`, `metFrom`/`rxnFrom`, `protein`) → `notes`. Not invented as attributes (`met.deltaG`), since cobra only persists `annotation`/`notes` through copy/SBML/JSON/YAML. |
| Y4 | NEW | both 💡 | 💡 | **Upstream candidate: a first-class thermodynamics/confidence field.** `deltaG` and `confidence_score` live in `notes` because neither cobra nor SBML core has a home; if a standard slot (e.g. SBML fbc/groups or a cobra attribute) emerges, migrate them there. Also applies to MATLAB RAVEN's `metDeltaG`/`rxnConfidenceScores` consistency. |
| Y3 | NEW | both 💡 | 💡 | **Emit `!!omap`-tagged output matching cobra/Metabolic-Atlas exactly** for byte-stable diffs and guaranteed `cobra.io.load_yaml_model` interop. Current writer dumps a plain ruamel mapping (still re-readable, and round-trips through ravengem); aligning the exact tag/key order is a later refinement. |

## changeRxns

RAVEN `core/changeRxns.m` — change reaction equations. Ported as
`change_reaction_equations` ([manipulation/change.py](src/ravengem/manipulation/change.py)).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| C1 | EFFICIENCY | ravengem 🔨 | 🔨 | **Edit the reaction in place** instead of RAVEN's copy-all-fields → `removeReactions` → `addRxns` → `permuteModel` round-trip. cobra mutates the same `Reaction` object, so every other field and the model order are preserved for free, with no O(n) re-sort. (Not a MATLAB back-port — the round-trip is inherent to the struct layout there.) |

## getIndexes

RAVEN `core/getIndexes.m` — resolve a list of IDs / logical mask / index vector into positional
indexes (or a logical array) for `rxns` / `mets` / `genes` / `metNames` / `metcomps` (and GECKO
`ec.*` fields).

**Decision (ravengem): do NOT port the function.** cobra is object-oriented, so the central
index-resolver that RAVEN's struct-of-parallel-arrays design requires is largely unnecessary.
cobra's `DictList` already covers the use cases more idiomatically — `get_by_any` (mixed
id/object/index → objects), `get_by_id` (O(1)), `query` (name/substring/regex), `index` (position),
list comprehensions for filtering. Porting a 1-based-index resolver would be redundant and
un-Pythonic. **Only** the `name[comp]` composite resolver is kept (G7), as a small internal helper.

The improvement insights below still hold for **MATLAB RAVEN**, where the function remains — flagged
as upstream-only back-port candidates.

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| G1 | EFFICIENCY | MATLAB RAVEN only | 💡 | **Hash-based lookup instead of per-query linear scan.** RAVEN loops `find(strcmp(obj(i), searchIn))` per query → O(n·m). Build a `containers.Map` `{id: position}` once (O(n)), look up in O(1); `metcomps` likewise. (Moot for ravengem — cobra's `DictList` is already hashed.) |
| G5 | ERGONOMICS (bug) | MATLAB RAVEN only | 💡 | **Disambiguate the `[1 1 1]` mask-vs-index bug.** RAVEN's `if all(objects)` conflates a logical all-true mask with the index vector `[1 1 1]` (its own comment: "This gets weird if it's all 1"). Test `islogical(objects)` explicitly instead of `all(objects)`. (Moot for ravengem — input kind is decided by dtype.) |
| G7 | NEW | ravengem helper | 💡 | **Extract a reusable `name[comp]` parser/resolver.** The composite-id parsing buried in `getIndexes`'s `metcomps` branch is the one capability cobra lacks. Expose as a standalone `parse_name_comp` / `resolve_metabolite_by_name_comp`, reused by `addRxns`/`addTransport`/`mergeModels`. This is the only piece carried into ravengem. |

**Obsoleted by cobra (no action — these were earlier ravengem proposals now covered by `DictList`):**
predictable return type, return objects-not-positions, configurable missing-object policy across a
batch, and substring/regex matching — all already provided by `get_by_any` / `get_by_id` / `query`.

---

## standardizeGrRules

RAVEN `core/standardizeGrRules.m` — normalize grRule syntax + flag rules not in simple
OR-of-AND-complex (DNF) form (`findPotentialErrors`).

**Decision (ravengem): port the lint half only.** cobra auto-normalizes a GPR on assignment
(`"(G1 AND G2)  OR  G3"` is stored as `"(G1 and G2) or G3"`), so the normalization half is
redundant. The non-DNF lint has no cobra equivalent and was ported as `find_non_dnf_grrules`/`is_dnf`
([utils/gpr.py](src/ravengem/utils/gpr.py)).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| S1 | ERGONOMICS | ravengem 🔨 + MATLAB RAVEN 💡 | 🔨 | **Return structured lint results, don't just print.** RAVEN's `findPotentialErrors` only emits a `warning()` string; you can't act on it programmatically. ravengem returns a list of `GPRIssue(reaction_id, gpr, reason)`. **MATLAB back-port:** return the `indexes2check`/messages as a struct array (it already computes `indexes2check` — just surface it cleanly instead of only warning). |
| S2 | EFFICIENCY (robustness) | ravengem 🔨 + MATLAB RAVEN 💡 | 🔨 | **Detect non-DNF via the boolean AST, not substring search.** RAVEN scans for the `) and (`, `) and`, `and (` substrings, which is brittle (sensitive to spacing/bracketing and to gene IDs containing those characters). ravengem walks cobra's GPR AST (`is_dnf`: no OR beneath any AND), which is exact. **MATLAB back-port:** parse the rule to a tree rather than string-matching. |
