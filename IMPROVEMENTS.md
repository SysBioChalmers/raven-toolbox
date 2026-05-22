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
