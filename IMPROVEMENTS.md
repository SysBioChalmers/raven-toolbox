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
