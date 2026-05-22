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
