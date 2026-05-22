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

| # | Cat | Status | Improvement |
|---|---|---|---|
| G1 | EFFICIENCY | 💡 | **Hash-based lookup instead of per-query linear scan.** RAVEN loops `find(strcmp(obj(i), searchIn))` for each query → O(n·m). Build a `{id: position}` map once (O(n)) and look up in O(1); the `metcomps` path likewise gets a `{name[comp]: position}` map. In ravengem this is free — cobra's `DictList.get_by_id` is already hashed; the port should lean on it rather than re-scan. **MATLAB back-port:** build a `containers.Map` once at the top. |
| G2 | ERGONOMICS | 💡 | **Single, predictable return type.** RAVEN returns a double array, *except* `metnames` returns a cell array, and a length-1 cell is auto-unwrapped — callers can't predict the shape. ravengem: always return the same container type for a given call; `name`-type lookups (1→many) return a list-per-query consistently, never a bare scalar. |
| G3 | ERGONOMICS | 💡 | **Configurable missing-object policy** (`on_missing="raise"\|"warn"\|"skip"\|None`). RAVEN always hard-errors on a missing object, forcing try/catch around bulk resolution. Default stays `raise` (preserve behavior); add the softer modes. **MATLAB back-port:** an optional `onMissing` arg. |
| G4 | ERGONOMICS | 💡 | **Return objects, not just positions.** In cobra you almost always want the `Reaction`/`Metabolite`/`Gene` object, not a 1-based index. ravengem's primary resolver returns objects (or `None` for missing under the soft policy); a thin `..._index()` variant returns positions only when genuinely needed. This makes the function Pythonic rather than a port of MATLAB's index-centric idiom. |
| G5 | ERGONOMICS | 💡 | **Disambiguate the `[1 1 1]` mask-vs-index bug.** RAVEN's `if all(objects)` branch conflates a logical all-true mask with the numeric index vector `[1 1 1]` (its own comment: "This gets weird if it's all 1"). Decide input kind by **dtype** (bool mask vs integer positions), not by value, so an index vector of all-ones resolves correctly. **MATLAB back-port:** test `islogical(objects)` explicitly instead of `all(objects)`. |
| G6 | NEW | 💡 | **Optional substring/regex/case-insensitive matching** (e.g. `match="exact"\|"icontains"\|"regex"`), mirroring cobra `DictList.query`. RAVEN only does exact `strcmp`. Useful for interactive exploration and name-based curation; default stays exact. |
| G7 | NEW | 💡 | **Generalize `metcomps` to a reusable `name[comp]` parser.** The composite-id parsing is buried inside `getIndexes`; expose it as a small standalone helper (`parse_name_comp`) reused by `addRxns`/`addTransport`/`mergeModels`, which all do the same name+compartment matching. |

**Net design for ravengem:** a `get_indexes(model, objects, type, *, return_logical=False, on_missing="raise", match="exact")` that resolves objects via cobra's hashed `DictList`, returns objects by default (positions on request), with predictable shapes and a configurable missing policy. The enduring value over plain cobra is the **tolerant mixed-input resolution** (id / name / index / mask / `name[comp]`) and the **uniform contract**, not the indexing per se.
