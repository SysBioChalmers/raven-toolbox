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
- **REMOVAL** — functionality that should be dropped (here *and* in MATLAB RAVEN) because it does
  more harm than good.

Status legend: 💡 proposed · 🔨 implemented in ravengem · ⬆️ upstreamed to MATLAB RAVEN ·
🗑️ dropped (and to remove from MATLAB RAVEN)

---

## R-MetaCyc: drop MetaCyc reconstruction (REMOVAL — also remove from MATLAB RAVEN)

**Decision 2026-05-24:** MetaCyc-based reconstruction is **not ported to ravengem** and should be
**removed from MATLAB RAVEN**. Status: 🗑️.

**What RAVEN does:** `getMetaCycModelForOrganism` builds a draft by BLAST/DIAMOND of the query
proteome against `protseq.fsa` — MetaCyc's **single representative protein sequence per enzyme**
(~11.6k sequences) — keeping each gene's best hit above a bitscore/positives cutoff and assigning
the linked reaction. With one representative per enzyme there is no profile to tell true family
members from look-alikes.

**Evidence (this repo, real MetaCyc + KEGG 118 data):** a leave-organism-out precision/recall test
(query each representative against the others, excluding its own organism; ground truth =
MetaCyc's own MONOMER→reaction):

| bitscore (ppos≥45) | reaction precision | EC-family precision | EC recall |
|---|---|---|---|
| 50 | 0.34 | 0.55 | 0.33 |
| 100 *(RAVEN default)* | 0.36 | 0.59 | 0.32 |
| 200 | 0.40 | 0.62 | 0.26 |
| 300 | 0.44 | 0.65 | 0.22 |

At the default cutoff **~64 % of reaction assignments are wrong** (~41 % wrong even at EC-family
level); **no cutoff rescues precision** — tightening to bitscore 300 reaches only ~44 %/65 % while
recall halves. Real proteomes (with non-enzyme decoys, not in this test) would be worse. Test
scripts/artifacts: `/home/eduardk/metacyc_test/` (not committed).

**Why drop rather than fix:** the low precision is intrinsic to MetaCyc's one-representative-per-
enzyme data (can't build KEGG-quality HMMs from it). Accurate gene-calling already exists via KEGG
HMMs (3b) and homology-to-template-models (3a). MetaCyc's genuine value (extra reactions/pathways/
compound structures) does not justify a separate, data-heavy, low-precision track.

**MATLAB RAVEN removal list** (`external/metacyc/`): `getMetaCycModelForOrganism.m`,
`getModelFromMetaCyc.m`, `getRxnsFromMetaCyc.m`, `getMetsFromMetaCyc.m`, `getEnzymesFromMetaCyc.m`,
`linkMetaCycKEGGRxns.m`, `addSpontaneousRxns.m`, and data `metaCycEnzymes.mat` / `metaCycMets.mat`
/ `metaCycRxns.mat` / `protseq.fsa`; plus any `combineMetaCycKEGGModels` and MetaCyc references in
tutorials/tests/docs. (`addSpontaneousRxns` could be kept as a small standalone helper if wanted —
it is only incidentally in the MetaCyc folder.)

---

## getModelFromHomology (Phase 3a — implemented)

Design + rationale in [docs/plan_get_model_from_homology.md](docs/plan_get_model_from_homology.md);
implemented in `reconstruction/homology/homology.py`. *Logic* improvements over RAVEN's algorithm
(RAVEN's own comments flag several of these spots as uncertain).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| H1 | ERGONOMICS | ravengem 🔨 + MATLAB RAVEN 💡 | 🔨 | Split the overloaded `strictness` 1/2/3 into two orthogonal params: `bidirectional` (reciprocal hits) and `best_hits_only`. RBH = both true. `strictness=` kept as a compat alias. |
| H2 | EFFICIENCY (robustness) | ravengem 🔨 + MATLAB RAVEN 💡 | 🔨 | Rewrite GPRs on the **cobra GPR AST**, not `regexprep` string substitution — eliminates partial-match hazards and the `OLD_… or` regex cleanup pass RAVEN needs. |
| H3 | ERGONOMICS (correctness) | ravengem 🔨 | 🔨 | Explicit `complex_policy` (default **`flag`** = RAVEN-compatible `OLD_`; plus `keep`/`drop`) for AND-subunits lacking an ortholog, via correct OR/AND **AST** semantics. |
| H4 | (correctness) | both 🔨/💡 | 🔨 | Best-hit selection by **bitscore** (db-size-independent, the RBH standard); `score="evalue"` optional. |
| H5 | EFFICIENCY | ravengem 🔨 | 🔨 | DataFrame ortholog map (pandas merge + dict) replaces `allGenes`/`allTo`/`allFrom` sparse-matrix `sub2ind` index juggling. |
| H6 | NEW | ravengem 🔨 | 🔨 | Structured provenance: `HomologyResult.gene_map` + per-reaction `notes['homology_source']`. |

## KEGG download / dump parsing / HMM build (Phase 3b.1 / 3b.2 / 3b.3 — implemented)

`fetch_keggdb.sh` → `reconstruction/kegg/download.py` (3b.1); parsing core of
`getRxnsFromKEGG` / `getMetsFromKEGG` / `getGenesFromKEGG` / `getModelFromKEGG`
→ `reconstruction/kegg/parse.py` (3b.2); `constructMultiFasta` + the
cluster/align/train stages of `getKEGGModelForOrganism` → `reconstruction/kegg/hmm.py`
and `taxonomy.py` (3b.3). Maintainer-side, build-time tooling (PLAN.md §2.3b).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| K1 | EFFICIENCY (robustness) | ravengem 🔨 + MATLAB RAVEN 💡 | 🔨 | Read each reaction's equation from its **own `EQUATION` field**, not from `reaction.lst` matched by line order. RAVEN reads `reaction.lst` line *i* into reaction *i*, assuming the two files stay perfectly aligned — brittle. **MATLAB back-port:** parse the `EQUATION` field already present in `reaction`. |
| K2 | ERGONOMICS (correctness) | ravengem 🔨 + MATLAB RAVEN 💡 | 🔨 | Undefined-stoichiometry terms (`n C00001`, `(n+1) C00002`) keep their **real compound id** with coefficient 1 and the reaction is *flagged*, instead of minting `"n C00001"` pseudo-metabolites later renamed `undefined_N`. Cleaner metabolite graph; flag still drives the `keep*` filters. |
| K3 | ERGONOMICS | ravengem 🔨 + MATLAB RAVEN 💡 | 🔨 | Reaction quality labels become a tidy boolean **`rxn_flags` table** (spontaneous/undefined-stoich/incomplete/general) instead of free-text appended to `rxnNotes`, so downstream filters join on a column rather than substring-matching notes. |
| K4 | EFFICIENCY | ravengem 🔨 | 🔨 | **Gene-free reference model** + separate `organism_gene_ko` table (the big one), instead of RAVEN's giant `rxnGeneMat` baked into the global model. Per-organism GPRs are built only at runtime (3b.4/3b.5), keeping the published artefact small. |
| K5 | EFFICIENCY (portability) | ravengem 🔨 | 🔨 | **KEGG download in pure Python stdlib** (`urllib`/`tarfile`/`gzip`/`netrc`), porting `fetch_keggdb.sh`. Drops the script's `wget`/`tar`/`gunzip` (and Cygwin-on-Windows) requirement, so it runs unchanged on Linux/macOS/Windows; tar extraction uses the `data` filter (no path traversal); same `~/.netrc` credential hygiene. The arrange step is split out (`extract_kegg_dump`) so it's network-free and unit-tested. |
| K6 | EFFICIENCY | ravengem 🔨 + MATLAB RAVEN 💡 | 🔨 | **Per-KO multi-FASTA via a stdlib offset index** (`_index_fasta` → seek), replacing `constructMultiFasta`'s Java-`Hashtable` byte scan with 5M-element preallocation. One streaming pass, only wanted ids retained; no MATLAB/Java heap tuning. |
| K7 | EFFICIENCY | ravengem 🔨 + MATLAB RAVEN 💡 | 🔨 | **Concatenate per-KO HMMs and `hmmpress` into one pressed library**, so the query path (3b.5) runs a single `hmmscan` against the database instead of RAVEN's thousands of per-KO `hmmsearch` invocations. |
| K8 | EFFICIENCY (scope) | ravengem 🔨 | 🔨 | **Drop the `getPhylDist` distance matrix.** Its only uses in RAVEN were per-organism HMM-sequence subsampling (`maxPhylDist`/`nSequences`) and the kingdom filter. Our fixed prok90/euk90 libraries (3b.3) remove the subsampling rationale, and domain mode (3b.4) uses the taxonomy domain classification directly — so the O(n²) matrix is never built. Simpler, faster, less code. |
| K9 | EFFICIENCY (memory) | ravengem 🔨 | 🔨 | **Stream `organism_gene_ko` straight to gzipped TSV** in `parse_kegg_dump` instead of building it in memory. Real KEGG has **9.05M** gene↔KO associations; the in-memory DataFrame build OOMs in a few GB. Streaming runs the full parse in **82 s at 0.9 GB peak**. (Found by validating against a real KEGG FTP dump.) |
| K10 | EFFICIENCY (size) | ravengem 🔨 | 🔨 | **Reference model as gzipped RAVEN/cobra YAML** (`reference_model.yml.gz`) rather than SBML: RAVEN-native, MATLAB-readable, and ~1.1 MB vs ~30 MB SBML for the real 12k-reaction model. Made `io/yaml.py` gzip-aware on a `.gz` suffix (general-purpose). |
| K11 | ERGONOMICS | ravengem 🔨 | 🔨 | **`ensure_data`** (`data.py`): version-pinned registry that fetches/verifies/caches the published KEGG artefacts under `~/.cache/ravengem/data/`, mirroring `ensure_binary`. End users get a draft model with no KEGG access and no manual data handling — the `…_from_artefacts` entry points auto-fetch when no local dir is supplied. |
| K12 | EFFICIENCY | ravengem 🔨 + MATLAB RAVEN 💡 | 🔨 | **Fast MAFFT (FFT-NS-2) for HMM training** instead of RAVEN's `--auto`, which selects slow iterative refinement (`dvtditr`) on medium/large KOs — observed ~2.5 min/KO (days for a domain) on real KEGG 118. FFT-NS-2 (`--retree 2 --maxiterate 0`) is seconds/KO and ample for profile-HMM building. **PartTree cutover is residue-based and memory-auto-tuned**: MAFFT memory tracks residues (count × length), not sequence count, so a count threshold let long-protein KOs (K00901: 2,788 seqs, 2.55 M residues) OOM under FFT-NS-2 — measured ~5 GB MAFFT RSS with FFT-NS-2 vs **0.69 GB with PartTree** for the same alignment. The cutover uses residues only, and the budget is **derived from available RAM by inverting an empirically-measured FFT-NS-2 memory model**: peak RSS is super-linear, `RSS_GB ≈ 1.32·R² + 1.84·R` (R = M residues; measured on real KEGG sequences: 0.25/0.5/1.0/1.5 M → 0.67/1.25/3.16/5.73 GB). `_auto_residue_budget` solves that for the residue count fitting in 0.7 × (total − 2.5 GB overhead) — ~1.09 M on a 7.6 GB box, ~2.1 M @ 16 GB, ~5 M @ 64 GB — and **warns on low-memory hosts**. (A naive linear estimate gave 1.9 M here, which would have OOM'd.) Override via `parttree_residues` / `--parttree-residues`. Back-portable to RAVEN. |
| K13 | EFFICIENCY | ravengem 🗑️ | 🗑️ | ~~Per-KO sequence cap (`max_sequences`)~~ — **removed.** Briefly added as a count-based cap, but the residue-based PartTree cutover (K12) bounds MAFFT memory without dropping any sequences, so the cap was redundant complexity. All deduplicated sequences are kept. |

## FSEOF (Phase 5 — implemented, redesigned)

RAVEN `core/FSEOF.m` → `analysis/fseof.py` (`fseof`). User was unhappy with RAVEN's
output; redesigned substantially.

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| FS1 | CORRECTNESS | ravengem 🔨 + MATLAB RAVEN 💡 | 🔨 | **Robust trend via linear regression** (slope + correlation) over the whole scan, instead of RAVEN's strict step-by-step monotonicity that discards a target on a single noisy step (LP alternative optima). |r| is a quality score for filtering. pFBA per step keeps the scan stable. |
| FS2 | NEW | ravengem 🔨 + MATLAB RAVEN 💡 | 🔨 | **Reports knockdown/knockout targets**, not just amplification. RAVEN only flags reactions whose \|flux\| *rises* with the enforced product; reactions driven *toward zero* — the down-regulation/deletion candidates, arguably the most actionable — are classified here (`knockdown`/`knockout`). |
| FS3 | ERGONOMICS | ravengem 🔨 | 🔨 | **Gene-level aggregation** (`gene_targets`) mapping reaction targets to genes, plus the **full flux scan** retained — all as DataFrames, vs RAVEN's printed TSV + endpoint-only slope. |
| FS4 | CORRECTNESS | ravengem 🔨 | 🔨 | Slope is the **regression slope** consistent with the selection criterion, not RAVEN's endpoint difference that disagreed with its own monotonicity filter. |

## reporterMetabolites (Phase 5 — implemented)

RAVEN `core/reporterMetabolites.m` → `analysis/reporter.py` (`reporter_metabolites`).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| RM1 | EFFICIENCY + CORRECTNESS | ravengem 🔨 + MATLAB RAVEN 💡 | 🔨 | **Exact closed-form background correction** instead of RAVEN's 100 000-random-set Monte Carlo *per neighbour-count*. RAVEN samples with replacement from the scored-gene Z pool, so a random aggregate `Σz/√n` provably has mean `√n·μ` and std `σ` — the corrected score is exactly `(metZ − √n·μ)/σ`. Removes the slow sampling **and** its run-to-run randomness (deterministic results); back-portable to RAVEN. |
| RM2 | ERGONOMICS | ravengem 🔨 | 🔨 | Returns a sorted **DataFrame** per test (`all`/`up`/`down`) and takes gene→p-value / gene→fold-change **dicts**, vs RAVEN's parallel arrays + struct array + print/file side-effects. Neighbour genes come from cobra's metabolite→reaction→gene graph (no `rxnGeneMat`). |

## runINIT (Phase 4c — MILP core implemented)

RAVEN `INIT/runINIT.m` → `init/init.py` (`run_init`).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| I1 | ERGONOMICS | ravengem 🔨 | 🔨 | **Clean optlang reformulation** of the INIT MILP instead of RAVEN's hand-built sparse `prob.A`/`blc`/`buc`/`vartype` arrays + fake "FAKEFORPM" metabolites. Standard include-indicator form `eps·x ≤ v ≤ ub·x` with objective `max Σ score·x + prod_weight·Σ sink`. Far more readable/reviewable; functional equivalence is the bar (PLAN §0). |
| I2 | ERGONOMICS | ravengem 🔨 | 🔨 | **`no_rev_loops` as a single `x_fwd + x_rev ≤ 1`** per reversible reaction, replacing RAVEN's auxiliary A/B/C metabolites with int1/int2 reactions and `C ub=-1` construction. Same effect (no spurious forward/back connectivity loop), a fraction of the machinery. |
| I3 | ERGONOMICS | ravengem 🔨 | 🔨 | **`present_mets` producibility via a small LP feasibility test** (sum of compartment-form drains ≥ 1), instead of mutating the live MILP's RHS one metabolite at a time. |
| I4 | CORRECTNESS | ravengem 🔨 + MATLAB RAVEN 💡 | 🔨 | **MILP big-M is each reaction's own `ub`** (`v ≤ ub·x`), not RAVEN's fixed 1000; and `eps`/`prod_weight` are exposed parameters. RAVEN's hard-coded 1/1000/0.1/0.5 only suit ±1000-bounded models with O(1) scores — flagged as scale-dependent and tunable (don't blindly trust them). |
| I5 | ERGONOMICS | ravengem 🔨 | 🔨 | **Predictor-agnostic scoring**: `get_init_model` takes gene *or* reaction scores; gene scoring is generic (`gene_scores_from_expression` for the common RNA-seq path), so single-cell/HPA are just alternative upstream sources feeding the same gene→score table — rather than RAVEN's HPA/array-specific structs baked into `getINITModel`. |

## parseTaskList / checkTasks (Phase 4a — implemented)

RAVEN `core/parseTaskList.m` + `core/checkTasks.m` → `tasks/tasklist.py` + `tasks/check.py`.

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| T1 | ERGONOMICS | ravengem 🔨 | 🔨 | **Structured `Task` dataclass + `TaskResult`** instead of RAVEN's parallel-array struct and a printed report; programmatic access to per-task pass/fail/feasibility/error. |
| T2 | ERGONOMICS | ravengem 🔨 | 🔨 | **TSV-first task files** (stdlib `csv`); `.xlsx` still supported but via the lazy `[excel]` extra — no hard Excel dependency just to read a task list. |
| T3 | EFFICIENCY | ravengem 🔨 | 🔨 | Inputs/outputs imposed directly on cobra's **metabolite mass-balance constraint bounds** (the analogue of RAVEN's two-column `model.b`), and existing boundary reactions are auto-closed — so a model with open exchanges is handled correctly (RAVEN assumes a closed model and silently misbehaves otherwise). |

## fillGaps (Phase 4b — implemented)

RAVEN `core/fillGaps.m`. Only the **connectivity** mode is ported, as
`connect_blocked_reactions` ([gapfilling/fill.py](src/ravengem/gapfilling/fill.py)) —
MILP via cobra/optlang (GLPK). RAVEN's other mode (fill to make the objective feasible)
is `cobra.flux_analysis.gapfill` and is **cheatsheeted, not re-wrapped** (PLAN §1).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| GF1 | NEW (vs cobra) | ravengem 🔨 | 🔨 | **Connectivity gap-fill has no cobra equivalent**: add the minimum-penalty set of template reactions so *blocked* draft reactions can carry flux (cobra's `gapfill` only fills toward the objective). Ported as `connect_blocked_reactions` — a name that avoids confusion with `cobra.gapfill` and says what it does, vs RAVEN's overloaded `fillGaps(useModelConstraints=...)` boolean. |
| GF2 | ERGONOMICS | ravengem 🔨 | 🔨 | **Templates matched by `name[comp]`** (via `add_reactions_from_model`), so a template in a different identifier namespace than the draft still contributes — as RAVEN's name-based merge does. (For the targeted `cobra.gapfill` path, ids must be aligned first, since cobra matches by id — noted in the cheatsheet.) |

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

## simplifyModel

Gap modes ported in [manipulation/simplify.py](src/ravengem/manipulation/simplify.py).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| S3 | EFFICIENCY (scope) | ravengem 🔨 | 🔨 | **Only the cobra-absent modes ported as focused functions**, not a monolithic 8-flag `simplifyModel`. `deleteMinMax`→`find_blocked_reactions`, `deleteZeroInterval`→filter+prune, `deleteUnconstrained`→moot are cheatsheeted. dead-end / duplicate / constrain-reversible / group-linear are standalone, composable functions. |

## mergeModels

Ported as `merge_models` ([manipulation/merge.py](src/ravengem/manipulation/merge.py)).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| M1 | EFFICIENCY (scope) | ravengem 🔨 | 🔨 | **~560 lines of struct field-padding + manual S-matrix assembly dropped.** On `cobra.Model` the merge is just: unify metabolites by `name[comp]`, re-add reactions remapped to the merged metabolites, let cobra rebuild S and create genes. |
| M2 | ERGONOMICS | ravengem 🔨 | 🔨 | **Provenance via `notes['origin']`** (one place) instead of three parallel `rxnFrom`/`metFrom`/`geneFrom` fields. `match_by="name"|"id"` keyword replaces RAVEN's `metParam` string. |

## checkModelStruct

Ported (curation subset) as `check_model` ([utils/validate.py](src/ravengem/utils/validate.py)).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| V1 | EFFICIENCY (scope) | ravengem 🔨 | 🔨 | **Drop the struct/type/duplicate-ID/`lb>ub`/`rev` checks** — cobra's object model enforces or precludes them (DictList forbids duplicate IDs, Reaction rejects `lb>ub`, no `rev` field). Only the curation checks cobra lacks survive. |
| V2 | ERGONOMICS | ravengem 🔨 + MATLAB RAVEN 💡 | 🔨 | **Return structured `ModelIssue`s, not printed warnings** (RAVEN prints / throws). Programmatically filterable by `category`. **MATLAB back-port:** return an issues struct array. |

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
| Y3 | NEW | ravengem 🔨 | 🔨 | **Emit cobra-native `!!omap` output** (via cobra's own dumper) — done, matching RAVEN `fa281a1`. Verified `cobra.io.load_yaml_model` reads the output. |
| Y5 | ERGONOMICS (correctness) | ravengem 🔨 | 🔨 | **Field placement realigned to `fa281a1`:** `smiles`/`ec-code` are in the cobra-owned `annotation` block (not top-level), `inchis` is top-level, and the top-level `notes` *string* (metNotes/rxnNotes) is handled rather than crashing a notes-as-dict assumption. |

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
