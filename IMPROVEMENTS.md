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

Status legend: 💡 proposed · 🔨 implemented in raven-python · ⬆️ upstreamed to MATLAB RAVEN ·
🗑️ dropped (and to remove from MATLAB RAVEN)

---

## R-MetaCyc: drop MetaCyc reconstruction (REMOVAL — also remove from MATLAB RAVEN)

**Decision 2026-05-24:** MetaCyc-based reconstruction is **not ported to raven-python** and should be
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

Design + rationale in [docs/plan_get_model_from_homology.md](https://github.com/SysBioChalmers/raven-python/blob/develop/docs/archive/plan_get_model_from_homology.md);
implemented in `reconstruction/homology/homology.py`. *Logic* improvements over RAVEN's algorithm
(RAVEN's own comments flag several of these spots as uncertain).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| H1 | ERGONOMICS | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | Split the overloaded `strictness` 1/2/3 into two orthogonal params: `bidirectional` (reciprocal hits) and `best_hits_only`. RBH = both true. `strictness=` kept as a compat alias. |
| H2 | EFFICIENCY (robustness) | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | Rewrite GPRs on the **cobra GPR AST**, not `regexprep` string substitution — eliminates partial-match hazards and the `OLD_… or` regex cleanup pass RAVEN needs. |
| H3 | ERGONOMICS (correctness) | raven-python 🔨 | 🔨 | Explicit `complex_policy` (default **`flag`** = RAVEN-compatible `OLD_`; plus `keep`/`drop`) for AND-subunits lacking an ortholog, via correct OR/AND **AST** semantics. |
| H4 | (correctness) | both 🔨/💡 | 🔨 | Best-hit selection by **bitscore** (db-size-independent, the RBH standard); `score="evalue"` optional. |
| H5 | EFFICIENCY | raven-python 🔨 | 🔨 | DataFrame ortholog map (pandas merge + dict) replaces `allGenes`/`allTo`/`allFrom` sparse-matrix `sub2ind` index juggling. |
| H6 | NEW | raven-python 🔨 | 🔨 | Structured provenance: `HomologyResult.gene_map` + per-reaction `notes['homology_source']`. |

## KEGG download / dump parsing / HMM build (Phase 3b.1 / 3b.2 / 3b.3 — implemented)

`fetch_keggdb.sh` → `reconstruction/kegg/download.py` (3b.1); parsing core of
`getRxnsFromKEGG` / `getMetsFromKEGG` / `getGenesFromKEGG` / `getModelFromKEGG`
→ `reconstruction/kegg/parse.py` (3b.2); `constructMultiFasta` + the
cluster/align/train stages of `getKEGGModelForOrganism` → `reconstruction/kegg/hmm.py`
and `taxonomy.py` (3b.3). Maintainer-side, build-time tooling (PLAN.md §2.3b).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| K1 | EFFICIENCY (robustness) | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | Read each reaction's equation from its **own `EQUATION` field**, not from `reaction.lst` matched by line order. RAVEN reads `reaction.lst` line *i* into reaction *i*, assuming the two files stay perfectly aligned — brittle. **MATLAB back-port:** parse the `EQUATION` field already present in `reaction`. |
| K2 | ERGONOMICS (correctness) | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | Undefined-stoichiometry terms (`n C00001`, `(n+1) C00002`) keep their **real compound id** with coefficient 1 and the reaction is *flagged*, instead of minting `"n C00001"` pseudo-metabolites later renamed `undefined_N`. Cleaner metabolite graph; flag still drives the `keep*` filters. |
| K3 | ERGONOMICS | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | Reaction quality labels become a tidy boolean **`rxn_flags` table** (spontaneous/undefined-stoich/incomplete/general) instead of free-text appended to `rxnNotes`, so downstream filters join on a column rather than substring-matching notes. |
| K4 | EFFICIENCY | raven-python 🔨 | 🔨 | **Gene-free reference model** + separate `organism_gene_ko` table (the big one), instead of RAVEN's giant `rxnGeneMat` baked into the global model. Per-organism GPRs are built only at runtime (3b.4/3b.5), keeping the published artefact small. |
| K5 | EFFICIENCY (portability) | raven-python 🔨 | 🔨 | **KEGG download in pure Python stdlib** (`urllib`/`tarfile`/`gzip`/`netrc`), porting `fetch_keggdb.sh`. Drops the script's `wget`/`tar`/`gunzip` (and Cygwin-on-Windows) requirement, so it runs unchanged on Linux/macOS/Windows; tar extraction uses the `data` filter (no path traversal); same `~/.netrc` credential hygiene. The arrange step is split out (`extract_kegg_dump`) so it's network-free and unit-tested. |
| K6 | EFFICIENCY | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **Per-KO multi-FASTA via a stdlib offset index** (`_index_fasta` → seek), replacing `constructMultiFasta`'s Java-`Hashtable` byte scan with 5M-element preallocation. One streaming pass, only wanted ids retained; no MATLAB/Java heap tuning. |
| K7 | EFFICIENCY | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **Concatenate per-KO HMMs into one library**, so the query path (3b.5) runs a single `hmmsearch` over the multi-profile library instead of RAVEN's thousands of per-KO `hmmsearch` invocations — same fast search direction, one process, no `hmmpress`. |
| K8 | EFFICIENCY (scope) | raven-python 🔨 | 🔨 | **Drop the `getPhylDist` distance matrix.** Its only uses in RAVEN were per-organism HMM-sequence subsampling (`maxPhylDist`/`nSequences`) and the kingdom filter. Our fixed prok90/euk90 libraries (3b.3) remove the subsampling rationale, and domain mode (3b.4) uses the taxonomy domain classification directly — so the *reconstruction* path never builds the O(n²) matrix. The matrix itself remains available on demand via `taxonomy.phyl_dist`, which regenerates RAVEN's `keggPhylDist` from the published `taxonomy` artefact for GECKO's organism-distance kcat selection (no `.mat` needed). |
| K9 | EFFICIENCY (memory) | raven-python 🔨 | 🔨 | **Stream `organism_gene_ko` to disk** in `parse_kegg_dump` instead of building it in memory. Real KEGG has **9.05M** gene↔KO associations; the in-memory DataFrame build OOMs in a few GB. Streaming (now via the external merge sort of K14) runs the full parse with flat, bounded peak memory. (Found by validating against a real KEGG FTP dump.) |
| K10 | EFFICIENCY (size) | raven-python 🔨 | 🔨 | **Reference model as gzipped RAVEN/cobra YAML** (`reference_model.yml.gz`) rather than SBML: RAVEN-native, MATLAB-readable, and ~1.1 MB vs ~30 MB SBML for the real 12k-reaction model. Made `io/yaml.py` gzip-aware on a `.gz` suffix (general-purpose). |
| K11 | ERGONOMICS | raven-python 🔨 | 🔨 | **`ensure_data`** (`data.py`): version-pinned registry that fetches/verifies/caches the published KEGG artefacts under `~/.cache/raven-python/data/`, mirroring `ensure_binary`. End users get a draft model with no KEGG access and no manual data handling — the `…_from_artefacts` entry points auto-fetch when no local dir is supplied. |
| K12 | EFFICIENCY | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **Fast MAFFT (FFT-NS-2) for HMM training** instead of RAVEN's `--auto`, which selects slow iterative refinement (`dvtditr`) on medium/large KOs — observed ~2.5 min/KO (days for a domain) on real KEGG 118. FFT-NS-2 (`--retree 2 --maxiterate 0`) is seconds/KO and ample for profile-HMM building. **PartTree cutover is residue-based and memory-auto-tuned**: MAFFT memory tracks residues (count × length), not sequence count, so a count threshold let long-protein KOs (K00901: 2,788 seqs, 2.55 M residues) OOM under FFT-NS-2 — measured ~5 GB MAFFT RSS with FFT-NS-2 vs **0.69 GB with PartTree** for the same alignment. The cutover is **length-aware and memory-auto-tuned**: FFT-NS-2 memory is driven by the progressive-alignment **DP cost ≈ n_seqs × mean_len²** (= residues²/n_seqs), *not* residue count — a few hundred long proteins cost far more than the same residues in many short ones. (First tried a residue-only model `RSS≈1.32R²+1.84R`; it then OOM'd on K12047 — 452 seqs but mean length 2082, 0.94 M residues — because long proteins blow the per-residue cost.) Calibrated `RSS_GB ≈ 4.2e-9 × (n_seqs × mean_len²)` across real KEGG KOs (250k/266→0.67 GB … 1.5M/1624→5.73 GB; K12047 cost 1.96e9 = the largest, hence its OOM). `_auto_cost_budget` switches to PartTree when the DP cost exceeds `0.65 × (total − 2.5 GB overhead) / 4.2e-9` (≈7.9e8 on a 7.6 GB box), **warns on low-memory hosts**, and `parttree_residues` overrides with a manual residue cutoff. Back-portable to RAVEN. |
| K13 | EFFICIENCY | raven-python 🗑️ | 🗑️ | ~~Per-KO sequence cap (`max_sequences`)~~ — **removed.** Briefly added as a count-based cap, but the residue-based PartTree cutover (K12) bounds MAFFT memory without dropping any sequences, so the cap was redundant complexity. All deduplicated sequences are kept. |
| K14 | EFFICIENCY (size) | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **Sort `organism_gene_ko` by `(organism, gene)` and store it gzipped** (`organism_gene_ko.tsv.gz`), cutting the dominant artefact **≈78 → 48 MB** by sorting alone. Gene IDs within an organism share long prefixes (locus tags, numeric runs), so sorting makes them adjacent and far more compressible. The sort is an **external merge sort** bounded to `chunk_rows` rows in memory (sorted runs spooled to gzipped temp files, merged with `heapq.merge`), so it keeps K9's flat memory profile. We first xz-compressed this file (≈27 MB, 2.9×) but switched to **gzip** (≈74 MB) so MATLAB reads it with built-in `gunzip` and no external `unxz`: the artefacts are shared with MATLAB RAVEN, and the once-per-release size cut wasn't worth a MATLAB toolchain dependency. All tables are now gzipped TSV, native on Windows/macOS/Linux. Sorted order also matches the by-organism query in `get_kegg_model_for_organism`, enabling a future `searchsorted` slice instead of loading all 9M rows. Back-portable to RAVEN. |
| K15 | ERGONOMICS (correctness) | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **Recalibrate the HMM-query KO-assignment defaults** (`assign_kos`): cut-off `1e-50 → 1e-30`, `min_score_ratio_g 0.8 → 0.9`; `min_score_ratio_ko` left at 0.3 but **documented as empirically inert**. Cross-validated the full 3b.5 pipeline against the true KEGG gene→KO annotation of four organisms across both libraries and the well-/lesser-studied axis — *S. cerevisiae*, *Cyanidioschyzon merolae* (red alga), *E. coli* K-12, *Mycoplasma genitalium* (minimal genome). Real annotations score overwhelmingly (median E ≈ 1e-100…1e-155; even the weakest 1% ≈ 1e-15…1e-36) while spurious hits cluster at ≈1e-8 — a ~20-order-of-magnitude gap. RAVEN's `1e-50` therefore sits **inside the true-positive tail** and silently drops real-but-divergent hits for no noise-rejection gain: gene→KO recall on *M. genitalium* was only 0.84 (reaction recall 0.87). At `1e-30` + `ratio_g=0.9`: *M. genitalium* recall **0.84→0.94** (rxn 0.87→0.97), *E. coli* 0.95→0.97 with **fewer** unannotated reactions (198→173, the tighter gene-ratio prunes spurious multi-KO genes), *S. cerevisiae*/*C. merolae* held or improved. The three sweep tables showed `min_score_ratio_ko` produced identical output at 0.0/0.3/0.5 across all four organisms — a magic-number knob that does nothing; `min_score_ratio_g` is the real precision lever. Full numbers in [docs/kegg_hmm_cutoff_calibration.md](https://github.com/SysBioChalmers/raven-python/blob/develop/docs/studies/kegg_hmm_cutoff_calibration.md) (reproduce with `scripts/analyze_hmm_cutoffs.py`). Back-portable to RAVEN. |

## FSEOF (Phase 5 — implemented, redesigned)

RAVEN `analysis/FSEOF.m` → `analysis/fseof.py` (`fseof`). User was unhappy with RAVEN's
output; redesigned substantially.

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| FS1 | CORRECTNESS | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **Robust trend via linear regression** (slope + correlation) over the whole scan, instead of RAVEN's strict step-by-step monotonicity that discards a target on a single noisy step (LP alternative optima). |r| is a quality score for filtering. pFBA per step keeps the scan stable. |
| FS2 | NEW | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **Reports knockdown/knockout targets**, not just amplification. RAVEN only flags reactions whose \|flux\| *rises* with the enforced product; reactions driven *toward zero* — the down-regulation/deletion candidates, arguably the most actionable — are classified here (`knockdown`/`knockout`). |
| FS3 | ERGONOMICS | raven-python 🔨 | 🔨 | **Gene-level aggregation** (`gene_targets`) mapping reaction targets to genes, plus the **full flux scan** retained — all as DataFrames, vs RAVEN's printed TSV + endpoint-only slope. |
| FS4 | CORRECTNESS | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | Slope is the **regression slope** consistent with the selection criterion, not RAVEN's endpoint difference that disagreed with its own monotonicity filter. **Proposed back-port:** record the enforced target-flux level at each iteration into a vector and replace both the reported slope `abs(fseof.results(num,iterations) - fseof.results(num,1)) / abs(targetMax - targetMax/iterations)` and the `targets.slope` field with `polyfit(enforcedFlux, fseof.results(num,:), 1)` slope across all iterations. The label classification (amplify/knockdown/knockout, see FS2) should also use the slope of `|flux|` instead of comparing endpoints. |

## randomSampling (Phase 5 — implemented)

RAVEN `analysis/randomSampling.m` → `analysis/sampling.py` (`random_sampling`). The
random-objective / extreme-point method of Bordel et al. (2010) — **not** what
`cobra.sampling` (OptGP/ACHR) does (those draw a near-uniform MCMC sample of the
polytope interior), so it is a genuine addition, and it was wrongly listed as
cobra-covered in the PLAN cheatsheet. Each sample maximises a small random linear
combination of reactions.

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| SAMP1 | EFFICIENCY | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **`good_reactions` (loop-free, flux-carrying objective candidates) via one `cobra` FVA pass** rather than RAVEN's hand-rolled per-reaction `parfor` that solves a separate LP maximising and minimising every reaction. FVA computes the same min/max ranges, optimised and optionally `loopless` (`cycleFreeFlux`), in far less code. |
| SAMP2 | ERGONOMICS | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **Reproducible** (`seed`) and **`n_objectives` is a parameter** — RAVEN has no seed and hard-codes 2 objective reactions (its docstring even claims 3, a transcription bug worth fixing upstream). |
| SAMP3 | ERGONOMICS | raven-python 🔨 | 🔨 | **Output is a samples × reactions DataFrame** (the `cobra.sampling` layout, directly usable with pandas/`analyzeSampling`-style stats) plus the reusable `good_reactions` list — instead of a reactions × samples matrix and a parallel `goodRxns` index vector that the caller must re-thread. |
| SAMP4 | CORRECTNESS | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **`replace_max_bound` (RAVEN's `replaceBoundsWithInf`) defaults off and is scoped to sampling only.** It applies *after* `good_reactions` is found (FVA cannot evaluate `inf` bounds — it errors as 'unbounded'), and it can open unbounded loop directions; documented to pair with `min_flux`. RAVEN ran goodRxns detection on the inf-replaced model, conflating "loop reaction" with "unbounded objective". |

## reporterMetabolites (Phase 5 — implemented)

RAVEN `analysis/reporterMetabolites.m` → `analysis/reporter.py` (`reporter_metabolites`).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| RM1 | EFFICIENCY + CORRECTNESS | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **Exact closed-form background correction** instead of RAVEN's 100 000-random-set Monte Carlo *per neighbour-count*. RAVEN samples with replacement from the scored-gene Z pool, so a random aggregate `Σz/√n` provably has mean `√n·μ` and std `σ` — the corrected score is exactly `(metZ − √n·μ)/σ`. Removes the slow sampling **and** its run-to-run randomness (deterministic results); back-portable to RAVEN. |
| RM2 | ERGONOMICS | raven-python 🔨 | 🔨 | Returns a sorted **DataFrame** per test (`all`/`up`/`down`) and takes gene→p-value / gene→fold-change **dicts**, vs RAVEN's parallel arrays + struct array + print/file side-effects. Neighbour genes come from cobra's metabolite→reaction→gene graph (no `rxnGeneMat`). |

## runINIT (Phase 4c — MILP core implemented)

RAVEN `INIT/runINIT.m` → `init/init.py` (`run_init`).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| I1 | ERGONOMICS | raven-python 🔨 | 🔨 | **Clean optlang reformulation** of the INIT MILP instead of RAVEN's hand-built sparse `prob.A`/`blc`/`buc`/`vartype` arrays + fake "FAKEFORPM" metabolites. Standard include-indicator form `eps·x ≤ v ≤ ub·x` with objective `max Σ score·x + prod_weight·Σ sink`. Far more readable/reviewable; functional equivalence is the bar (PLAN §0). |
| I2 | ERGONOMICS | raven-python 🔨 | 🔨 | **`no_rev_loops` as a single `x_fwd + x_rev ≤ 1`** per reversible reaction, replacing RAVEN's auxiliary A/B/C metabolites with int1/int2 reactions and `C ub=-1` construction. Same effect (no spurious forward/back connectivity loop), a fraction of the machinery. |
| I3 | ERGONOMICS | raven-python 🔨 | 🔨 | **`present_mets` producibility via a small LP feasibility test** (sum of compartment-form drains ≥ 1), instead of mutating the live MILP's RHS one metabolite at a time. |
| I4 | CORRECTNESS | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **MILP big-M is each reaction's own `ub`** (`v ≤ ub·x`), not RAVEN's fixed 1000; and `eps`/`prod_weight` are exposed parameters. RAVEN's hard-coded 1/1000/0.1/0.5 only suit ±1000-bounded models with O(1) scores — flagged as scale-dependent and tunable (don't blindly trust them). |
| I5 | ERGONOMICS | raven-python 🔨 | 🔨 | **Predictor-agnostic scoring**: `get_init_model` takes gene *or* reaction scores; gene scoring is generic (`gene_scores_from_expression` for the common RNA-seq path), so single-cell/HPA are just alternative upstream sources feeding the same gene→score table — rather than RAVEN's HPA/array-specific structs baked into `getINITModel`. |

## ftINIT (Phase 4d — in progress)

RAVEN `INIT/ftINITInternalAlg.m` (+ orchestration) → `init/ftinit.py` (`run_ftinit`),
`tasks/check.py` (`find_task_essential_reactions`). See docs/ftinit_review_and_plan.md.

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| FT1 | EFFICIENCY | raven-python 🔨 | 🔨 | **Essential-reaction discovery via one FVA pass** (`find_task_essential_reactions`) — a reaction is essential for a task iff its FVA range excludes 0 (= RAVEN's constrain-to-0→infeasible), restricted to the flux-carrying candidates of a `pfba` solution. Replaces `getEssentialRxns`' per-reaction knockout loop. |
| FT2 | ERGONOMICS | raven-python 🔨 | 🔨 | **Clean optlang reformulation** of the 6-category MILP instead of RAVEN's hand-built block `prob.a` with the `pi/ni/ei/vprb/vnrvm…` figure. Positive-score reactions keep the continuous-indicator trick (no binary — the ftINIT speedup), encoded as `net_flux ≥ force_on·y`; only negative scores and reversible-direction get binaries. |
| FT3 | CORRECTNESS | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **Big-M is each reaction's own bound**, not RAVEN's fixed 100/1000; `force_on`/`force_on_ess` exposed. Scale-dependent, calibrated in 4d.7 (cf. I4). |
| FT4 | (caveat) | RAVEN parity | — | **No loopless constraint** (matches RAVEN): the bare MILP can "include" an internal thermodynamically-infeasible cycle if it carries positive net score. Loop-free models rely on the staged pipeline + exchange handling and, at genome scale, real exchanges making cycles non-optimal. A loopless option could be added later. |
| FT9 | CORRECTNESS (robustness) | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **Per-reaction essential forcing, clamped to capacity** (4d.7 calibration). Staging fixes previous-step reactions as essential; RAVEN forces them at `min(0.99·\|prev flux\|, 0.1)` so none is forced above the flux it carried (`essential_force`). On top, the forced magnitude is **clamped to the reaction's bound** so a low-capacity essential never produces an `lb>ub` error — it is forced to its capacity instead (RAVEN would error). |
| FT-met | (deferred) | — | ⬜ | **Metabolomics production bonus (4d.6) deferred.** The linear merge eliminates degree-2 detected metabolites, so it needs RAVEN's producer-group-mapping + `mon`/`vnrbm`/`vnrvm`/`vnim` negative-producer force-flux block — the most intricate MILP in ftINIT, for its least-used input. `ftinit(metabolomics=…)` raises `NotImplementedError`. |
| FT8 | ERGONOMICS | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **`remove_low_score_genes`** (port of `removeLowScoreGenes`): prune negative-scoring genes from GPRs via a recursive walk of **cobra's GPR AST** — isozymes (OR) dropped keeping ≥1 (least-negative if all negative), complex subunits (AND) kept intact, nested groups handled — instead of RAVEN's regex subset-extraction with `#n#` placeholders. The all-negative tie-break is **deterministic** (reproducible) vs RAVEN's random. Wired into `ftinit(gene_scores=…)` (with `fill_tasks` for gap-filling) to complete the pipeline. Matches RAVEN's three docstring example cases. |
| FT7 | EFFICIENCY | raven-python 🔨 | 🔨 | **Task gap-filling** (`fill_tasks`, port of `ftINITFillGapsForAllTasks` + the fill MILP): per task, only if it is **infeasible** in the current model (a cheap LP gates the MILP), add the minimum-cost set of reference reactions (cost = `−min(score, −0.1)`) to satisfy the task's ranged metabolite bounds; accumulate across tasks; exchange reactions excluded as candidates. Built on the shared `apply_task_constraints` (factored out of `check_tasks`) + an optlang on/off MILP, instead of RAVEN's hand-built `[S pos neg int b var]` block and the custom `rxnScores`-field-on-the-model hack. Matches RAVEN tinitTests T0003 (gap at R7 → R7 added back). |
| FT6 | ERGONOMICS | raven-python 🔨 | 🔨 | **Staged pipeline** (`prep_init_model` → `PrepData`, `get_init_steps`/`InitStep`, `ftinit`): one-time reaction classification (`classify_reactions` = the `toIgnore` masks) + essential discovery + linear merge bundled into `PrepData`, reused per sample; the staged `'1+1'`/`'2+1'`/`full` schedule runs `run_ftinit` per step with previous-step reactions fixed as essential **in their flux direction** (carried as `essential_directions`, no model-flipping). RAVEN's Gurobi-specific per-step MIPGap retry schedule is dropped (solver-agnostic; tuning → 4d.7). Matches RAVEN tinitTests T0001 (no tasks → {R1,R4,R6,R8,R9,R10}; with R7/R10 spontaneous → {R1,R2,R4,R6,R7,R8}) and T0002 (task → essentials {R1,R7}, model {R1,R2,R4,R6,R7,R8,R9,R10}). |
| FT5 | EFFICIENCY | raven-python 🔨 | 🔨 | **Linear merge** (`merge_linear` + `group_rxn_scores`, port of `mergeLinear`/`groupRxnScores`): contracts degree-2-metabolite reaction chains (one producer, one consumer; reversibles included) into single reactions, losslessly shrinking the MILP (~⅓ fewer reactions on Human-GEM). Clean Python on a working representation with reversibility ≡ `lb<0` (RAVEN's `rev1·rev2` falls out of the most-constraining bound recompute); drops genes and objective on the reduced model (it exists only to feed the MILP, which scores via `group_rxn_scores`). Matches RAVEN's tinitTests T0004 exactly (group ids, merged bounds/reversibility, flipped reactions, grouped scores incl. the 0→0.01 handling). |

## parseTaskList / checkTasks (Phase 4a — implemented)

RAVEN `tasks/parseTaskList.m` + `tasks/checkTasks.m` → `tasks/tasklist.py` + `tasks/check.py`.

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| T1 | ERGONOMICS | raven-python 🔨 | 🔨 | **Structured `Task` dataclass + `TaskResult`** instead of RAVEN's parallel-array struct and a printed report; programmatic access to per-task pass/fail/feasibility/error. |
| T2 | ERGONOMICS | raven-python 🔨 | 🔨 | **TSV-first task files** (stdlib `csv`); `.xlsx` still supported but via the lazy `[excel]` extra — no hard Excel dependency just to read a task list. |
| T3 | EFFICIENCY | raven-python 🔨 | 🔨 | Inputs/outputs imposed directly on cobra's **metabolite mass-balance constraint bounds** (the analogue of RAVEN's two-column `model.b`), and existing boundary reactions are auto-closed — so a model with open exchanges is handled correctly (RAVEN assumes a closed model and silently misbehaves otherwise). |

## fillGaps (Phase 4b — implemented)

RAVEN `gapfilling/fillGaps.m`. Only the **connectivity** mode is ported, as
`connect_blocked_reactions` ([gapfilling/fill.py](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/gapfilling/fill.py)) —
MILP via cobra/optlang (GLPK). RAVEN's other mode (fill to make the objective feasible)
is `cobra.flux_analysis.gapfill` and is **cheatsheeted, not re-wrapped** (PLAN §1).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| GF1 | NEW (vs cobra) | raven-python 🔨 | 🔨 | **Connectivity gap-fill has no cobra equivalent**: add the minimum-penalty set of template reactions so *blocked* draft reactions can carry flux (cobra's `gapfill` only fills toward the objective). Ported as `connect_blocked_reactions` — a name that avoids confusion with `cobra.gapfill` and says what it does, vs RAVEN's overloaded `fillGaps(useModelConstraints=...)` boolean. |
| GF2 | ERGONOMICS | raven-python 🔨 | 🔨 | **Templates matched by `name[comp]`** (via `add_reactions_from_model`), so a template in a different identifier namespace than the draft still contributes — as RAVEN's name-based merge does. (For the targeted `cobra.gapfill` path, ids must be aligned first, since cobra matches by id — noted in the cheatsheet.) |

## addRxns

RAVEN `manipulation/addRxns.m` — add reactions from equation strings (or mets+coeffs), auto-creating
metabolites/genes. Ported as `add_reactions_from_equations`
([manipulation/add.py](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/manipulation/add.py)).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| A1 | ERGONOMICS | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **Readable matching mode instead of the `eqnType` integer.** RAVEN takes `eqnType=1/2/3` (match by id / by name in a given compartment / `name[comp]`), which is opaque at call sites. raven-python uses `mets_by="id"\|"name"` and auto-detects `name[comp]` per token. **MATLAB back-port:** accept a string keyword. |
| A2 | ERGONOMICS (bug-class) | raven-python 🔨 | 🔨 | **Error on duplicate reaction IDs explicitly.** RAVEN errors; cobra's `add_reactions` *silently ignores* a duplicate. raven-python keeps RAVEN's stricter behaviour (raise) rather than cobra's silent drop. |
| A3 | EFFICIENCY (reuse) | raven-python 🔨 | 🔨 | **Delegate equation/arrow/coefficient parsing and gene/met creation to cobra** (`build_reaction_from_string` semantics, GPR auto-creation) instead of re-implementing RAVEN's `constructS`/`addGenesRaven`. Only the genuinely cobra-absent pieces (name matching, compartment for new mets, strict policies) are hand-written. |
| A4 | NEW | both 💡 | 💡 | **Infer compartment from a structured metabolite ID** (e.g. `atp_c` → `c`) as an alternative to requiring `compartment`. Not yet implemented; would reduce boilerplate for SBML-style IDs. Revisit alongside `addMets`. |

## changeGrRules

Ported as `change_gene_reaction_rules` ([manipulation/change.py](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/manipulation/change.py)).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| G8 | EFFICIENCY (reuse) | raven-python 🔨 | 🔨 | **changeGrRules: delegate gene creation + normalization to cobra.** RAVEN calls `getGenesFromGrRules` + `addGenesRaven` + `standardizeGrRules` + rebuilds `rxnGeneMat`; cobra does all of that automatically on `gene_reaction_rule =`. The port keeps only the batch loop and the append (`(old) or (new)`) option. |

## simplifyModel

Gap modes ported in [manipulation/simplify.py](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/manipulation/simplify.py).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| S3 | EFFICIENCY (scope) | raven-python 🔨 | 🔨 | **Only the cobra-absent modes ported as focused functions**, not a monolithic 8-flag `simplifyModel`. `deleteMinMax`→`find_blocked_reactions`, `deleteZeroInterval`→filter+prune, `deleteUnconstrained`→moot are cheatsheeted. dead-end / duplicate / constrain-reversible / group-linear are standalone, composable functions. |

## mergeModels

Ported as `merge_models` ([manipulation/merge.py](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/manipulation/merge.py)).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| M1 | EFFICIENCY (scope) | raven-python 🔨 | 🔨 | **~560 lines of struct field-padding + manual S-matrix assembly dropped.** On `cobra.Model` the merge is just: unify metabolites by `name[comp]`, re-add reactions remapped to the merged metabolites, let cobra rebuild S and create genes. |
| M2 | ERGONOMICS | raven-python 🔨 | 🔨 | **Provenance via `notes['origin']`** (one place) instead of three parallel `rxnFrom`/`metFrom`/`geneFrom` fields. `match_by="name"|"id"` keyword replaces RAVEN's `metParam` string. |

## checkModelStruct

Ported (curation subset) as `check_model` ([utils/validate.py](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/utils/validate.py)).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| V1 | EFFICIENCY (scope) | raven-python 🔨 | 🔨 | **Drop the struct/type/duplicate-ID/`lb>ub`/`rev` checks** — cobra's object model enforces or precludes them (DictList forbids duplicate IDs, Reaction rejects `lb>ub`, no `rev` field). Only the curation checks cobra lacks survive. |
| V2 | ERGONOMICS | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **Return structured `ModelIssue`s, not printed warnings** (RAVEN prints / throws). Programmatically filterable by `category`. **MATLAB back-port:** return an issues struct array. |

## setParam / getElementalBalance

Ported as `set_parameters` ([manipulation/parameters.py](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/manipulation/parameters.py))
and `get_elemental_balance` ([utils/balance.py](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/utils/balance.py)).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| ~~P1~~ | ERGONOMICS | — | ↩ revised | A 6-mode keyword `set_parameters` was built then **trimmed** (review: not Pythonic — it re-wrapped cobra one-liners for `lb`/`ub`/`eq`/`obj`/`unc`). Only the `var` ±% band, which cobra has no idiom for, is kept as `set_variance_bounds`; the rest are documented as cobra idioms in the §1 cheatsheet. |
| B1 | ERGONOMICS (correctness) | raven-python 🔨 | 🔨 | **getElementalBalance: report `unknown` for missing formulas.** cobra's `check_mass_balance` silently treats a metabolite with no formula as contributing nothing, so the reaction can read as (un)balanced on incomplete data. raven-python flags those as `unknown` rather than guessing — preserving RAVEN's distinction (its `-1` status). |
| B2 | CORRECTNESS | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **Empty-stoichiometry reactions report `unknown`, not vacuous `balanced`.** A reaction with no metabolites used to fall through to `balanced` (any-over-empty is False; the zero-imbalance dict is empty). Same bug in MATLAB `queries/getElementalBalance.m`: with no entries in `model.S(:,j)`, `balanceStatus` stays NaN through both loops and the final `isnan→1` step labels it `balanced`. **Proposed back-port:** add an `emptyRxns = full(sum(model.S~=0,1))==0` mask and `balanceStatus(emptyRxns) = min(-1, balanceStatus(emptyRxns))` before the `isnan→1` line, so empty reactions are marked "missing information". |

## getRxnsInComp / getMetsInComp — not ported

Briefly ported, then **removed** (user review): too thin over cobra (`metabolite.compartment` /
`reaction.compartments` one-liners). Mapped in the §1 migration cheatsheet instead. Reconsider only
if a downstream consumer needs the `include_partial` (fully-contained vs touching) distinction in
several places — and ask before re-adding (see process note: argue pros/cons for marginal WRAPs).

Ported as `remove_metabolites` / `remove_genes`
([manipulation/remove.py](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/manipulation/remove.py)). `removeReactions` was **not**
ported: with orphan cleanup kept coupled (decision: don't separate metabolites from genes), it is
identical to `cobra.Model.remove_reactions(remove_orphans=...)`.

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| ~~R1~~ | ERGONOMICS | — | ❌ rejected | Separable orphan-met vs orphan-gene cleanup was considered, then **dropped** by decision — keep them coupled like cobra. (Removed the `remove_reactions` wrapper entirely as a result.) |
| R2 | EFFICIENCY (reuse) | raven-python 🔨 | 🔨 | **GPR rewriting delegated to cobra's AST**, not RAVEN's `eval` of a `&&`/`\|\|`-substituted rule string. cobra's `remove_genes` already gives correct AND/OR semantics (removing one gene of `A and B` empties the rule; of `A or B` keeps the other). **MATLAB back-port:** replace `canRxnCarryFlux`'s `eval` with a parsed boolean tree (safer, no eval). |
| R3 | ERGONOMICS | raven-python 🔨 | 🔨 | **`blocked_reactions` policy as a clear keyword** (`remove`/`constrain`/`keep`) instead of RAVEN's `removeBlockedRxns` boolean — and `keep` (rewrite GPR, leave bounds) is a third option RAVEN lacks. |
| R4 | (review) | raven-python ⚠️ | 💡 | **`remove_metabolites` is a deletion candidate.** Its only value over cobra is `by_name` cross-compartment deletion, likely rarely used; revisit and possibly drop the wrapper. |

## readYAMLmodel / writeYAMLmodel

RAVEN `io/readYAMLmodel.m` + `writeYAMLmodel.m` (+ private legacy parser). Ported as
`read_yaml_model`/`write_yaml_model` ([io/yaml.py](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/io/yaml.py)).

**Lens correction (no separate legacy parser).** RAVEN ships a 462-line `parseYAMLLegacy.m` for the
`!!omap` dialect, and geckopy refuses it ("re-save from MATLAB"). But `!!omap` is **cobra's own YAML
format**: `cobra.io.load_yaml_model` reads a real yeast-GEM.yml (4102 rxns) directly. So the
raven-python-unique capability the PLAN imagined (a legacy reader) is unnecessary; the real cobra-absent
value is preserving `metaData` identity and RAVEN-only per-entry fields, which is what was built.

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| Y1 | EFFICIENCY (scope) | raven-python 🔨 | 🔨 | **Drop the bespoke legacy YAML parser; delegate to cobra's `!!omap` loader.** ~460 lines of RAVEN parsing not reimplemented. raven-python reads old RAVEN/Human-GEM YAML in pure Python with no MATLAB needed (geckopy can't — it tells users to re-save from MATLAB). |
| Y2 | ERGONOMICS (data loss) | raven-python 🔨 | 🔨 | **Don't silently drop model identity/provenance or RAVEN-only fields.** A plain `cobra.io.load_yaml_model` of a RAVEN file yields `model.id is None` and discards `smiles`/`deltaG`/`confidence_score`/etc. raven-python preserves them. **Routed by meaning** (not blindly to notes): chemical-structure identifiers `smiles`/`inchis` → cobra `annotation` (the MIRIAM-style store other tools read); genuinely non-standard data (`deltaG`, `confidence_score`, `metFrom`/`rxnFrom`, `protein`) → `notes`. Not invented as attributes (`met.deltaG`), since cobra only persists `annotation`/`notes` through copy/SBML/JSON/YAML. |
| Y4 | NEW | both 💡 | 💡 | **Upstream candidate: a first-class thermodynamics/confidence field.** `deltaG` and `confidence_score` live in `notes` because neither cobra nor SBML core has a home; if a standard slot (e.g. SBML fbc/groups or a cobra attribute) emerges, migrate them there. Also applies to MATLAB RAVEN's `metDeltaG`/`rxnConfidenceScores` consistency. |
| Y3 | NEW | raven-python 🔨 | 🔨 | **Emit cobra-native `!!omap` output** (via cobra's own dumper) — done, matching RAVEN `fa281a1`. Verified `cobra.io.load_yaml_model` reads the output. |
| Y5 | ERGONOMICS (correctness) | raven-python 🔨 | 🔨 | **Field placement realigned to `fa281a1`:** `smiles`/`ec-code` are in the cobra-owned `annotation` block (not top-level), `inchis` is top-level, and the top-level `notes` *string* (metNotes/rxnNotes) is handled rather than crashing a notes-as-dict assumption. |

## changeRxns

RAVEN `manipulation/changeRxns.m` — change reaction equations. Ported as
`change_reaction_equations` ([manipulation/change.py](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/manipulation/change.py)).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| C1 | EFFICIENCY | raven-python 🔨 | 🔨 | **Edit the reaction in place** instead of RAVEN's copy-all-fields → `removeReactions` → `addRxns` → `permuteModel` round-trip. cobra mutates the same `Reaction` object, so every other field and the model order are preserved for free, with no O(n) re-sort. (Not a MATLAB back-port — the round-trip is inherent to the struct layout there.) |

## getIndexes

RAVEN `queries/getIndexes.m` — resolve a list of IDs / logical mask / index vector into positional
indexes (or a logical array) for `rxns` / `mets` / `genes` / `metNames` / `metcomps` (and GECKO
`ec.*` fields).

**Decision (raven-python): do NOT port the function.** cobra is object-oriented, so the central
index-resolver that RAVEN's struct-of-parallel-arrays design requires is largely unnecessary.
cobra's `DictList` already covers the use cases more idiomatically — `get_by_any` (mixed
id/object/index → objects), `get_by_id` (O(1)), `query` (name/substring/regex), `index` (position),
list comprehensions for filtering. Porting a 1-based-index resolver would be redundant and
un-Pythonic. **Only** the `name[comp]` composite resolver is kept (G7), as a small internal helper.

The improvement insights below still hold for **MATLAB RAVEN**, where the function remains — flagged
as upstream-only back-port candidates.

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| G1 | EFFICIENCY | MATLAB RAVEN only | 💡 | **Hash-based lookup instead of per-query linear scan.** RAVEN loops `find(strcmp(obj(i), searchIn))` per query → O(n·m). Build a `containers.Map` `{id: position}` once (O(n)), look up in O(1); `metcomps` likewise. (Moot for raven-python — cobra's `DictList` is already hashed.) |
| G5 | ERGONOMICS (bug) | MATLAB RAVEN only | 💡 | **Disambiguate the `[1 1 1]` mask-vs-index bug.** RAVEN's `if all(objects)` conflates a logical all-true mask with the index vector `[1 1 1]` (its own comment: "This gets weird if it's all 1"). Test `islogical(objects)` explicitly instead of `all(objects)`. (Moot for raven-python — input kind is decided by dtype.) |
| G7 | NEW | raven-python helper | 💡 | **Extract a reusable `name[comp]` parser/resolver.** The composite-id parsing buried in `getIndexes`'s `metcomps` branch is the one capability cobra lacks. Expose as a standalone `parse_name_comp` / `resolve_metabolite_by_name_comp`, reused by `addRxns`/`addTransport`/`mergeModels`. This is the only piece carried into raven-python. |

**Obsoleted by cobra (no action — these were earlier raven-python proposals now covered by `DictList`):**
predictable return type, return objects-not-positions, configurable missing-object policy across a
batch, and substring/regex matching — all already provided by `get_by_any` / `get_by_id` / `query`.

---

## standardizeGrRules

RAVEN `manipulation/standardizeGrRules.m` — normalize grRule syntax + flag rules not in simple
OR-of-AND-complex (DNF) form (`findPotentialErrors`).

**Decision (raven-python): port the lint half only.** cobra auto-normalizes a GPR on assignment
(`"(G1 AND G2)  OR  G3"` is stored as `"(G1 and G2) or G3"`), so the normalization half is
redundant. The non-DNF lint has no cobra equivalent and was ported as `find_non_dnf_grrules`/`is_dnf`
([utils/gpr.py](https://github.com/SysBioChalmers/raven-python/blob/develop/src/raven_python/utils/gpr.py)).

| # | Cat | Target | Status | Improvement |
|---|---|---|---|---|
| S1 | ERGONOMICS | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **Return structured lint results, don't just print.** RAVEN's `findPotentialErrors` only emits a `warning()` string; you can't act on it programmatically. raven-python returns a list of `GPRIssue(reaction_id, gpr, reason)`. **MATLAB back-port:** return the `indexes2check`/messages as a struct array (it already computes `indexes2check` — just surface it cleanly instead of only warning). |
| S2 | EFFICIENCY (robustness) | raven-python 🔨 + MATLAB RAVEN 💡 | 🔨 | **Detect non-DNF via the boolean AST, not substring search.** RAVEN scans for the `) and (`, `) and`, `and (` substrings, which is brittle (sensitive to spacing/bracketing and to gene IDs containing those characters). raven-python walks cobra's GPR AST (`is_dnf`: no OR beneath any AND), which is exact. **MATLAB back-port:** parse the rule to a tree rather than string-matching. |
