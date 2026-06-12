# MATLAB RAVEN — proposed back-ports

This is a consolidated list of fixes and improvements that the Python port surfaced
and that are worth carrying upstream into MATLAB RAVEN. Each item names a RAVEN
file/function, briefly diagnoses the current behaviour, and proposes the minimal
MATLAB-side patch. Items are sourced from [IMPROVEMENTS.md](improvements.md)
(the `MATLAB RAVEN 💡` rows) plus two new items from the section-F quality sweep
in [docs/known_issues.md](known_issues.md).

The matching Python-side implementation is linked for context — in most cases the
fix is identical in spirit; the patch shape just differs.

> Status legend used below: 💡 *proposed back-port* · 🐛 *bug* · ⚡ *efficiency* ·
> 🧹 *ergonomics / readability*.

---

## `reconstruction/homology/getModelFromHomology.m`

Implemented in raventoolbox as [`reconstruction.homology.get_model_from_homology`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/reconstruction/homology/homology.py).

* **H1** 🧹 — Split the overloaded `strictness` 1/2/3 into two orthogonal
  options: `bidirectional` (reciprocal hits) and `bestHitsOnly`. Reciprocal
  best hits = both true. Keep `strictness=` as a compat alias. Easier to read
  at call sites and easier to combine.
* **H2** 🐛 — Rewrite GPRs by walking a parsed boolean tree (the same AST
  RAVEN already builds for other checks), not via `regexprep` substring
  substitution. The current regex approach has partial-match hazards on gene
  IDs that contain `OLD_` or share prefixes, and needs a separate
  cleanup pass to fix degenerate `(OLD_… or …)` rules. Doing it on the tree
  is exact and removes the cleanup step.
* **H4** 🐛 — Default best-hit selection to **bitscore**, not E-value. Bitscore
  is database-size-independent (the RBH standard); the current E-value
  preference can flip the chosen template when the database grows. Keep
  E-value as an option.

## `reconstruction/kegg/getKEGGModelForOrganism.m` and the KEGG pipeline

Implemented across [`reconstruction.kegg`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/reconstruction/kegg/).

* **K1** 🐛 — In the KEGG flat-file parser, read each reaction's equation from
  its **own `EQUATION` field**, not by matching line *i* of `reaction.lst`
  against reaction *i*. The current line-matched approach silently corrupts
  the model if `reaction.lst` and `reaction` ever drift out of perfect order.
  The `EQUATION` field is already in the entry — just read it.
* **K2** 🧹 — Undefined-stoichiometry terms (`n C00001`, `(n+1) C00002`) keep
  the **real compound id** with coefficient 1 and *flag the reaction* as
  undefined-stoichiometry, instead of minting `"n C00001"` pseudo-metabolites
  that get renamed `undefined_N` later. Cleaner metabolite graph; the flag
  still drives the `keep*` filters.
* **K3** 🧹 — Replace the free-text reaction-quality strings appended to
  `rxnNotes` (`spontaneous`, `undefined stoichiometry`, `incomplete reaction`,
  `general reaction`) with a structured boolean **`rxn_flags`** table.
  Downstream filters join on a column instead of substring-matching notes.
* **K6** ⚡ — Per-KO multi-FASTA construction via a single streaming pass
  with an offset index (`seek`-based), replacing the Java `Hashtable`
  byte-scan in `constructMultiFasta`. Drops the 5M-element preallocation
  and the Java heap tuning.
* **K7** ⚡ — Concatenate all per-KO HMMs and `hmmpress` them into one
  pressed library, so the query path runs a single `hmmscan` against the
  database instead of thousands of per-KO `hmmsearch` invocations. Same
  result, orders-of-magnitude less I/O.
* **K12** ⚡ — Use fast MAFFT (`FFT-NS-2`, i.e. `--retree 2 --maxiterate 0`)
  for HMM training instead of `--auto`. `--auto` selects slow iterative
  refinement on medium/large KOs (observed ~2.5 min/KO, days for a domain
  on real KEGG); `FFT-NS-2` is seconds/KO and ample for profile-HMM
  building. Switch to `PartTree` when the alignment DP cost
  ≈ `n_seqs × mean_len²` exceeds the host's available memory budget (the
  Python implementation calibrated this as ≈4.2e-9 × DP-cost in GB and
  switches when DP-cost > 0.65 × (RAM − 2.5 GB) / 4.2e-9). Same protein
  set, no sequences dropped.
* **K14** ⚡ — Sort `organism_gene_ko` by `(organism, gene)` and store it
  gzipped (`organism_gene_ko.tsv.gz`). Sorting makes gene IDs within an organism
  adjacent (shared locus-tag/numeric prefixes), so they compress far better
  (~78 → 48 MB) and the order matches the by-organism query. Stored as gzip (not
  xz) so MATLAB reads it with built-in `gunzip`, since the same artefact is shared
  with raven-python; xz would be ~3× smaller but needs an external `unxz`.
* **K15** 🐛 — Recalibrate the HMM-query KO-assignment defaults: change
  `cutoff` from `1e-50` to `1e-30` and `min_score_ratio_g` from `0.8` to
  `0.9`. `min_score_ratio_ko` can stay at `0.3` but is empirically inert
  (identical output at 0.0/0.3/0.5 across the validation organisms).
  Cross-validated against the true KEGG gene→KO annotation of *S. cerevisiae*,
  *Cyanidioschyzon merolae*, *E. coli* K-12, and *Mycoplasma genitalium*:
  the current `1e-50` cutoff sits inside the *true-positive* tail
  (median true E ≈ 1e-100…1e-155; spurious ≈ 1e-8), silently dropping
  divergent real hits. At the proposed values, *M. genitalium* gene→KO recall
  rose from 0.84 → 0.94 (reaction recall 0.87 → 0.97) with no precision loss.
  Full numbers in [docs/kegg_hmm_cutoff_calibration.md](../studies/kegg_hmm_cutoff_calibration.md).

## `analysis/FSEOF.m`

Implemented in raventoolbox as [`analysis.fseof`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/analysis/fseof.py).

* **FS1** 🐛 — Replace the strict step-by-step monotonicity gate (a target is
  discarded if any single step's flux fails to exceed the previous) with a
  linear-regression filter: a reaction is a target if its flux track is
  reasonably correlated with the enforced product flux (|r| above a
  threshold). One noisy step from LP alternative optima no longer kills a
  valid target. Run `pFBA` per step to keep the scan stable.
* **FS2** 🧹 — Report **knockdown** and **knockout** targets as first-class
  outputs, not just amplification. The currently-flagged "amplify" set are
  reactions whose `|flux|` rises with the enforced product; the reactions
  driven *toward zero* (knock-down / knock-out candidates — arguably the most
  actionable for strain design) are silently dropped.
* **FS4** 🐛 — Replace the reported "Slope" — currently
  `abs(fseof.results(num, iterations) - fseof.results(num, 1)) / abs(targetMax - targetMax/iterations)`,
  i.e. an endpoint difference — with the **regression slope across all
  iterations**:

  ```matlab
  enforcedFlux = (1:iterations) * targetMax / iterations;   % record per-iter
  ...
  for num = 1:length(model.rxns)
      p = polyfit(enforcedFlux, fseof.results(num,:), 1);
      slope(num) = p(1);
  end
  ```

  An endpoint difference disagrees with the monotonicity filter on traces
  whose first and last values happen to be similar; the regression slope
  matches the selection criterion. The amplify/knockdown/knockout label
  (after FS2) should use the slope of `|flux|`, not endpoint comparison.

## `analysis/randomSampling.m`

Implemented in raventoolbox as [`analysis.random_sampling`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/analysis/sampling.py).

* **SAMP1** ⚡ — Compute `goodRxns` (loop-free, flux-carrying objective
  candidates) via a **single FVA pass**, not the current per-reaction `parfor`
  that solves a separate LP minimising and maximising every reaction. FVA
  computes the same min/max ranges, optionally `loopless` (`cycleFreeFlux`),
  in far less code.
* **SAMP2** 🐛 — Add a `seed` argument for reproducibility, and make
  `nObjectives` a parameter. The current code hard-codes 2 objective
  reactions per sample (the docstring even claims 3 — a transcription bug
  worth fixing alongside).
* **SAMP4** 🐛 — Run `goodRxns` detection **before** `replaceBoundsWithInf`
  applies, not after. FVA errors as "unbounded" on infinite bounds, and
  running goodRxns on the inf-replaced model conflates "loop reaction" with
  "unbounded objective". Scope `replaceBoundsWithInf` to sampling only, after
  the goodRxns set is found.

## `analysis/reporterMetabolites.m`

Implemented in raventoolbox as [`analysis.reporter_metabolites`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/analysis/reporter.py).

* **RM1** ⚡🐛 — Replace the per-neighbour-count Monte Carlo background
  correction (currently 100 000 random sets drawn *for each distinct
  neighbour count*) with the exact closed-form expression. RAVEN samples with
  replacement from the scored-gene Z pool, so a random aggregate Z = Σz/√n
  has mean `√n · μ` and standard deviation `σ` provably — the corrected
  metabolite Z is just `(metZ − √n·μ) / σ`. This removes both the slow
  sampling step and its run-to-run randomness (results become deterministic).

## `INIT/runINIT.m`

Implemented in raventoolbox as [`init.run_init`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/init.py).

* **I4** 🐛 — Drop the hard-coded big-M (1000) in the MILP and use **each
  reaction's own upper bound** instead: `v ≤ ub · x`. Expose `eps` and
  `prod_weight` as parameters. The current 1/1000/0.1/0.5 quadruple only
  fits ±1000-bounded models with O(1) scores — flagged in RAVEN as
  scale-dependent and tunable, but no API knobs are exposed.

## `INIT/ftINIT.m` and the ftINIT pipeline

Implemented across [`init.ftinit`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/ftinit.py), [`init.taskfill`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/taskfill.py), [`init.merge`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/merge.py), [`init.prep`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/prep.py).

* **FT3** 🐛 — Same big-M issue as `runINIT` — use each reaction's own bound
  as big-M instead of the fixed 100/1000. Expose `force_on` and `force_on_ess`
  as parameters (the calibrated values for genome-scale Human-GEM are
  documented in [docs/init_param_calibration.md](../studies/init_param_calibration.md)).
* **FT9** 🐛 — Per-reaction essential forcing must be **clamped to the
  reaction's bound**. The staged pipeline fixes previous-step reactions as
  essential at `min(0.99 · |prev flux|, 0.1)`; if a low-capacity reaction
  cannot carry that magnitude, the current code produces `lb > ub` and
  errors. Clamp the forced magnitude to `rxn.ub` (in the positive direction)
  / `rxn.lb` (in the negative direction) so a low-capacity essential is
  forced *to* its capacity instead.
* **FT8** 🧹 — Rewrite `removeLowScoreGenes` on the parsed boolean tree
  instead of regex subset-extraction with `#n#` placeholders. Same logic
  (isozyme OR-branches dropped while keeping at least one; complex
  AND-subunits kept intact; nested groups handled), but exact and tractable
  to test. Make the all-negative tie-break **deterministic** (sort + take
  the least-negative) instead of random.

## `manipulation/addRxns.m`

Implemented in raventoolbox as [`manipulation.add_reactions_from_equations`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/add.py).

* **A1** 🧹 — Accept a string keyword (`metsBy = 'id'` / `'name'`) instead of
  the opaque `eqnType = 1 / 2 / 3` integer. Call-sites become self-documenting.

## `queries/checkModelStruct.m` (or its curation subset)

Implemented in raventoolbox as [`utils.check_model`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/utils/validate.py).

* **V2** 🧹 — Return a struct array of issues (one per finding, with
  `category` / `target` / `message` fields) instead of printing warnings or
  throwing on the first problem. Programmatically filterable.

## `queries/getElementalBalance.m`

Implemented in raventoolbox as [`utils.get_elemental_balance`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/utils/balance.py).

* **B2** 🐛 — A reaction with no metabolites (an empty `S(:,j)`) currently
  falls through to `balanceStatus = 1` ("balanced") because both pre-loops
  leave `balanceStatus` at NaN and the final `isnan→1` step maps NaN to
  balanced. Fix:

  ```matlab
  % Right before the final `balanceStatus(isnan(...)) = 1;` line:
  emptyRxns = full(sum(model.S ~= 0, 1)) == 0;
  balanceStatus(emptyRxns) = min(-1, balanceStatus(emptyRxns));
  ```

  Empty reactions then report `-1` ("missing information") rather than
  vacuously balanced.

## `manipulation/findPotentialErrors.m` (GPR linting)

Implemented in raventoolbox as [`utils.find_non_dnf_grrules`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/utils/gpr.py).

* **S1** 🧹 — Return the `indexes2check` array (which the function already
  computes) plus per-reaction reason strings as a struct array, instead of
  only emitting a `warning()`. Callers cannot act on the warning text.
* **S2** 🐛 — Detect non-DNF GPRs via the parsed boolean tree, not substring
  search for `) and (`, `) and`, `and (`. The substring approach is sensitive
  to spacing and bracketing and to gene IDs containing those characters; the
  tree walk (no OR beneath any AND) is exact.

## `queries/getIndexes.m`

* **G1** ⚡ — Build a `containers.Map` `{id: position}` once per call and use
  hashed lookup (O(1)) instead of looping `find(strcmp(obj(i), searchIn))`
  per query (O(n · m)). Same for the `metcomps` branch.
* **G5** 🐛 — `if all(objects)` conflates a logical all-true mask with the
  index vector `[1 1 1]` (the function's own comment notes "this gets weird
  if it's all 1"). Test `islogical(objects)` explicitly instead.

## Cross-cutting

* **Y4** 🧹 — A first-class home for `metDeltaG` / `rxnConfidenceScores` (in
  RAVEN) and `deltaG` / `confidence_score` (in cobra/raventoolbox `.notes`).
  Neither cobra nor SBML core has a standard slot yet; if one emerges (SBML
  fbc/groups, or a cobra attribute) both projects should migrate. For now
  this is a watching brief.

## R-MetaCyc — removal candidate

Already flagged in `IMPROVEMENTS.md` as a removal target on both sides — the
MetaCyc reconstruction path is dropped in raventoolbox and proposed for removal
from MATLAB RAVEN (`external/metacyc/*`). Pasting here for visibility because
the actual removal still needs to land upstream. See
[IMPROVEMENTS.md § R-MetaCyc](improvements.md) for the full rationale.
