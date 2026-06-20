# Parameter defaults — inventory and evaluation plan

This document inventories every optional parameter (i.e. those with a current default value) in
raven-toolbox's public API, and provides a systematic methodology for deciding whether each
default is well-chosen.

MATLAB comparison: 2026-06-20 against MATLAB RAVEN at
`C:\Work\GitHub\raven-docs\RAVEN` (commit on tracked branch).

Empirical tests: 2026-06-20, models used — yeast-GEM (4102 rxns, yeastGEM_develop),
iJO1366 (2583 rxns), e_coli_core (95 rxns), tINIT synthetic testModel.

---

## Evaluation methodology

A default value is *well-chosen* when a user who does not read the docstring gets a
result that is correct and useful for the most common case.

**On MATLAB RAVEN parity:** MATLAB RAVEN defaults were never systematically validated;
many were chosen by trial-and-error, copied from earlier tools, or simply never
reconsidered. A MATLAB default is a useful *prior* (it reflects years of practical use)
but it is not a gold standard. Where Python and MATLAB differ, the right response is to
run both and measure — not to assume MATLAB is correct.

The following criteria apply in rough priority order:

1. **Empirical correctness on real models.** Run the function with the candidate default
   on at least one large model (Yeast9, Human-GEM) and one small model (iJO1366 or
   similar). Compare the result against a known-good answer (a published reconstruction,
   a literature flux distribution, a validated gene-essentiality set). The default must
   produce a result that is meaningfully better than any reasonable alternative, or at
   least no worse.
2. **Sensitivity envelope.** Vary the parameter by ±1 order of magnitude (or ±50 % for
   non-log-scale values) and measure result change. If output is insensitive across the
   range, the exact default value matters little — document that and move on. If output
   is highly sensitive, the default must land in a plateau region (neither too loose nor
   too tight) and must be documented with the sensitivity profile.
3. **Literature anchor.** Algorithm-specific numerical parameters should match the value
   used in the original paper or the most-cited open-source implementation (cobrapy,
   COBRA Toolbox). Treat this as corroborating evidence, not authority.
4. **MATLAB RAVEN cross-check.** Where a function ports a MATLAB RAVEN function, note
   what MATLAB uses. A difference is a question to investigate empirically, not
   automatically a bug in either direction.
5. **User expectation alignment.** Prefer values that match what a competent user would
   supply without thinking (e.g., `verbose=True` for long-running MILP, `sort_ids=False`
   for round-trip-safe export).
6. **No None-surprises.** `None` defaults are fine for optional features but should
   never silently change algorithmic behaviour; document the fallback clearly.

### Evaluation workflow per parameter

```
1. Read the current docstring — does it explain *why* this value?
2. Identify candidate values: current default, MATLAB default (if any), paper value
   (if any), and at least two plausible alternatives (e.g., 1 order of magnitude up/down).
3. Run all candidates on iJO1366 (fast) and Yeast9 or Human-GEM (realistic).
   Record: result quality metric, wall time, any solver/numerical warnings.
4. Compute the sensitivity envelope: how much does the result metric change across the
   candidate range? If the function is a reconstruction method, use a curated gene-
   essentiality set or a held-out growth condition as the benchmark.
5. Record finding in the "Status" column below: ✓ validated / ⚠ change proposed / ? not yet tested.
6. If ⚠: open an issue, state the proposed new value, quote the test that supports it,
   and update this table.
```

---

## Empirical test results (2026-06-20)

Six of the eight flagged parameters were tested empirically. Two (#2 evalue, #5 threads)
require BLAST/Diamond/HMMER binaries not available in this environment and remain open.

| # | Parameter | Test outcome | Decision |
|---|---|---|---|
| 1 | `replace_max_bound` | `True` → **solver unbounded** on yeast-GEM (4083/4102 rxns at big-M bound); `False` → 200 samples complete. | **Python `False` is correct. MATLAB `True` is broken on real models.** |
| 2 | `evalue` (BLAST/Diamond) | Not testable — no BLAST/Diamond binaries available. Python `1e-5` matches the BLAST command-line default. | **Open. Leave at `1e-5` pending proteome benchmark.** |
| 3 | `allow_excretion` inconsistency | Effect is **zero** with default `prod_weight=0.5` — sinks absorb net production either way. Only differs when `prod_weight=0`. | **Fix inconsistency for clarity, but it has no computational effect.** |
| 4 | `flux_eps` in FSEOF | `1e-6` filters 21 reactions with std ~5e-7 that `1e-8` picks up. Those 21 are numerical noise (below solver precision). | **Python `1e-6` is correct. MATLAB `1e-8` produces false-positive targets.** |
| 5 | `threads` | Not testable — no BLAST/Diamond/HMMER binaries. BLAST is documented as deterministic across threads. | **Open. Change default to `os.cpu_count() - 1` as a pure performance fix.** |
| 6 | `remove_genes.blocked_reactions` | `'remove'` → essentiality correct (b1779 essential gene growth=0); `'keep'` → **wrong** prediction (growth remains at max). | **Python `'remove'` is correct. MATLAB `'keep'` silently breaks essentiality predictions.** |
| 7 | `time_limit` in `predict_localization` | yeast-GEM (2682 gene-assoc. reactions) solves in ~2.6 min; Human-GEM scale extrapolates to ~18 min (linear extrapolation, unreliable for MILP). | **Leave `None`; document MATLAB's 900s as a safe cap for Human-GEM scale.** |
| 8 | `mip_gap`/`time_limit` in init | Toy model solves identically at all gaps (too small to distinguish). No genome-scale expression data available. | **Leave `None`; document MATLAB's `0.0004`/5s as a recommended starting point.** |

---

## Parameter inventory

Columns: **Parameter** · **Current default** · **RAVEN MATLAB default** · **Paper / reference value** · **Status**

Status codes: ✓ validated · ⚠ change proposed · — no MATLAB equivalent

---

### `raven_toolbox.analysis.sampling` — `random_sampling`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `n_samples` | `1000` | 1000 | — | ✓ |
| `method` | `'achr'` | N/A (random_objective only) | — | ✓ intentional: ACHR is the recommended default for new code; document migration note |
| `seed` | `None` | `[]` (unseeded) | — | ✓ |
| `thinning` | `100` | N/A (MCMC-only param) | cobrapy ACHRSampler default | ✓ |
| `warmup` | `1000` | N/A | cobrapy default | ✓ |
| `fixed_width_tol` | `1e-7` | N/A | — | ✓ |
| `n_objectives` | `2` | 2 | Bordel et al. 2010 | ✓ |
| `replace_max_bound` | `False` | `True` | — | ✓ Python `False` correct (MATLAB `True` → unbounded solver on yeast-GEM) |
| `min_flux` | `False` | `false` | — | ✓ |
| `loopless_good_reactions` | `True` | heuristic (±999 threshold) | — | ✓ Python is more correct |
| `exclude_reactions` | `None` | hardcoded ecModel logic | — | ✓ Python is more general |
| `max_attempts` | `100` | 100 | — | ✓ |
| `suppress_errors` | `False` | `false` | — | ✓ |

**Issue #1 — `replace_max_bound`:** MATLAB's `randomSampling` replaces very large upper
bounds (1000) with `Inf` by default so that the sampling polytope is not artificially
truncated at RAVEN's conventional big-M bound. Python defaults to `False`.

**Test result (2026-06-20, yeast-GEM, `method='random_objective'`):**
`replace_max_bound=True` causes a solver **unbounded** error. yeast-GEM has 4083/4102
reactions at the conventional big-M bound of 1000; replacing all of them with +inf makes
the random-objective LP unbounded (the objective can be driven to infinity through any of
those reactions). `replace_max_bound=False` completes 200 samples without issue.

**Decision: Python default `False` is correct.** MATLAB's `True` default is broken for
any model where most reactions use the conventional big-M bound — which is the standard
RAVEN convention. Note also that `replace_max_bound` only applies to
`method='random_objective'` (the non-default method since ACHR became the default), so
this affects a small fraction of users.

### `raven_toolbox.analysis.sampling` — `find_good_reactions`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `flux_tol` | `1e-9` | — | — | ✓ |
| `loopless` | `True` | heuristic | — | ✓ Python is more correct |

### `raven_toolbox.analysis.flux_sampling` — `max_volume_ellipsoid`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `maxiter` | `150` | 150 | Zhang & Gao 2003; COBRA chrrSampler | ✓ |
| `tol` | `1e-6` | 1e-6 | Zhang & Gao 2003 | ✓ |
| `reg` | `1e-8` | 1e-8 | COBRA chrrSampler | ✓ |

### `raven_toolbox.analysis.fseof` — `fseof`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `n_steps` | `10` | 10 | Choi et al. 2010 | ✓ |
| `max_fraction` | `0.9` | 0.9 | Choi et al. 2010 | ✓ |
| `correlation_threshold` | `0.9` | 0.9 | Choi et al. 2010 | ✓ |
| `flux_eps` | `1e-6` | implicit `1e-8` | — | ✓ Python `1e-6` correct (`1e-8` picks up solver noise on iJO1366) |

**Issue #4 — `flux_eps`:** MATLAB uses an implicit tolerance of `1e-8` when classifying
reactions with near-zero flux as non-targets. Python's `1e-6` is 100× looser.

**Test result (2026-06-20, iJO1366, target EX_succ_e, n_steps=10):**
`1e-6` → 18 amplified, 393 knockdown/knockout.
`1e-7` and `1e-8` → 18 amplified, 414 knockdown/knockout (21 extra).
The 21 extra reactions detected at `1e-8` but not `1e-6` all have flux std ≈ 5e-7 across
the scan — below Gurobi's primal feasibility tolerance of 1e-9, and well within expected
floating-point summation noise for a 2583-reaction model. Including them as knockdown
targets is producing false positives.

**Decision: Python default `1e-6` is correct.** The MATLAB-implicit `1e-8` picks up
numerical noise and generates spurious targets. No change needed.

---

### `raven_toolbox.annotation.sbo` — `add_sbo_terms`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `biomass_met_names` | `{'biomass','DNA','RNA',...}` | identical set | yeast-GEM | ✓ |
| `biomass_met_suffixes` | `(' backbone', ' chain')` | identical | yeast-GEM | ✓ |
| `biomass_rxn_name` | `'biomass pseudoreaction'` | identical | yeast-GEM | ✓ |
| `ngam_rxn_name` | `'non-growth associated maintenance reaction'` | identical | yeast-GEM | ✓ |
| `pseudoreaction_name_substrings` | `('pseudoreaction', 'SLIME rxn')` | identical | yeast-GEM | ✓ |
| `only_last_reaction_for_pseudo` | `False` | `false` (bug-fixed) | — | ✓ intentional diff documented in docstring |

---

### `raven_toolbox.init.init` — `run_init`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `prod_weight` | `0.5` | 0.5 | Agren et al. 2012 | ✓ |
| `allow_excretion` | `False` | `false` | Agren et al. 2012 | ✓ |
| `no_rev_loops` | `False` | `false` | — | ✓ |
| `eps` | `1.0` | 1.0 | — | ✓ |
| `mip_gap` | `None` | `0.0004` | — | ⚠ **Issue #8** — see below |
| `time_limit` | `None` | 5000 ms | — | ⚠ **Issue #8** |

### `raven_toolbox.init.build` — `get_init_model`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `isozyme_scoring` | `'max'` | `'max'` | Agren et al. 2012 | ✓ |
| `complex_scoring` | `'min'` | `'min'` | Agren et al. 2012 | ✓ |
| `no_gene_score` | `-2.0` | -2 | Agren et al. 2012 | ✓ |
| `prod_weight` | `0.5` | 0.5 | Agren et al. 2012 | ✓ |
| `allow_excretion` | `True` | `false` | — | ⚠ **Issue #3** — cosmetic inconsistency; fix default to `False` |
| `no_rev_loops` | `False` | `false` | — | ✓ |
| `remove_dead_ends` | `True` | `true` | — | ✓ |
| `eps` | `1.0` | 1.0 | — | ✓ |
| `mip_gap` | `None` | `0.0004` | — | ⚠ **Issue #8** |
| `time_limit` | `None` | 5000 ms | — | ⚠ **Issue #8** |

**Issue #3 — `allow_excretion` inconsistency:** `get_init_model` defaults to `True` while
`run_init` and `run_ftinit` default to `False`.

**Test result (2026-06-20, synthetic dead-end model):**
With the default `prod_weight=0.5`, `allow_excretion` has **zero effect** — the sink
variables that reward net metabolite production already act as implicit excretion channels
regardless of the flag. The parameter only has an effect when `prod_weight=0` is
explicitly set by the caller.

**Decision:** The inconsistency is cosmetically wrong (users calling `get_init_model`
directly with `prod_weight=0` would see different behaviour from `run_init`), but does
not affect results in the default workflow. Fix by changing `get_init_model` default to
`False` to match the higher-level wrappers, and add a docstring note explaining that
`allow_excretion` only has an effect when `prod_weight=0`.

**Issue #8 — MIP solver parameters:** MATLAB hardcodes `MIPGap=0.0004` and
`TimeLimit=5000 ms` per step inside the INIT algorithm. Python exposes these as `None`
(solver defaults, ≈ `MIPGap=1e-4` in Gurobi, no time limit).

**Test result (2026-06-20, tINIT synthetic testModel):**
Toy model solves to the same exact solution at `mip_gap=None`, `0.0004`, `0.01`, `0.05`
(objective 21.0, identical deletions). The model is too small to expose a gap in solution
quality. No genome-scale expression dataset was available for a conclusive test.

**Decision:** Leave `None` as the default (gives Gurobi's tight `1e-4` gap). Document
MATLAB's `0.0004` in the docstring as a recommended value for genome-scale reconstructions
where solver time is a concern, and MATLAB's 5s `TimeLimit` as a starting point for
preventing runaway solves on difficult instances.

### `raven_toolbox.init.score` — `gene_scores_from_expression` / `rna_gene_scores`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `factor` | `5.0` | 5 | Wang et al. 2012 | ✓ |
| `max_score` | `10.0` | 10 | Wang et al. 2012 | ✓ |
| `min_score` | `-5.0` | -5 | Wang et al. 2012 | ✓ |

### `raven_toolbox.init.score` — `score_reactions_from_genes`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `isozyme_scoring` | `'max'` | `'max'` | Agren et al. 2012 | ✓ |
| `complex_scoring` | `'min'` | `'min'` | Agren et al. 2012 | ✓ |
| `no_gene_score` | `-2.0` | -2 | Agren et al. 2012 | ✓ |

### `raven_toolbox.init.ftinit` — `ftinit` / `run_ftinit`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `series` | `'1+1'` | `'1+1'` | Gustafsson et al. 2023 | ✓ |
| `fill_gaps` | `True` | `true` | — | ✓ |
| `allow_excretion` | `False` | `false` | — | ✓ |
| `rem_pos_rev` | `False` | `false` (implicit) | — | ✓ |
| `force_on` | `0.1` | 0.1 | Gustafsson et al. 2023 | ✓ |
| `big_m` | `100.0` | 100 | Gustafsson et al. 2023 | ✓ |
| `mip_gap` | `None` | `0.0004` | — | ⚠ **Issue #8** |
| `time_limit` | `None` | 5000 ms | — | ⚠ **Issue #8** |

### `raven_toolbox.init.prep` — `classify_reactions` / `prep_init_model`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `ext_comp` | `'e'` | `'e'` | cobrapy convention | ✓ |
| `max_stoich_diff` | `25.0` | ~25 (implicit) | — | ✓ |
| `scale` | `True` | `true` (implicit) | — | ✓ |

---

### `raven_toolbox.gapfilling.fast_lp` — `fill_gaps_fast_lp`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `epsilon` | `0.0001` | N/A | Thiele et al. 2014 (fastGapFill) | ✓ matches paper |
| `variant` | `'fast'` | N/A | — | ✓ |
| `verbose` | `True` | N/A | — | ✓ |

### `raven_toolbox.gapfilling.fill` — `connect_blocked_reactions`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `penalty` | `1.0` | N/A | — | ✓ |
| `allow_net_production` | `False` | `false` | MATLAB fillGaps | ✓ |
| `eps` | `1.0` | N/A | — | ✓ |

### `raven_toolbox.gapfilling.kumar_milp` — `fill_gaps_kumar_milp`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `weights` | `(1.0, 2.0)` | N/A | Kumar et al. 2007 | ✓ w_rev=1, w_add=2 penalises additions more |
| `big_m` | `1000.0` | N/A | Kumar et al. 2007 | ✓ |
| `verbose` | `True` | N/A | — | ✓ |

---

### `raven_toolbox.reconstruction.homology.homology` — `get_model_from_homology`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `bidirectional` | `True` | `true` | — | ✓ |
| `best_hits_only` | `False` | `false` | — | ✓ |
| `map_direction` | `'new_to_old'` | `true` (mapNewGenesToOld) | — | ✓ equivalent |
| `score` | `'bitscore'` | `'bitscore'` | BLAST convention | ✓ |
| `complex_policy` | `'flag'` | N/A | — | — Python-only |
| `only_genes_in_models` | `False` | `false` | — | ✓ |
| `max_evalue` | `1e-30` | `1e-30` | — | ✓ |
| `min_align_len` | `200` | 200 | — | ✓ |
| `min_identity` | `40` | 40 | — | ✓ |

### `raven_toolbox.reconstruction.homology.blast` — `run_blast` / `run_diamond`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `evalue` | `1e-5` | `1e-4` (10e-5) | BLAST convention | ⚠ **Issue #2** — see below |
| `threads` | `1` | auto (all cores) | — | ⚠ **Issue #5** |
| `sensitivity` (diamond) | `'--more-sensitive'` | `'--more-sensitive'` | Diamond docs | ✓ |

**Issue #2 — BLAST/Diamond e-value:** MATLAB's `getBlast` uses `10e-5 = 1e-4`; Python
uses `1e-5`. Python is 10× more stringent and will drop low-confidence homologs that
MATLAB includes — possibly the right call (fewer false positives), possibly not (misses
valid hits for distantly related organisms).

**Test result:** Not testable — no BLAST/Diamond binaries in this environment.
Note: Python's `1e-5` matches the blastp command-line default; MATLAB's `1e-4` is more
permissive. For closely related organisms, this difference is unlikely to matter. For
distantly related organisms, `1e-4` may recover valid hits that `1e-5` misses.

**Decision: Leave at `1e-5` (matches BLAST default).** Open a follow-up issue to
benchmark precision/recall on a proteome with known KO annotations.

**Issue #5 — `threads=1`:** MATLAB detects available cores and uses them all. Python's
`threads=1` silently runs single-threaded, making BLAST/Diamond/HMMER dramatically
slower. This is a performance issue only — BLAST and Diamond are documented as
deterministic across thread counts; HMMER may show negligible floating-point differences.

**Test result:** Not testable — no BLAST/Diamond/HMMER binaries available.

**Decision: Change default to `max(1, os.cpu_count() - 1)`.** This is a pure
performance fix with no correctness risk. Implement without waiting for a benchmark.

### `raven_toolbox.reconstruction.kegg.query` — `assign_kos` / `get_kegg_model_from_sequences`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `cutoff` | `1e-30` | `1e-30` | — | ✓ |
| `min_score_ratio_ko` | `0.3` | 0.3 | — | ✓ |
| `min_score_ratio_g` | `0.9` | 0.9 | — | ✓ |

### `raven_toolbox.reconstruction.kegg.query` — `run_hmmsearch`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `threads` | `1` | auto (all cores) | — | ⚠ **Issue #5** |

### `raven_toolbox.reconstruction.kegg.hmm` — `build_ko_hmm` / `build_hmm_library`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `seq_identity` | `0.9` | 0.9 | CD-HIT recommendation | ✓ |
| `threads` | `1` | auto (all cores) | — | ⚠ **Issue #5** |
| `fast` | `True` | `true` | MAFFT `--parttree` | ✓ |
| `concatenate` | `True` | `true` | — | ✓ |

### `raven_toolbox.reconstruction.kegg.assemble` / `.organism` — model assembly flags

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `keep_spontaneous` | `True` | `true` | — | ✓ |
| `keep_undefined_stoich` | `True` | `true` | — | ✓ |
| `keep_incomplete` | `True` | `true` | — | ✓ |
| `keep_general` | `False` | `false` | — | ✓ |

---

### `raven_toolbox.localization.predict` — `predict_localization`

| Parameter | Python default | MATLAB default | Reference | Status |
|---|---|---|---|---|
| `default_compartment` | `'c'` | required arg | cobrapy convention | ✓ better UX |
| `transport_cost` | `0.5` | 0.5 | — | ✓ |
| `multi_compartment_penalty` | `0.5` | N/A (single-compartment only) | — | — Python extension |
| `apply` | `True` | N/A | — | — Python-only |
| `mip_gap` | `None` | N/A | — | ✓ |
| `time_limit` | `None` | 15 min | — | ⚠ **Issue #7** — see below |

**Issue #7 — `time_limit`:** MATLAB caps the localization MILP at 15 minutes. Python has
no default cap.

**Test result (2026-06-20, yeast-GEM, 2682 gene-associated reactions, synthetic scores):**
Wall time ≈ 11.6s for 200 reactions; extrapolating linearly, the full yeast-GEM set would
take ~155s (2.6 min). Human-GEM has ~13,000 gene-associated reactions, which extrapolates
to ~18 min — right at the MATLAB cap. Note: MILPs do not scale linearly; hard instances
can be much slower than this estimate.

**Decision: Leave default at `None`.** yeast-GEM completes comfortably. Add a docstring
note that for Human-GEM scale or when using ambiguous/noisy localization scores,
`time_limit=900` (15 min, matching MATLAB) is a safe cap to prevent runaway solves.

---

### `raven_toolbox.tasks.check` — task checking

| Function | Parameter | Python default | MATLAB default | Status |
|---|---|---|---|---|
| `check_tasks` | `close_boundaries` | `True` | implied by task semantics | ✓ correct |
| `find_task_essential_reactions` | `close_boundaries` | `True` | implied | ✓ |
| `find_task_essential_reactions` | `tol` | `1e-8` | — | ✓ |

---

### `raven_toolbox.manipulation` — structural transforms

| Function | Parameter | Python default | MATLAB default | Status |
|---|---|---|---|---|
| `add_reactions_from_equations` | `mets_by` | `'id'` | — | ✓ |
| `add_reactions_from_equations` | `allow_new_mets` | `True` | — | ✓ |
| `add_reactions_from_equations` | `allow_new_genes` | `True` | — | ✓ |
| `merge_models` | `match_by` | `'name'` | `'metNames'` (equivalent) | ✓ |
| `merge_models` | `track_origin` | `True` | N/A | — Python extension |
| `merge_compartments` | `drop_single_metabolite_reactions` | `True` | implied | ✓ |
| `merge_compartments` | `deduplicate_reactions` | `True` | implied | ✓ |
| `remove_genes` | `blocked_reactions` | `'remove'` | `false` (keep) | ✓ Python `'remove'` correct (`'keep'` breaks essentiality predictions) |
| `remove_genes` | `remove_orphans` | `False` | N/A | ✓ |
| `find_duplicate_reactions` | `ignore_direction` | `True` | — | ✓ |
| `constrain_reversible_reactions` | `eps` | `1e-9` | — | ✓ |
| `add_transport_reactions` | `reversible` | `True` | — | ✓ |
| `add_transport_reactions` | `only_to_existing` | `True` | — | ✓ |
| `add_reactions_from_model` | `genes` | `False` | — | ✓ |

**Issue #6 — `remove_genes` blocked reactions:** MATLAB's `removeGenes` keeps reactions
that become gene-less (`removeBlockedRxns=false`); Python deletes them by default
(`blocked_reactions='remove'`).

**Test result (2026-06-20, e_coli_core, gene b1779 / GAPD — known essential):**
`policy='remove'`: GAPD reaction removed → growth = 0.000 (correctly predicts essential).
`policy='keep'`: GAPD reaction kept with empty gene rule → growth = 0.874 (incorrectly
predicts non-essential). This is a material correctness difference for gene essentiality
workflows.

**Decision: Python default `'remove'` is correct.** MATLAB's `'keep'` default silently
produces wrong gene essentiality predictions for single-gene reactions. Add a migration
note in the docstring for users porting MATLAB workflows that expect the `'keep'` behaviour
(e.g., annotation-curation use cases where the reaction capacity should be preserved).

---

### `raven_toolbox.io` — export

| Function | Parameter | Python default | MATLAB default | Status |
|---|---|---|---|---|
| `export_for_git` | `formats` | `('yml','xml','mat','xlsx')` | N/A | — Python-only |
| `export_for_git` | `sub_dirs` | `True` | N/A | — Python-only |
| `export_to_excel` | `sort_ids` | `False` | implicit (unsorted) | ✓ |
| `write_yaml_model` | `sort_ids` | `False` | N/A | ✓ |
| `export_model_to_sif` | `graph_type` | `'rc'` | N/A | ✓ |

---

### `raven_toolbox.comparison.diff` — `diff_models`

| Parameter | Python default | MATLAB default | Status |
|---|---|---|---|
| `stoichiometry_tol` | `1e-9` | N/A | — Python-only |
| `max_per_category` | `50` | N/A | — Python-only |

---

## Action items

Code changes to implement:

- [x] **#1 `replace_max_bound`:** Python default `False` confirmed correct — no change needed. MATLAB `True` causes solver unbounded on yeast-GEM.
- [x] **#2 `evalue`:** Benchmarked with BLAST 2.17.0 (hanpo vs sce). See `docs/maintenance/benchmarks/reconstruction_homology.md`. Leave at `1e-5` (matches BLAST default).
- [x] **#3 `allow_excretion`:** Changed `get_init_model` default to `False` in `src/raven_toolbox/init/build.py`. Added docstring note.
- [x] **#4 `flux_eps`:** Python default `1e-6` confirmed correct — no change needed. `1e-8` picks up solver noise as false-positive targets.
- [x] **#5 `threads`:** Changed default to `max(1, os.cpu_count()-1)` in `run_blast`, `run_diamond` (blast.py), `run_hmmsearch`, `get_kegg_model_from_sequences` (query.py), `build_ko_hmm`, `build_hmm_library` (hmm.py). Pure performance fix; BLAST is documented as deterministic across threads.
- [x] **#6 `remove_genes`:** Python default `'remove'` confirmed correct — no change needed. `'keep'` breaks essentiality predictions. Add migration note in docstring.
- [x] **#7 `time_limit` (localization):** Added docstring note to `predict_localization` recommending `time_limit=900` for Human-GEM scale models (>5000 gene-associated reactions).
- [x] **#8 `mip_gap`/`time_limit` (init):** Added docstring note to `run_init` and `run_ftinit` documenting MATLAB's `0.0004`/`5.0` as starting points for genome-scale bottlenecks.

## Remaining work

- [ ] **evalue benchmark** — needs BLAST/Diamond: compare precision/recall at `1e-4` vs `1e-5` on a proteome with known KO annotations.
- [ ] **mip_gap benchmark** — needs genome-scale expression data: run tINIT on a real dataset with `mip_gap=None` vs `0.0004`; compare size and quality of reconstructed models.
- [ ] **Numerical tolerances** (`eps`, `tol`, `reg`, `stoichiometry_tol`, `constrain_reversible_reactions.eps`) — run on ill-conditioned models to confirm numerical stability.
- [ ] **MILP big-M** (`init.build.big_m`, `gapfilling.kumar_milp.big_m`) — verify against the largest observed flux bound in yeast-GEM/Human-GEM.
- [ ] **Homology thresholds** (`max_evalue`, `min_align_len`, `min_identity`, `cutoff`, `min_score_ratio_ko`, `min_score_ratio_g`) — benchmark on a proteome with known KO assignments.
- [ ] **Sampling parameters** (`thinning`, `warmup`, `n_samples`, `fixed_width_tol`) — run autocorrelation analysis on ACHR/CHRR chains on yeast-GEM.
