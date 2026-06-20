# Parameter defaults — inventory and evaluation plan

This document inventories every optional parameter (i.e. those with a current default value) in
raven-toolbox's public API, and provides a systematic methodology for deciding whether each
default is well-chosen.

Last evaluated: 2026-06-20 against MATLAB RAVEN at
`C:\Work\GitHub\raven-docs\RAVEN` (commit on tracked branch).

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

## Questions requiring empirical testing

The 2026-06-20 comparison pass flagged the following parameters as needing real-world
tests before a decision can be made. Where Python and MATLAB differ, this does **not**
mean MATLAB is correct — both should be tested. Each is expanded in the relevant section
below.

| # | Priority | Parameter | Question |
|---|---|---|---|
| 1 | **High** | `random_sampling.replace_max_bound` | Python `False`, MATLAB `True`. Does replacing big-M bounds with ±inf produce meaningfully different or more representative samples on ecModels? Test on ecYeast9 with ACHR/CHRR. |
| 2 | **High** | `run_blast / run_diamond.evalue` | Python `1e-5`, MATLAB `1e-4`. Does the stricter Python cutoff discard valid hits? Run on a proteome with known KO assignments, count true/false positives at both thresholds. |
| 3 | **Medium** | `get_init_model.allow_excretion` | `True` here, `False` in `run_init` / `run_ftinit`. Internal inconsistency — decide on one value backed by reconstruction quality on a test proteome. |
| 4 | **Medium** | `fseof.flux_eps` | Python `1e-6`, MATLAB implicit `1e-8`. Does the looser cutoff produce spurious low-flux targets? Run FSEOF on iJO1366 and compare the target reaction lists. |
| 5 | **Medium** | `run_blast / run_diamond / run_hmmsearch.threads` | Python `1`, MATLAB auto-detects. This is a performance issue, not a correctness one — fix is to default to available cores, but first establish that results are reproducible across thread counts. |
| 6 | **Medium** | `remove_genes.blocked_reactions` | Python removes, MATLAB keeps. Which default produces more useful models? Test by removing a known non-essential gene from iJO1366 and checking network integrity. |
| 7 | **Low** | `predict_localization.time_limit` | MATLAB caps at 15 min. Does omitting the cap cause unbounded solves on large models? Profile on Human-GEM. |
| 8 | **Low** | `run_init / run_ftinit` MIP gap | MATLAB `0.0004`, Python `None` (solver default ≈ `1e-4`). Does the MATLAB value give meaningfully different reconstructions? Compare on a medium-size proteome. |

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
| `replace_max_bound` | `False` | `True` | — | ⚠ **Issue #1** — see below |
| `min_flux` | `False` | `false` | — | ✓ |
| `loopless_good_reactions` | `True` | heuristic (±999 threshold) | — | ✓ Python is more correct |
| `exclude_reactions` | `None` | hardcoded ecModel logic | — | ✓ Python is more general |
| `max_attempts` | `100` | 100 | — | ✓ |
| `suppress_errors` | `False` | `false` | — | ✓ |

**Issue #1 — `replace_max_bound`:** MATLAB's `randomSampling` replaces very large upper
bounds (1000) with `Inf` by default so that the sampling polytope is not artificially
truncated at RAVEN's conventional big-M bound. Python defaults to `False` (no
replacement), which means enzyme-constrained or bounded models will have their sampling
polytope clipped at the big-M boundary — producing biased samples without any warning.
**Open question:** Does `True` actually produce better samples on real models? Run on
ecYeast9 with both values; compare sample variance and feasibility rate. See action item #1.

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
| `flux_eps` | `1e-6` | implicit `1e-8` | — | ⚠ **Issue #4** — see below |

**Issue #4 — `flux_eps`:** MATLAB uses an implicit tolerance of `1e-8` when classifying
reactions with near-zero flux as non-targets. Python's `1e-6` is 100× looser and may
classify reactions carrying negligible flux as targets or knockdowns.
**Open question:** Do the two values produce different target lists on realistic models?
Run FSEOF on iJO1366 and Yeast9 at `1e-6`, `1e-7`, `1e-8` and compare. See action item #4.

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
| `allow_excretion` | `True` | `false` | — | ⚠ **Issue #3** — see below |
| `no_rev_loops` | `False` | `false` | — | ✓ |
| `remove_dead_ends` | `True` | `true` | — | ✓ |
| `eps` | `1.0` | 1.0 | — | ✓ |
| `mip_gap` | `None` | `0.0004` | — | ⚠ **Issue #8** |
| `time_limit` | `None` | 5000 ms | — | ⚠ **Issue #8** |

**Issue #3 — `allow_excretion` inconsistency:** `get_init_model` defaults to `True` while
`run_init` and `run_ftinit` default to `False`. Since `run_init` calls `get_init_model`
and passes its own `allow_excretion=False` through, the `get_init_model` default is
overridden in normal usage — but direct calls to `get_init_model` behave differently.
This is confusing and error-prone.
**Open question:** Which value produces better reconstructions? Run `get_init_model`
directly on a test proteome with both values; compare feasibility and gene-essentiality
agreement. Once settled, align both wrappers. See action item #3.

**Issue #8 — MIP solver parameters:** MATLAB hardcodes `MIPGap=0.0004` and
`TimeLimit=5000 ms` per step inside the INIT algorithm. Python exposes these as `None`
(solver defaults), which on most solvers means `MIPGap=1e-4` and no time limit. The
MATLAB time limit prevents runaway solves on difficult models.
**Open question:** Do MATLAB's specific values give meaningfully different reconstruction
quality? Benchmark on a medium-size proteome. See action item #8.

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
**Open question:** On a proteome with known KO assignments, do the dropped hits at
`1e-5` correspond to true positives or noise? Report precision/recall at both thresholds.
See action item #2.

**Issue #5 — `threads=1`:** MATLAB detects available cores and uses them all. Python's
`threads=1` silently runs single-threaded, making BLAST/Diamond/HMMER dramatically
slower. This is a performance issue, not a correctness issue.
**Next step:** Verify results are identical across thread counts (BLAST is deterministic;
HMMER bit-scores may differ marginally between thread counts — check). If deterministic,
change default to `max(1, os.cpu_count() - 1)`. See action item #5.

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
no default cap, meaning large models can run indefinitely.
**Open question:** How long does localization actually take on Human-GEM? Profile and
check whether the MATLAB 15-min cap is ever hit. If solves routinely complete in well
under 15 min, leaving `None` is fine; if not, add a cap. See action item #7.

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
| `remove_genes` | `blocked_reactions` | `'remove'` | `false` (keep) | ⚠ **Issue #6** — see below |
| `remove_genes` | `remove_orphans` | `False` | N/A | ✓ |
| `find_duplicate_reactions` | `ignore_direction` | `True` | — | ✓ |
| `constrain_reversible_reactions` | `eps` | `1e-9` | — | ✓ |
| `add_transport_reactions` | `reversible` | `True` | — | ✓ |
| `add_transport_reactions` | `only_to_existing` | `True` | — | ✓ |
| `add_reactions_from_model` | `genes` | `False` | — | ✓ |

**Issue #6 — `remove_genes` blocked reactions:** MATLAB's `removeGenes` keeps reactions
that become gene-less (`removeBlockedRxns=false`); Python deletes them by default
(`blocked_reactions='remove'`). Both behaviours are defensible.
**Open question:** Which produces models with better predictive accuracy for gene
essentiality? Test on iJO1366 with a known essential-gene set; compare predictions
under `'remove'` vs `'keep'`. See action item #6.

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

Tests to run before deciding on defaults (in priority order):

- [ ] **#1** `replace_max_bound`: run `random_sampling` on ecYeast9 with `False` and `True`; compare sample variance and distribution of objective function values. Decide based on results.
- [ ] **#2** `evalue` (BLAST/Diamond): compare reconstruction coverage at `1e-5` vs `1e-4` on a proteome with a published KO annotation (e.g., a well-studied yeast or bacterium). Report precision/recall at each threshold.
- [ ] **#3** `allow_excretion` in `get_init_model`: run `run_init` with an expression dataset; compare model size and flux feasibility with `True` vs `False`. Pick the value that gives more biologically reasonable models, then make both wrappers consistent.
- [ ] **#4** `flux_eps` in `fseof`: run FSEOF on iJO1366 at `1e-6`, `1e-7`, `1e-8`; compare target reaction lists and check whether reactions at `1e-6`–`1e-8` range are biologically meaningful.
- [ ] **#5** `threads`: verify results are identical across thread counts (BLAST is deterministic, HMMER may vary); then change default to `max(1, os.cpu_count() - 1)`.
- [ ] **#6** `remove_genes` blocked policy: test on a model where a non-essential gene deletion is known; check whether removing vs keeping blocked reactions changes growth predictions.
- [ ] **#7** `time_limit` in `predict_localization`: run on Human-GEM and measure wall time; if solve completes in <15 min reliably, leave at `None` and document; otherwise add a default cap.
- [ ] **#8** `mip_gap` / `time_limit` in init: compare reconstruction outputs on a medium proteome with MATLAB's `0.0004`/`5000 ms` vs solver defaults; report solution quality and solve time.

## Evaluation checklist — remaining work

- [ ] **Numerical tolerances** (`eps`, `tol`, `flux_eps`, `reg`, `stoichiometry_tol`, `constrain_reversible_reactions.eps`) — run on ill-conditioned models to confirm they do not cause numeric issues.
- [ ] **MILP big-M** (`init.build.big_m`, `gapfilling.kumar_milp.big_m`) — verify against the largest observed flux bound in Yeast9/Human-GEM.
- [ ] **Homology thresholds** (`max_evalue`, `min_align_len`, `min_identity`, `cutoff`, `min_score_ratio_ko`, `min_score_ratio_g`) — benchmark on a proteome with known KO assignments.
- [ ] **KEGG assembly flags** (`keep_spontaneous`, `keep_undefined_stoich`, `keep_incomplete`, `keep_general`) — measure fraction of reactions retained/dropped; confirm MATLAB parity numerically.
- [ ] **Sampling parameters** (`thinning`, `warmup`, `n_samples`, `fixed_width_tol`) — run autocorrelation analysis on CHRR and ACHR chains on ecYeast9.
