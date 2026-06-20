# Parameter defaults — inventory and evaluation plan

This document inventories every optional parameter (i.e. those with a current default value) in
raven-toolbox's public API, and provides a systematic methodology for deciding whether each
default is well-chosen.

Last evaluated: 2026-06-20 against MATLAB RAVEN at
`C:\Work\GitHub\raven-docs\RAVEN` (commit on tracked branch).

---

## Evaluation methodology

A default value is *well-chosen* when a user who does not read the docstring gets a
result that is correct and useful for the most common case. The following criteria
apply in rough priority order:

1. **MATLAB RAVEN parity.** Where a function ports a MATLAB RAVEN function, the Python
   default must reproduce RAVEN's behaviour unless there is a documented, intentional
   improvement.
2. **Literature anchor.** Algorithm-specific numerical parameters (tolerances, scoring
   weights, iteration counts) should match the value used in the original paper or the
   most-cited open-source implementation (cobrapy, COBRA Toolbox).
3. **Fail-safe on realistic models.** Run the function with the default on at least one
   large model (Yeast9, Human-GEM) and one small model (iJO1366 or similar). The default
   must not crash, produce NaN/Inf, or return an obviously wrong answer.
4. **Sensitivity envelope.** If varying the parameter by ±50 % changes the result by
   more than 5–10 %, the default should be documented with a stronger note, and the
   sensitivity test result should be recorded here.
5. **User expectation alignment.** Prefer values that match what a competent user would
   supply without thinking (e.g., `verbose=True` for long-running MILP, `sort_ids=False`
   for round-trip-safe export).
6. **No None-surprises.** `None` defaults are fine for optional features but should
   never silently change algorithmic behaviour; document the fallback clearly.

### Evaluation workflow per parameter

```
1. Read the current docstring — does it explain *why* this value?
2. Find the MATLAB RAVEN equivalent (if any) and compare values.
3. Find the source paper / reference implementation and compare values.
4. Run the benchmark suite (or a quick FBA sanity check) with default and ±50 % variants.
5. Record finding in the "Status" column below: ✓ validated / ⚠ needs change / ? unknown.
6. If ⚠: open an issue, propose a new default and rationale, update this table.
```

---

## Summary of issues found

Issues identified during the 2026-06-20 evaluation pass. Each is expanded in the
relevant section below.

| # | Severity | Parameter | Finding |
|---|---|---|---|
| 1 | **High** | `random_sampling.replace_max_bound` | Python default `False`; MATLAB default `True` — migration hazard |
| 2 | **High** | `run_blast / run_diamond.evalue` | Python `1e-5`; MATLAB `1e-4` — Python 10× more stringent, misses valid hits |
| 3 | **Medium** | `get_init_model.allow_excretion` | `True` here, `False` in `run_init` and `run_ftinit` — internal inconsistency |
| 4 | **Medium** | `fseof.flux_eps` | Python `1e-6`; MATLAB implicit `1e-8` — Python includes near-zero "targets" |
| 5 | **Medium** | `run_blast / run_diamond / run_hmmsearch.threads` | Python `1`; MATLAB auto-detects cores — silent performance trap |
| 6 | **Medium** | `remove_genes.blocked_reactions` | Python removes blocked reactions; MATLAB keeps them |
| 7 | **Low** | `predict_localization.time_limit` | MATLAB caps at 15 min; Python has no default cap |
| 8 | **Low** | `run_init / run_ftinit` MIP gap | MATLAB hardcodes `0.0004`; Python `mip_gap=None` (solver default) |

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
**Proposed fix:** change Python default to `True` and document why, matching MATLAB.

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
**Proposed fix:** tighten to `1e-8` (or document the intentional difference if the looser
value was chosen for numerical robustness with certain solvers).

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
**Proposed fix:** change `get_init_model` default to `False` to match the higher-level
wrappers.

**Issue #8 — MIP solver parameters:** MATLAB hardcodes `MIPGap=0.0004` and
`TimeLimit=5000 ms` per step inside the INIT algorithm. Python exposes these as `None`
(solver defaults), which on most solvers means `MIPGap=1e-4` and no time limit. The
MATLAB time limit prevents runaway solves on difficult models. Document the MATLAB
values explicitly as recommended starting points.

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
uses `1e-5`. Python is 10× more stringent, meaning it will silently drop valid
low-confidence homologs that MATLAB would include. This changes reconstruction results
for distantly related organisms.
**Proposed fix:** change Python default to `1e-4` to match MATLAB, or document the
intentional tightening with a clear rationale and migration note.

**Issue #5 — `threads=1`:** MATLAB detects available cores and uses them all. Python's
`threads=1` silently runs single-threaded, which makes BLAST/Diamond/HMMER
dramatically slower on multi-core machines. This is a silent performance trap with no
warning.
**Proposed fix:** default to `max(1, os.cpu_count() - 1)` or add a package-level
`raven_toolbox.config.threads` setting, and warn if `threads=1` is used on a machine
with more cores.

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
**Proposed fix:** document the MATLAB value in the docstring as a recommended starting
point; optionally default to `900` (seconds) for large models.

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
(`blocked_reactions='remove'`). Both behaviours are defensible but the difference will
silently change model size for users porting MATLAB workflows.
**Proposed fix:** add a migration note in the docstring and consider whether `'keep'`
should be the default to match MATLAB.

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

Issues to resolve (in priority order):

- [ ] **#1** `replace_max_bound`: change default to `True`; add migration note.
- [ ] **#2** `evalue` (BLAST/Diamond): align with MATLAB (`1e-4`) or document the tightening.
- [ ] **#3** `allow_excretion` in `get_init_model`: change to `False`.
- [ ] **#4** `flux_eps` in `fseof`: tighten to `1e-8` or document the looser value.
- [ ] **#5** `threads`: default to `os.cpu_count() - 1` or add package-level config.
- [ ] **#6** `remove_genes` blocked policy: add migration note; reconsider default.
- [ ] **#7** `time_limit` in `predict_localization`: document MATLAB's 15-min cap.
- [ ] **#8** `mip_gap` / `time_limit` in init: document MATLAB's `0.0004` / 5 s as guidance.

## Evaluation checklist — remaining work

- [ ] **Numerical tolerances** (`eps`, `tol`, `flux_eps`, `reg`, `stoichiometry_tol`, `constrain_reversible_reactions.eps`) — run on ill-conditioned models to confirm they do not cause numeric issues.
- [ ] **MILP big-M** (`init.build.big_m`, `gapfilling.kumar_milp.big_m`) — verify against the largest observed flux bound in Yeast9/Human-GEM.
- [ ] **Homology thresholds** (`max_evalue`, `min_align_len`, `min_identity`, `cutoff`, `min_score_ratio_ko`, `min_score_ratio_g`) — benchmark on a proteome with known KO assignments.
- [ ] **KEGG assembly flags** (`keep_spontaneous`, `keep_undefined_stoich`, `keep_incomplete`, `keep_general`) — measure fraction of reactions retained/dropped; confirm MATLAB parity numerically.
- [ ] **Sampling parameters** (`thinning`, `warmup`, `n_samples`, `fixed_width_tol`) — run autocorrelation analysis on CHRR and ACHR chains on ecYeast9.
