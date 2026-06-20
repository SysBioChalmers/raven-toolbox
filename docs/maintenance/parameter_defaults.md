# Parameter defaults — inventory and evaluation plan

This document inventories every optional parameter (i.e. those with a current default value) in
raven-toolbox's public API, and provides a systematic methodology for deciding whether each
default is well-chosen.

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

## Parameter inventory

Columns: **Parameter** · **Current default** · **RAVEN MATLAB default** · **Paper/reference value** · **Status**

Status codes: ✓ validated · ⚠ change proposed · ? not yet evaluated

---

### `raven_toolbox.analysis.sampling` — `random_sampling`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `n_samples` | `1000` | 1000 | — | ? |
| `method` | `'achr'` | `'achr'` | — | ? |
| `seed` | `None` | — | — | ? |
| `thinning` | `100` | 100 | Kaufman & Smith 1998 | ? |
| `warmup` | `1000` | 1000 | — | ? |
| `fixed_width_tol` | `1e-7` | — | — | ? |
| `n_objectives` | `2` | 2 | Bordel et al. 2010 | ? |
| `replace_max_bound` | `False` | False | — | ? |
| `min_flux` | `False` | False | — | ? |
| `loopless_good_reactions` | `True` | True | — | ? |
| `max_attempts` | `100` | 100 | — | ? |
| `suppress_errors` | `False` | False | — | ? |

### `raven_toolbox.analysis.sampling` — `find_good_reactions`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `flux_tol` | `1e-9` | — | — | ? |
| `loopless` | `True` | True | — | ? |

### `raven_toolbox.analysis.flux_sampling` — `max_volume_ellipsoid`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `maxiter` | `150` | 150 | Zhang & Gao 2003; COBRA chrrSampler | ? |
| `tol` | `1e-6` | 1e-6 | Zhang & Gao 2003 | ? |
| `reg` | `1e-8` | 1e-8 | COBRA chrrSampler | ? |

### `raven_toolbox.analysis.fseof` — `fseof`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `n_steps` | `10` | 10 | Choi et al. 2010 | ? |
| `max_fraction` | `0.9` | 0.9 | Choi et al. 2010 | ? |
| `correlation_threshold` | `0.9` | 0.9 | Choi et al. 2010 | ? |
| `flux_eps` | `1e-6` | 1e-6 | — | ? |

---

### `raven_toolbox.annotation.sbo` — `add_sbo_terms`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `biomass_met_names` | `{'biomass','DNA','RNA','protein','carbohydrate','lipid','cofactor','ion'}` | yeast-GEM `addSBOterms.m` | yeast-GEM | ? |
| `biomass_met_suffixes` | `(' backbone', ' chain')` | yeast-GEM | yeast-GEM | ? |
| `biomass_rxn_name` | `'biomass pseudoreaction'` | yeast-GEM | yeast-GEM | ? |
| `ngam_rxn_name` | `'non-growth associated maintenance reaction'` | yeast-GEM | yeast-GEM | ? |
| `pseudoreaction_name_substrings` | `('pseudoreaction', 'SLIME rxn')` | yeast-GEM | yeast-GEM | ? |
| `only_last_reaction_for_pseudo` | `False` | True (bug) | — | ⚠ bug-compat flag — default intentionally differs from MATLAB |

---

### `raven_toolbox.init.init` — `run_init`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `prod_weight` | `0.5` | 0.5 | Agren et al. 2012 | ? |
| `allow_excretion` | `False` | False | Agren et al. 2012 | ? |
| `no_rev_loops` | `False` | False | — | ? |

### `raven_toolbox.init.build` — `get_init_model`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `isozyme_scoring` | `'max'` | `'max'` | Agren et al. 2012 | ? |
| `complex_scoring` | `'min'` | `'min'` | Agren et al. 2012 | ? |
| `no_gene_score` | `-2.0` | -2.0 | Agren et al. 2012 | ? |
| `prod_weight` | `0.5` | 0.5 | Agren et al. 2012 | ? |
| `allow_excretion` | `True` | True | — | ⚠ differs from `run_init` default (`False`) — inconsistency to resolve |
| `no_rev_loops` | `False` | False | — | ? |
| `remove_dead_ends` | `True` | True | — | ? |
| `eps` | `1.0` | 1.0 | — | ? |

### `raven_toolbox.init.score` — `gene_scores_from_expression` / `rna_gene_scores`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `factor` | `5.0` | 5.0 | Wang et al. 2012 (tINIT) | ? |
| `max_score` | `10.0` | 10.0 | Wang et al. 2012 | ? |
| `min_score` | `-5.0` | -5.0 | Wang et al. 2012 | ? |

### `raven_toolbox.init.score` — `score_reactions_from_genes`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `isozyme_scoring` | `'max'` | `'max'` | Agren et al. 2012 | ? |
| `complex_scoring` | `'min'` | `'min'` | Agren et al. 2012 | ? |
| `no_gene_score` | `-2.0` | -2.0 | Agren et al. 2012 | ? |

### `raven_toolbox.init.ftinit` — `ftinit` / `run_ftinit`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `series` | `'1+1'` | `'1+1'` | Gustafsson et al. 2023 | ? |
| `fill_gaps` | `True` | True | — | ? |
| `allow_excretion` | `False` | False | — | ? |
| `rem_pos_rev` | `False` | False | — | ? |

### `raven_toolbox.init.prep` — `prep_init_model` / `classify_reactions`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `ext_comp` | `'e'` | `'e'` | cobrapy convention | ? |
| `max_stoich_diff` | `25.0` | 25.0 | — | ? |
| `scale` | `True` | True | — | ? |

---

### `raven_toolbox.gapfilling.fast_lp` — `fill_gaps_fast_lp`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `epsilon` | `0.0001` | 0.0001 | Thiele et al. 2014 (fastGapFill) | ? |
| `variant` | `'fast'` | `'fast'` | — | ? |
| `verbose` | `True` | True | — | ? |

### `raven_toolbox.gapfilling.fill` — `connect_blocked_reactions`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `penalty` | `1.0` | 1.0 | — | ? |
| `allow_net_production` | `False` | False | — | ? |
| `eps` | `1.0` | 1.0 | — | ? |

### `raven_toolbox.gapfilling.kumar_milp` — `fill_gaps_kumar_milp`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `weights` | `(1.0, 2.0)` | — | Kumar et al. 2007 | ? |
| `big_m` | `1000.0` | 1000 | Kumar et al. 2007 | ? |
| `verbose` | `True` | — | — | ? |

---

### `raven_toolbox.reconstruction.homology.homology` — `get_model_from_homology`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `bidirectional` | `True` | True | — | ? |
| `best_hits_only` | `False` | False | — | ? |
| `map_direction` | `'new_to_old'` | `'new_to_old'` | — | ? |
| `score` | `'bitscore'` | `'bitscore'` | BLAST convention | ? |
| `complex_policy` | `'flag'` | `'flag'` | — | ? |
| `only_genes_in_models` | `False` | False | — | ? |
| `max_evalue` | `1e-30` | 1e-30 | — | ? |
| `min_align_len` | `200` | 200 | — | ? |
| `min_identity` | `40` | 40 | — | ? |

### `raven_toolbox.reconstruction.homology.blast` — `run_blast` / `run_diamond`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `evalue` | `1e-5` | 1e-5 | BLAST convention | ? |
| `threads` | `1` | 1 | — | ⚠ single-threaded is conservative; consider `os.cpu_count()` or a user-settable global |
| `sensitivity` | `'--more-sensitive'` (diamond) | — | Diamond docs | ? |

### `raven_toolbox.reconstruction.kegg.query` — `assign_kos` / `get_kegg_model_from_sequences`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `cutoff` | `1e-30` | 1e-30 | — | ? |
| `min_score_ratio_ko` | `0.3` | 0.3 | — | ? |
| `min_score_ratio_g` | `0.9` | 0.9 | — | ? |

### `raven_toolbox.reconstruction.kegg.hmm` — `build_ko_hmm` / `build_hmm_library`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `seq_identity` | `0.9` | 0.9 | CD-HIT recommendation | ? |
| `threads` | `1` | 1 | — | ⚠ same single-thread issue as BLAST |
| `fast` | `True` | True | MAFFT `--parttree` flag | ? |
| `concatenate` | `True` | True | — | ? |

### `raven_toolbox.reconstruction.kegg.assemble` / `.organism` / `.query` — model assembly

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `keep_spontaneous` | `True` | True | — | ? |
| `keep_undefined_stoich` | `True` | True | — | ? |
| `keep_incomplete` | `True` | True | — | ? |
| `keep_general` | `False` | False | — | ? |

---

### `raven_toolbox.localization.predict` — `predict_localization`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `default_compartment` | `'c'` | `'c'` | cobrapy convention | ? |
| `transport_cost` | `0.5` | 0.5 | — | ? |
| `multi_compartment_penalty` | `0.5` | 0.5 | — | ? |
| `apply` | `True` | True | — | ? |

---

### `raven_toolbox.manipulation` — structural transforms

| Function | Parameter | Default | Status |
|---|---|---|---|
| `add_reactions_from_equations` | `mets_by` | `'id'` | ? |
| `add_reactions_from_equations` | `allow_new_mets` | `True` | ? |
| `add_reactions_from_equations` | `allow_new_genes` | `True` | ? |
| `merge_compartments` | `drop_single_metabolite_reactions` | `True` | ? |
| `merge_compartments` | `deduplicate_reactions` | `True` | ? |
| `merge_models` | `match_by` | `'name'` | ? |
| `remove_genes` | `blocked_reactions` | `'remove'` | ? |
| `remove_genes` | `remove_orphans` | `False` | ? |
| `find_duplicate_reactions` | `ignore_direction` | `True` | ? |
| `constrain_reversible_reactions` | `eps` | `1e-9` | ? |
| `add_transport_reactions` | `reversible` | `True` | ? |
| `add_transport_reactions` | `only_to_existing` | `True` | ? |
| `add_transport_reactions` | `id_prefix` | `'tr_'` | ? |
| `add_reactions_from_model` | `genes` | `False` | ? |

---

### `raven_toolbox.io` — export

| Function | Parameter | Default | Status |
|---|---|---|---|
| `export_for_git` | `formats` | `('yml','xml','mat','xlsx')` | ? |
| `export_for_git` | `sub_dirs` | `True` | ? |
| `export_to_excel` | `sort_ids` | `False` | ? |
| `write_yaml_model` | `sort_ids` | `False` | ? |
| `export_model_to_sif` | `graph_type` | `'rc'` | ? |

---

### `raven_toolbox.comparison.diff` — `diff_models`

| Parameter | Default | RAVEN MATLAB | Reference | Status |
|---|---|---|---|---|
| `stoichiometry_tol` | `1e-9` | — | — | ? |
| `max_per_category` | `50` | — | — | ? |

---

### `raven_toolbox.tasks.check` — task checking

| Function | Parameter | Default | Status |
|---|---|---|---|
| `check_tasks` | `close_boundaries` | `True` | ? |
| `find_task_essential_reactions` | `close_boundaries` | `True` | ? |
| `find_task_essential_reactions` | `tol` | `1e-8` | ? |

---

## Known inconsistencies to resolve

| # | Description | Functions involved |
|---|---|---|
| 1 | `allow_excretion` default is `True` in `get_init_model` but `False` in `run_init` and `run_ftinit` — likely a bug in `get_init_model`. | `init.build`, `init.init`, `init.ftinit` |
| 2 | `threads=1` is safe but slow for `run_blast`, `run_diamond`, `build_ko_hmm`, `build_hmm_library`. Consider defaulting to `os.cpu_count()` or a package-level config. | `reconstruction.homology.blast`, `reconstruction.kegg.hmm` |
| 3 | `verbose=True` on gap-filling functions is noisy in library use. Consistent with RAVEN MATLAB but worth reconsidering if raven-toolbox is used as a library. | `gapfilling.fast_lp`, `gapfilling.kumar_milp`, `gapfilling.topological` |

---

## Evaluation checklist (per parameter group)

- [ ] **Numerical tolerances** (`eps`, `tol`, `flux_eps`, `reg`, `stoichiometry_tol`, `constrain_reversible_reactions.eps`) — confirm units, confirm they are not conflated with FBA solver tolerances, test on ill-conditioned models.
- [ ] **MILP big-M** (`init.build.big_m`, `gapfilling.kumar_milp.big_m`) — a too-small big-M silently makes constraints inactive; test against the largest observed flux bound in Yeast9/Human-GEM.
- [ ] **Scoring weights** (`no_gene_score`, `prod_weight`, `factor`, `max_score`, `min_score`, `weights`) — compare to original tINIT paper values and to any published sensitivity studies.
- [ ] **Homology thresholds** (`max_evalue`, `min_align_len`, `min_identity`, `cutoff`, `min_score_ratio_ko`, `min_score_ratio_g`) — benchmark against a proteome with known KO assignments.
- [ ] **KEGG assembly flags** (`keep_spontaneous`, `keep_undefined_stoich`, `keep_incomplete`, `keep_general`) — measure the fraction of reactions retained/dropped in a test reconstruction; confirm MATLAB RAVEN parity.
- [ ] **Sampling parameters** (`thinning`, `warmup`, `n_samples`, `fixed_width_tol`) — run autocorrelation analysis on CHRR and ACHR chains to validate that defaults give well-mixed chains on ecYeast9.
