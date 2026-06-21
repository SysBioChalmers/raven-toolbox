# Parameter benchmark index — master to-do list

All parameters with non-trivial defaults, grouped by function. Each entry shows the
current Python default, the MATLAB RAVEN equivalent (where one exists), and the
recommended action based on empirical benchmarks.

Status codes: **✓ keep** (tested, correct) · **⚠ change** (requires code edit) ·
**? untested** (benchmark not yet run) · **— Python-only** (no MATLAB counterpart)

Benchmark date: 2026-06-20. Models: yeast-GEM 4102 rxns, iJO1366 2583 rxns,
e_coli_core 95 rxns, synthetic toy models. Binaries: BLAST 2.17.0.

---

## Flux sampling

| Function | Parameter | Python | MATLAB | Action |
|---|---|---|---|---|
| `random_sampling` | `n_samples` | `1000` | 1000 | ✓ keep |
| `random_sampling` | `method` | `'achr'` | `'random_objective'` | ✓ keep (ACHR is preferred) |
| `random_sampling` | `seed` | `None` | unseeded | ✓ keep |
| `random_sampling` | `thinning` | `100` | N/A | ⚠ keep value but add docstring warning: yeast-GEM gives ESS≈12 from 300 samples; ~12 effective samples. Use `n_samples≥2600` or switch sampler for genome scale. |
| `random_sampling` | `warmup` | `1000` | N/A | ✓ keep (cobrapy default) |
| `random_sampling` | `n_objectives` | `2` | 2 | ✓ keep (Bordel 2010) |
| `random_sampling` | `replace_max_bound` | `False` | `True` | ✓ keep `False` (MATLAB `True` → solver unbounded) |
| `random_sampling` | `min_flux` | `False` | `false` | ✓ keep |
| `random_sampling` | `loopless_good_reactions` | `True` | heuristic ±999 | ✓ keep (more correct) |
| `random_sampling` | `max_attempts` | `100` | 100 | ✓ keep |
| `find_good_reactions` | `flux_tol` | `1e-9` | — | ✓ keep |
| `max_volume_ellipsoid` | `maxiter` | `150` | 150 | ✓ keep (Zhang & Gao 2003) |
| `max_volume_ellipsoid` | `tol` | `1e-6` | 1e-6 | ✓ keep |
| `max_volume_ellipsoid` | `reg` | `1e-8` | 1e-8 | ✓ keep |

**Benchmark file:** [sampling.md](sampling.md)

---

## FSEOF

| Function | Parameter | Python | MATLAB | Action |
|---|---|---|---|---|
| `fseof` | `n_steps` | `10` | 10 | ✓ keep (Choi 2010; n=5 gives marginally fewer targets) |
| `fseof` | `max_fraction` | `0.9` | 0.9 | ✓ keep (Choi 2010; 0.5 picks up spurious targets) |
| `fseof` | `correlation_threshold` | `0.9` | 0.9 | ✓ keep (Choi 2010; 0.7 adds 4 spurious amplification targets on iJO1366) |
| `fseof` | `flux_eps` | `1e-6` | implicit `1e-8` | ✓ keep `1e-6` (MATLAB `1e-8` picks up solver noise as false-positive targets) |

**Benchmark file:** [fseof.md](fseof.md)

---

## INIT (tINIT/ftINIT)

| Function | Parameter | Python | MATLAB | Action |
|---|---|---|---|---|
| `run_init` | `prod_weight` | `0.5` | 0.5 | ✓ keep (Agren 2012) |
| `run_init` | `allow_excretion` | `False` | `false` | ✓ keep |
| `run_init` | `eps` | `1.0` | 1.0 | ✓ keep |
| `run_init` | `mip_gap` | `None` | `0.0004` | ✓ keep `None`; document MATLAB value as guidance |
| `run_init` | `time_limit` | `None` | 5000 ms | ✓ keep `None`; document MATLAB value as guidance |
| `get_init_model` | `allow_excretion` | `True` | `false` | ⚠ **change to `False`** (inconsistency with `run_init`) |
| `get_init_model` | `eps` | `1.0` | 1.0 | ✓ keep |
| `run_ftinit` | `series` | `'1+1'` | `'1+1'` | ✓ keep (Gustafsson 2023) |
| `run_ftinit` | `force_on` | `0.1` | 0.1 | ✓ keep |
| `run_ftinit` | `big_m` | `100.0` | 100 | ✓ keep (intentional LP tightener; see `init.md`) |
| `run_ftinit` | `mip_gap` | `None` | `0.0004` | ✓ keep `None`; document MATLAB value |
| `run_ftinit` | `time_limit` | `None` | 5000 ms | ✓ keep `None`; document MATLAB value |
| `gene_scores_from_expression` | `factor` | `5.0` | 5 | ✓ keep (Wang 2012) |
| `gene_scores_from_expression` | `max_score` | `10.0` | 10 | ✓ keep |
| `gene_scores_from_expression` | `min_score` | `-5.0` | -5 | ✓ keep |
| `score_reactions_from_genes` | `isozyme_scoring` | `'max'` | `'max'` | ✓ keep |
| `score_reactions_from_genes` | `complex_scoring` | `'min'` | `'min'` | ✓ keep |
| `score_reactions_from_genes` | `no_gene_score` | `-2.0` | -2 | ✓ keep |
| `classify_reactions` | `ext_comp` | `'e'` | `'e'` | ✓ keep |
| `classify_reactions` | `max_stoich_diff` | `25.0` | ~25 | ✓ keep |

**Benchmark file:** [init.md](init.md)

---

## Gapfilling

| Function | Parameter | Python | MATLAB | Action |
|---|---|---|---|---|
| `fill_gaps_fast_lp` | `epsilon` | `0.0001` | N/A | ✓ keep (fastGapFill paper) |
| `connect_blocked_reactions` | `penalty` | `1.0` | N/A | ✓ keep |
| `connect_blocked_reactions` | `allow_net_production` | `False` | `false` | ✓ keep |
| `connect_blocked_reactions` | `eps` | `1.0` | N/A | ✓ keep; edge case: lower if nutrient supply < 1 mmol/gDW/h |
| `fill_gaps_kumar_milp` | `weights` | `(1.0, 2.0)` | N/A | ✓ keep (Kumar 2007; reversal preferred over addition at 2:1 ratio, confirmed) |
| `fill_gaps_kumar_milp` | `big_m` | `1000.0` | N/A | ✓ keep (matches RAVEN ±1000 bounds; increase for enzyme-constrained models) |

**Benchmark file:** [gapfilling.md](gapfilling.md)

---

## Homology-based reconstruction

| Function | Parameter | Python | MATLAB | Action |
|---|---|---|---|---|
| `run_blast` | `evalue` | `1e-5` | `1e-4` | ✓ keep `1e-5`; see `reconstruction_homology.md` for BLAST test |
| `run_blast` | `threads` | `1` | all cores | ⚠ **change to `max(1, os.cpu_count() - 1)`** (performance; BLAST deterministic across threads) |
| `run_diamond` | `evalue` | `1e-5` | `1e-4` | ✓ keep `1e-5` (same rationale) |
| `run_diamond` | `threads` | `1` | all cores | ⚠ **change to `max(1, os.cpu_count() - 1)`** |
| `run_diamond` | `sensitivity` | `'--more-sensitive'` | `'--more-sensitive'` | ✓ keep |
| `get_model_from_homology` | `bidirectional` | `True` | `true` | ✓ keep |
| `get_model_from_homology` | `max_evalue` | `1e-30` | `1e-30` | ? untested |
| `get_model_from_homology` | `min_align_len` | `200` | 200 | ? untested |
| `get_model_from_homology` | `min_identity` | `40` | 40 | ? untested |

**Benchmark file:** [reconstruction_homology.md](reconstruction_homology.md)

---

## KEGG-based reconstruction

| Function | Parameter | Python | MATLAB | Action |
|---|---|---|---|---|
| `run_hmmsearch` | `threads` | `1` | all cores | ⚠ **change to `max(1, os.cpu_count() - 1)`** |
| `build_ko_hmm` | `seq_identity` | `0.9` | 0.9 | ✓ keep (CD-HIT recommendation) |
| `build_ko_hmm` | `threads` | `1` | all cores | ⚠ **change** (same as above) |
| `assign_kos` | `cutoff` | `1e-30` | `1e-30` | ? untested |
| `assign_kos` | `min_score_ratio_ko` | `0.3` | 0.3 | ? untested |
| `assign_kos` | `min_score_ratio_g` | `0.9` | 0.9 | ? untested |
| `get_kegg_model_*` | `keep_spontaneous` | `True` | `true` | ✓ keep |
| `get_kegg_model_*` | `keep_undefined_stoich` | `True` | `true` | ✓ keep |
| `get_kegg_model_*` | `keep_incomplete` | `True` | `true` | ✓ keep |
| `get_kegg_model_*` | `keep_general` | `False` | `false` | ✓ keep |

**Benchmark file:** [reconstruction_kegg.md](reconstruction_kegg.md)

---

## Localization

| Function | Parameter | Python | MATLAB | Action |
|---|---|---|---|---|
| `predict_localization` | `default_compartment` | `'c'` | required arg | ✓ keep (better UX) |
| `predict_localization` | `transport_cost` | `0.5` | 0.5 | ✓ keep |
| `predict_localization` | `time_limit` | `None` | 900 s | ✓ keep `None`; document 900 s cap for Human-GEM scale |
| `predict_localization` | `mip_gap` | `None` | N/A | ✓ keep |

**Benchmark file:** [localization.md](localization.md)

---

## Model manipulation

| Function | Parameter | Python | MATLAB | Action |
|---|---|---|---|---|
| `remove_genes` | `blocked_reactions` | `'remove'` | `'keep'` (false) | ✓ keep `'remove'` (MATLAB `'keep'` breaks essentiality predictions) |
| `remove_genes` | `remove_orphans` | `False` | N/A | ✓ keep |
| `constrain_reversible_reactions` | `eps` | `1e-9` | — | ✓ keep |
| `merge_models` | `match_by` | `'name'` | `'metNames'` | ✓ keep (equivalent) |
| `add_reactions_from_equations` | `mets_by` | `'id'` | — | ✓ keep |
| `find_duplicate_reactions` | `ignore_direction` | `True` | — | ✓ keep |
| `add_transport_reactions` | `reversible` | `True` | — | ✓ keep |
| `add_transport_reactions` | `only_to_existing` | `True` | — | ✓ keep |

**Benchmark file:** [manipulation.md](manipulation.md)

---

## Tasks

| Function | Parameter | Python | MATLAB | Action |
|---|---|---|---|---|
| `check_tasks` | `close_boundaries` | `True` | implied | ✓ keep |
| `find_task_essential_reactions` | `close_boundaries` | `True` | implied | ✓ keep |
| `find_task_essential_reactions` | `tol` | `1e-8` | — | ✓ keep |

**Benchmark file:** [tasks.md](tasks.md)

---

## IO / export

| Function | Parameter | Python | MATLAB | Action |
|---|---|---|---|---|
| `export_to_excel` | `sort_ids` | `False` | implicit unsorted | ✓ keep |
| `write_yaml_model` | `sort_ids` | `False` | N/A | ✓ keep |
| `export_model_to_sif` | `graph_type` | `'rc'` | N/A | ✓ keep |
| `export_for_git` | `formats` | all four | N/A | ✓ keep |

---

## Summary of required code changes

| Change | File | Priority |
|---|---|---|
| `get_init_model` `allow_excretion` default: `True` → `False` | `src/raven_toolbox/init/build.py` | Medium |
| `run_blast` / `run_diamond` `threads` default: `1` → `max(1, os.cpu_count()-1)` | `src/raven_toolbox/reconstruction/homology/blast.py` | High |
| `run_hmmsearch` `threads` default: `1` → `max(1, os.cpu_count()-1)` | `src/raven_toolbox/reconstruction/kegg/query.py` | High |
| `build_ko_hmm` `threads` default: `1` → `max(1, os.cpu_count()-1)` | `src/raven_toolbox/reconstruction/kegg/hmm.py` | High |
| Docstring: `time_limit` note in `predict_localization` | `src/raven_toolbox/localization/predict.py` | Low |
| Docstring: `mip_gap`/`time_limit` note in INIT functions | `src/raven_toolbox/init/init.py`, `ftinit.py`, `build.py` | Low |

## Parameters needing further benchmarks

- `get_model_from_homology` thresholds (`max_evalue`, `min_align_len`, `min_identity`): require a proteome with known KO/ortholog assignments and a reference reconstruction to compute precision/recall
- `assign_kos` score ratios (`min_score_ratio_ko`, `min_score_ratio_g`, `cutoff`): same requirement
- Sampling `thinning`/`warmup` autocorrelation on ACHR: yeast-GEM multi-chain analysis pending
- INIT `mip_gap` genome-scale: needs real expression data to distinguish solution quality at different gaps

---

```{toctree}
:hidden:

fseof
gapfilling
init
localization
manipulation
reconstruction_homology
reconstruction_kegg
sampling
tasks
```
