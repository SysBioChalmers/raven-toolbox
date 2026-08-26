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
| `fseof` | `flux_eps` | `1e-6` | implicit `1e-8` | ⚠ **unify at `1e-6`** (MATLAB changes — needs exposing as a tunable first) — measured: `1e-8` catches 21 solver-noise false positives (std≈5e-7, below accumulated Gurobi feasibility tol) |

**Benchmark file:** [fseof.md](fseof.md)

---

## INIT (tINIT/ftINIT)

| Function | Parameter | Python | MATLAB | Action |
|---|---|---|---|---|
| `run_init` | `prod_weight` | `0.5` | 0.5 | ✓ keep (Agren 2012) |
| `run_init` | `allow_excretion` | `False` | `false` | ✓ keep |
| `run_init` | `eps` | `1.0` | 1.0 | ✓ keep |
| `run_init` | `mip_gap` | `None` | `0.0004` | ⚠ not yet resolved — see [parity decisions](#cross-toolbox-parity-decisions) |
| `run_init` | `time_limit` | `None` | 5000 ms | ⚠ not yet resolved — see [parity decisions](#cross-toolbox-parity-decisions) |
| `get_init_model` | `allow_excretion` | `True` | `false` | ⚠ **change to `False`** (inconsistency with `run_init`) |
| `get_init_model` | `eps` | `1.0` | 1.0 | ✓ keep |
| `run_ftinit` | `series` | `'1+1'` | `'1+1'` | ✓ keep (Gustafsson 2023) |
| `run_ftinit` | `force_on` | `0.1` | 0.1 | ✓ keep |
| `run_ftinit` | `big_m` | `100.0` | 100 | ✓ keep (intentional LP tightener; see `init.md`) |
| `run_ftinit` | `mip_gap` | `None` | `0.0004` | ⚠ not yet resolved — see [parity decisions](#cross-toolbox-parity-decisions) |
| `run_ftinit` | `time_limit` | `None` | 5000 ms | ⚠ not yet resolved — see [parity decisions](#cross-toolbox-parity-decisions) |
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
| `run_blast` | `evalue` | `1e-4` | `1e-4` | ✓ implemented — unified with MATLAB |
| `run_blast` | `threads` | `max(1, cpu_count-1)` | all cores | ✓ implemented (BLAST deterministic across threads; ~1.9–4× speedup measured) |
| `run_diamond` | `evalue` | `1e-4` | `1e-4` | ✓ implemented — unified with MATLAB |
| `run_diamond` | `threads` | `max(1, cpu_count-1)` | all cores | ✓ implemented (same change, applied) |
| `run_diamond` | `sensitivity` | `'--more-sensitive'` | `'--more-sensitive'` | ✓ keep |
| `get_model_from_homology` | `bidirectional` | `True` | `true` | ✓ keep |
| `get_model_from_homology` | `max_evalue` | `1e-30` | `1e-30` | ✓ keep (confirmed inert 1e-4…1e-50; measured against KEGG+OMA) |
| `get_model_from_homology` | `min_align_len` | `100` | 200 | ✓ implemented on Python side — measured; MATLAB side still pending (see MATLAB parity table) |
| `get_model_from_homology` | `min_identity` | `40` | 40 | ✓ keep (confirmed optimal against KEGG+OMA, β=0.5) |

**Benchmark file:** [reconstruction_homology.md](reconstruction_homology.md)

---

## KEGG-based reconstruction

| Function | Parameter | Python | MATLAB | Action |
|---|---|---|---|---|
| `run_hmmsearch` | `threads` | `max(1, cpu_count-1)` | all cores | ✓ implemented |
| `build_ko_hmm` | `seq_identity` | `0.9` | 0.9 | ✓ keep (CD-HIT recommendation) |
| `build_ko_hmm` | `threads` | `max(1, cpu_count-1)` | all cores | ✓ implemented |
| `assign_kos` | `cutoff` | `1e-30` | `1e-50` | ⚠ **unify at `1e-30`** (MATLAB changes) — measured, see `kegg_hmm_cutoff_calibration.md` |
| `assign_kos` | `min_score_ratio_ko` | `0.3` | 0.3 | ✓ keep (confirmed empirically inert; retained for RAVEN parity) |
| `assign_kos` | `min_score_ratio_g` | `0.9` | 0.8 | ⚠ **unify at `0.9`** (MATLAB changes) — the real precision lever, measured |
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
| `predict_localization` | `time_limit` | `None` | 900 s | ⚠ **tentatively unify at `None`** (MATLAB changes) — Medium confidence; see [parity decisions](#cross-toolbox-parity-decisions) |
| `predict_localization` | `mip_gap` | `None` | N/A | ✓ keep |

**Benchmark file:** [localization.md](localization.md)

---

## Model manipulation

| Function | Parameter | Python | MATLAB | Action |
|---|---|---|---|---|
| `remove_genes` | `blocked_reactions` | `'remove'` | `'keep'` (false) | ⚠ **unify at `'remove'` semantics** (MATLAB changes) — measured on e_coli_core: `'keep'` gives false-positive growth after an essential gene is removed |
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

## Cross-toolbox parity decisions

Policy: where Python and MATLAB RAVEN disagree on a default, pick one value for
both rather than let them drift — except where the divergence is forced by a real
implementation difference (different solver stack, different algorithm, different
data schema), in which case forcing identical values would change correct
behaviour into incorrect behaviour. Every row below either says which side has to
change, or explains why it stays split.

### Unify — one value, both toolboxes converge

| Parameter | Python now | MATLAB now | Unify at | Side that changes | Confidence | Evidence |
|---|---|---|---|---|---|---|
| `run_blast`/`run_diamond` `evalue` | `1e-4` ✓ done | `1e-4` | **`1e-4`** | ~~Python~~ done | High | No downstream effect measured either way (`get_model_from_homology`'s `1e-30` dominates); matches `develop`'s independent PR #91 call |
| `get_model_from_homology.min_align_len` | `100` ✓ done | `200` | **`100`** | Python done; MATLAB pending | High | Directly measured — [homology_cutoff_calibration.md](../../studies/homology_cutoff_calibration.md) |
| `assign_kos.cutoff` | `1e-30` | `1e-50` | **`1e-30`** | MATLAB | High | Directly measured — [kegg_hmm_cutoff_calibration.md](../../studies/kegg_hmm_cutoff_calibration.md) |
| `assign_kos.min_score_ratio_g` | `0.9` | `0.8` | **`0.9`** | MATLAB | High | Directly measured, same study |
| `fseof.flux_eps` | `1e-6` | implicit `1e-8` | **`1e-6`** | MATLAB (needs exposing as a tunable first) | High | Measured — `1e-8` catches 21 solver-noise false positives (std≈5e-7, below Gurobi's feasibility tolerance accumulated genome-wide); see `fseof.md` |
| `remove_genes.blocked_reactions` | `'remove'` | `'keep'` | **`'remove'` semantics** | MATLAB | High | Measured on e_coli_core — `'keep'` gives false-positive growth after an essential gene is removed; see `manipulation.md` |
| `get_init_model.allow_excretion` | `True` (pending) | `False` | **`False`** | Python | High | Already established — zero effect at default `prod_weight`, pure inconsistency; see `manipulation.md` |
| `predict_localization.time_limit` | `None` | `900 s` | **`None`**, tentatively | MATLAB | Medium | Wall-clock caps are hardware-relative, not portable; Python's `None` already validated on the primary dev-scale model. No cross-solver test run for this function specifically. |

### Gated — not a simple value choice, resolve the underlying issue first

| Parameter | Why this isn't just "pick one value" |
|---|---|
| `run_init`/`run_ftinit` `mip_gap` | [`init_solver_benchmark.md`](../../studies/init_solver_benchmark.md) found GLPK doesn't converge in 1h+ and HiGHS doesn't work with cobra in this stack at all — Gurobi is the only viable genome-scale backend today. `mip_gap=None` resolves to Gurobi's own default (`~1e-4`), which is already *tighter* than MATLAB's `0.0004`, so there's no known correctness gap on the one backend that actually works. Revisit once the pending genome-scale expression-data test (below) exists, and once the GLPK/HiGHS upstream issues are fixed enough that solver choice is a real option. |
| `run_init`/`run_ftinit` `time_limit` | Same study found GLPK **does not honor `configuration.timeout` at all** — a `time_limit` value is silently ignored on that backend regardless of what default is picked, so unifying the number doesn't unify the behaviour. On Gurobi, `None` (uncapped) is what's actually been run successfully at genome scale ([humangem_validation.md](../../studies/humangem_validation.md)); MATLAB's `5000 ms` cap has not been shown to be necessary or sufficient there. Needs the genome-scale test below before picking a firm value. |

### Keep different — forcing identical values would break correctness

| Parameter | Python | MATLAB | Why |
|---|---|---|---|
| `random_sampling.method` | `'achr'` | `'random_objective'` | Different algorithms, not a flag. ACHR is methodologically preferred (uniform interior sampling vs vertex-biased). Unifying means porting ACHR into MATLAB RAVEN — a real implementation project, tracked separately, not a default change. |
| `random_sampling.replace_max_bound` | `False` | `True` | Applying MATLAB's `True` inside cobrapy/optlang makes the sampler unbounded on standard RAVEN-convention models (measured — see `sampling.md`). Whether MATLAB's own solver path handles `True` safely wasn't tested here; either way this is a solver-stack constraint, not a preference. |
| `random_sampling.loopless_good_reactions` | `True` (loopless FVA) | heuristic (exclude FVA ≥ 999) | Different techniques, not a value. Python's is more correct; MATLAB's heuristic over-excludes reactions that legitimately reach capacity. Porting the proper technique to MATLAB is a bigger project than a default flip. |
| `merge_models.match_by` | `'name'` | `'metNames'` | Same semantic field (metabolite display name), different field name — an artifact of cobrapy vs COBRA Toolbox schemas, not a tunable behaviour. |
| `predict_localization.default_compartment` | `'c'` | required arg (no default) | MATLAB has no default at all; Python's `'c'` is a convenience default for the near-universal correct choice, produces no output difference (a MATLAB user must already supply `'c'` explicitly in the common case). Optional, low-priority: MATLAB could add the same default. |
| `run_blast`/`run_diamond`/`run_hmmsearch`/`build_ko_hmm` `threads` | `max(1, cpu_count-1)` | all cores | Confirmed deterministic regardless of thread count — doesn't affect output, so this is a resource-policy choice (leave one core free), not a correctness-relevant value. Both are dynamic ("use available cores") in spirit. |

---

## Summary of required code changes

| Change | File | Priority |
|---|---|---|
| `get_init_model` `allow_excretion` default: `True` → `False` | `src/raven_toolbox/init/build.py` | Medium |
| ~~`get_model_from_homology` `min_align_len` default: `200` → `100`~~ | `src/raven_toolbox/reconstruction/homology/homology.py` | **Done** 2026-08-26 |
| ~~`run_blast`/`run_diamond` `evalue` default: `1e-5` → `1e-4`~~ | `src/raven_toolbox/reconstruction/homology/blast.py` | **Done** 2026-08-26 |
| Docstring: `time_limit` note in `predict_localization` | `src/raven_toolbox/localization/predict.py` | Low |
| Docstring: `mip_gap`/`time_limit` note in INIT functions | `src/raven_toolbox/init/init.py`, `ftinit.py`, `build.py` | Low |

## Changes needed in MATLAB RAVEN for parity

| Change | Priority | Evidence |
|---|---|---|
| KO-assignment cut-off in the `getKEGGModelForOrganism` pipeline (exact sub-function unconfirmed): `1e-50` → `1e-30` | High | `kegg_hmm_cutoff_calibration.md` |
| Same step, gene score ratio: `0.8` → `0.9` | High | Same study |
| `getModelFromHomology` `min_align_len` (confirmed function name, via PR #91's reference to its sibling `getBlast`): `200` → `100` | Medium | `homology_cutoff_calibration.md` |
| Gene-deletion blocked-reaction policy (likely `deleteGenes`, name not confirmed this session): default to removing single-gene reactions rather than keeping them with an empty gene rule | Medium | `manipulation.md` — measured false-positive growth on e_coli_core after essential-gene deletion |
| `FSEOF`: expose the flux noise-floor as a tunable and default it to `1e-6`, not the current implicit `1e-8` | Low | `fseof.md` — measured solver-noise false positives |
| Localization prediction (Python: `predict_localization`; MATLAB name not confirmed this session): consider dropping the 900 s default cap, or confirming it's still needed | Low (tentative) | Wall-clock caps aren't portable across hardware; not yet tested against MATLAB's own solver stack |

## Parameters needing further benchmarks

- Sampling `thinning`/`warmup` autocorrelation on ACHR: yeast-GEM multi-chain analysis pending
- INIT `mip_gap`/`time_limit` genome-scale, across solvers: gated on the GLPK-timeout and HiGHS-cobra upstream issues in `init_solver_benchmark.md`; needs real expression data to distinguish solution quality at different gaps once a second working backend exists

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
