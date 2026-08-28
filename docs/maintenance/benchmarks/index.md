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
| `random_sampling` | `thinning` | `100` | N/A | ⚠ keep value but warn hard: yeast-GEM gives ESS≈12 from 300 samples (single-chain), **and** the median reaction fails Gelman-Rubin R-hat>1.1 between independent chains (see `sampling_convergence_calibration.md`) — genome-scale default-settings output should be treated as unconverged for most reactions, not a minority caveat. Use `n_samples≥2600` or switch sampler for genome scale. |
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
| `run_init` | `mip_gap` | `None` | `0.0004` | ✓ keep `None`; docstring now cites measured genome-scale values, not MATLAB's untested one — see [parity decisions](#cross-toolbox-parity-decisions) |
| `run_init` | `time_limit` | `None` | 5000 ms | ✓ keep `None`; docstring now flags the measured >75 min uncapped risk on hard inputs — see [parity decisions](#cross-toolbox-parity-decisions) |
| `get_init_model` | `allow_excretion` | `False` | `false` | ✓ implemented (was `True`; fixed 2026-06-20, `6f3b57c`) |
| `get_init_model` | `eps` | `1.0` | 1.0 | ✓ keep |
| `run_ftinit` | `series` | `'1+1'` | `'1+1'` | ✓ keep (Gustafsson 2023) |
| `run_ftinit` | `force_on` | `0.1` | 0.1 | ✓ keep |
| `run_ftinit` | `big_m` | `100.0` | 100 | ✓ keep (intentional LP tightener; see `init.md`) |
| `run_ftinit` | `mip_gap` | `None` | `0.0004` | ✓ keep `None`; docstring now cites measured genome-scale values, not MATLAB's untested one — see [parity decisions](#cross-toolbox-parity-decisions) |
| `run_ftinit` | `time_limit` | `None` | 5000 ms | ✓ keep `None`; docstring now flags the measured >75 min uncapped risk on hard inputs — see [parity decisions](#cross-toolbox-parity-decisions) |
| `gene_scores_from_expression` | `factor` | `5.0` | 5 | ✓ keep (RAVEN's own formula; "Wang 2012" attribution checked and could not be confirmed, see `init.md`) |
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
| `predict_localization` | `time_limit` | `None` | `maxTime=15` (minutes; MATLAB's simulated-annealing budget, not a MILP cutoff) | ✓ keep `None` — not a value to unify, see [parity decisions](#cross-toolbox-parity-decisions) |
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
| `get_init_model.allow_excretion` | `False` ✓ done | `False` | **`False`** | ~~Python~~ done (2026-06-20) | High | Zero effect at default `prod_weight`, pure inconsistency; see `manipulation.md` |

### Resolved without unifying — `run_init`/`run_ftinit` `mip_gap`/`time_limit`

These were listed as gated pending a genome-scale test with real expression
data. That test already existed on this branch
([init_param_calibration.md](../../studies/init_param_calibration.md), Human-GEM
HCT116/Hart2015) — it just hadn't been connected to this question. It resolves
both parameters, but not by matching MATLAB's numbers:

- **`mip_gap`:** solve time is flat across the gap at single-step scale (model
  build dominates), so a tight gap is nearly free there; the full genome-scale
  pipeline *does* benefit from loosening to `0.01` (~37% faster, Jaccard 0.995).
  Neither number is MATLAB's `0.0004` specifically, and the right choice depends
  on single-step vs. full-pipeline use, so this stays a documented choice rather
  than a new hardcoded default. **Keep `None`** (Gurobi's own `~1e-4` default is
  at least as tight as the measured-good `0.001`); docstring updated with the
  measured numbers.
- **`time_limit`:** the more consequential finding. MATLAB's `5000 ms` cap is far
  too tight for genome scale (measured solves there took 42-901s+); Python's
  uncapped `None` has a real, measured failure mode — one severely-degraded-input
  case ran **>75 minutes** before being manually killed. The study's own working
  values are `120-600 s/step` (ftINIT) and `400 s` (tINIT). **Keep `None`** as the
  code default (still correct for clean-data / single-step use where it's never
  observed to run away), but the docstring now states the measured working values
  and the >75 min risk explicitly, so users degrading input quality know to set a
  cap rather than discovering the runaway case themselves.

This doesn't touch the separate, still-valid finding in
[`init_solver_benchmark.md`](../../studies/init_solver_benchmark.md) that GLPK
ignores `configuration.timeout` entirely and HiGHS doesn't work with cobra in
this stack — Gurobi remains the only backend genome-scale numbers above were
measured on, and the only one a `time_limit` value does anything on today.

### Keep different — forcing identical values would break correctness

| Parameter | Python | MATLAB | Why |
|---|---|---|---|
| `random_sampling.method` | `'achr'` | `'random_objective'` | Different algorithms, not a flag. ACHR is methodologically preferred (uniform interior sampling vs vertex-biased). Unifying means porting ACHR into MATLAB RAVEN — a real implementation project, tracked separately, not a default change. |
| `random_sampling.replace_max_bound` | `False` | `True` | Applying MATLAB's `True` inside cobrapy/optlang makes the sampler unbounded on standard RAVEN-convention models (measured — see `sampling.md`). Whether MATLAB's own solver path handles `True` safely wasn't tested here; either way this is a solver-stack constraint, not a preference. |
| `random_sampling.loopless_good_reactions` | `True` (loopless FVA) | heuristic (exclude FVA ≥ 999) | Different techniques, not a value. Python's is more correct; MATLAB's heuristic over-excludes reactions that legitimately reach capacity. Porting the proper technique to MATLAB is a bigger project than a default flip. |
| `merge_models.match_by` | `'name'` | `'metNames'` | Same semantic field (metabolite display name), different field name — an artifact of cobrapy vs COBRA Toolbox schemas, not a tunable behaviour. |
| `predict_localization.time_limit` | `None` | `maxTime=15` (minutes) | Confirmed by reading `core/predictLocalization.m` directly (2026-08-28): MATLAB solves the problem with simulated annealing (`maxTime` is its search budget — more time generally improves the heuristic answer, no optimality guarantee either way), while Python solves a deterministic MILP (`time_limit` is a solver cutoff — returns a proven-bounded incumbent). The same number plays a structurally different role in each; not a value to unify. Previously listed as a tentative unify-at-`None` candidate, which assumed both sides were solving the same kind of problem. |
| `predict_localization.default_compartment` | `'c'` | required arg (no default) | MATLAB has no default at all; Python's `'c'` is a convenience default for the near-universal correct choice, produces no output difference (a MATLAB user must already supply `'c'` explicitly in the common case). Optional, low-priority: MATLAB could add the same default. |
| `run_blast`/`run_diamond`/`run_hmmsearch`/`build_ko_hmm` `threads` | `max(1, cpu_count-1)` | all cores | Confirmed deterministic regardless of thread count — doesn't affect output, so this is a resource-policy choice (leave one core free), not a correctness-relevant value. Both are dynamic ("use available cores") in spirit. |

---

## Summary of required code changes

| Change | File | Priority |
|---|---|---|
| ~~`get_init_model` `allow_excretion` default: `True` → `False`~~ | `src/raven_toolbox/init/build.py` | **Done** 2026-06-20 |
| ~~`get_model_from_homology` `min_align_len` default: `200` → `100`~~ | `src/raven_toolbox/reconstruction/homology/homology.py` | **Done** 2026-08-26 |
| ~~`run_blast`/`run_diamond` `evalue` default: `1e-5` → `1e-4`~~ | `src/raven_toolbox/reconstruction/homology/blast.py` | **Done** 2026-08-26 |
| ~~Docstring: `time_limit` note in `predict_localization`~~ | `src/raven_toolbox/localization/predict.py` | **Done** 2026-06-20, `f06faa4` |
| ~~Docstring: `mip_gap`/`time_limit` note in INIT functions~~ | `src/raven_toolbox/init/init.py`, `ftinit.py` | **Done** 2026-08-26 — updated with measured genome-scale values instead of untested MATLAB ones |

## Changes needed in MATLAB RAVEN for parity

All rows below were re-checked directly against `core/*.m` on `origin/develop` in
the RAVEN MATLAB repo on 2026-08-28 (not just inferred from raven-toolbox's own
prior write-ups, which is how the FSEOF and localization rows below turned out
to need correcting from an earlier version of this table).

| Change | Priority | Evidence |
|---|---|---|
| `getKEGGModelForOrganism`'s `cutOff` parameter: `10^-50` → `10^-30` | High | `kegg_hmm_cutoff_calibration.md`; current value confirmed directly in `external/kegg/getKEGGModelForOrganism.m` |
| Same function, `minScoreRatioG`: `0.8` → `0.9` | High | Same study and same confirmation |
| `getModelFromHomology`'s `minLen`: `200` → `100` | Medium | `homology_cutoff_calibration.md`; confirmed directly in `core/getModelFromHomology.m`. A fix already exists as an uncommitted local branch (`fix/homology-minlen`) but is stale relative to current `origin/develop` (file moved `reconstruction/homology/` → `core/`) and hasn't been pushed. |
| `removeGenes`'s `removeBlockedRxns` parameter (default `false`, i.e. keep blocked reactions): default to `true` instead, matching raven-toolbox's `'remove'` | Medium | `manipulation.md` — measured false-positive growth on e_coli_core after essential-gene deletion. Confirmed directly in `core/removeGenes.m`: `removeBlockedRxns` defaults to `false`. |
| `FSEOF`: add an explicit noise floor — currently there is **no tolerance at all** (bare `>0`/`<0`/`>`/`<` float comparison), not merely a looser fixed one as an earlier version of this row said | Medium | `fseof.md` measured `1e-8` catching solver-noise false positives *in raven-toolbox*; confirmed directly in `core/FSEOF.m` that MATLAB has no threshold whatsoever, which was not previously checked and is likely a worse version of the same problem, though not independently measured on the MATLAB side |

Localisation's time-budget parameter was removed from this table: `core/predictLocalization.m`
confirmed MATLAB solves the problem with simulated annealing (`maxTime`, default 15 minutes,
a search budget), not the deterministic MILP raven-toolbox uses (`time_limit`, a solver
cutoff) — different algorithms, not a value to port. See the "Keep different" table above.

## Parameters needing further benchmarks

None remaining as a concrete action item. The one open research question —
sampling between-chain convergence — was investigated as far as it reasonably
goes without a genuinely new approach:

- Sampling `thinning`/`warmup` between-chain convergence (Gelman-Rubin R-hat):
  **investigated, no cheap fix found.** At genome scale (yeast-GEM, default
  ACHR settings) the *median* reaction fails the R-hat>1.1 threshold. Neither
  follow-up tried fixes it cheaply: reallocating the same step budget onto a
  bigger `thinning` makes no measurable difference; `method='chrr'` genuinely
  converges on a small model but its genome-scale cost (a fixed per-chain
  rounding step alone costing ~80 min in a trivial probe) makes it impractical
  as a drop-in fix today. See
  [sampling_convergence_calibration.md](../../studies/sampling_convergence_calibration.md)
  for the full investigation, including a possible future unblock (caching
  CHRR's rounding transform across calls — an engineering change, not a
  parameter default, so out of scope here).

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
