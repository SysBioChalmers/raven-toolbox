# RAVEN parameter defaults and how they were determined

Every RAVEN parameter with a non-trivial default has been through one of three
processes: matched to a **literature value** from the method's original paper,
**measured** empirically on real or synthetic models, or **kept for consistency**
with an existing convention (a widely-used underlying tool's own default, or a
value already established elsewhere in RAVEN). This page gives the short version
— current default and a one-line reason — grouped by capability, covering both
the MATLAB and Python (raven-toolbox) implementations. For the full measurements
and methodology, follow the linked study.

:::{admonition} Where the two implementations still differ
:class: note

RAVEN has existed as MATLAB since 2013 and gained a Python implementation
(raven-toolbox) from 2024; a handful of parameters have only been re-measured on
one side and haven't yet been ported to the other. Neither implementation's
un-measured default was assumed correct going in — where the two disagreed, the
answer was to measure, not to defer to whichever came first. Every row below
states today's actual value(s); where the two genuinely still differ, both are
given, with which side is pending. See
[the master parameter index](../maintenance/benchmarks/index.md#cross-toolbox-parity-decisions)
for the full cross-implementation accounting, including the cases that are
*deliberately* kept different because the two implementations rest on different
solvers or algorithms.
:::

---

## Flux sampling

RAVEN MATLAB: `randomSampling` (`method='achr'|'chrr'|'random_objective'`),
`sampleACHR`, `sampleCHRR`, `sampleMaxVolEllipse`, `sampleChebyshevCenter`.
raven-toolbox: `analysis.random_sampling`, `analysis.find_good_reactions`,
`analysis.max_volume_ellipsoid`. Both implementations expose all three sampling
methods through one entry point.

| Parameter | Default | How determined |
|---|---|---|
| `method` | `'achr'` (Python) / `'random_objective'` (MATLAB) | Literature: ACHR gives near-uniform interior sampling; the random-objective method draws polytope vertices instead. Both algorithms are implemented and selectable on both sides — only which one is the *default* differs. |
| `thinning` | `100` | Upstream convention (cobrapy's own ACHR default, which RAVEN's Python side inherited). **Measured insufficient at genome scale**: yeast-GEM gives only ~12 effective samples from 300 stored (ESS), and independent chains frequently disagree entirely (median Gelman-Rubin R-hat fails the convergence threshold). No cheap fix found — see the [convergence study](../studies/sampling_convergence_calibration.md). |
| `warmup` | `1000` | Upstream convention (cobrapy). |
| `n_objectives` | `2` | Literature: Bordel et al. 2010. |
| `replace_max_bound` | `False` | Measured: applying the alternative (`True`) inside the Python/cobrapy solver stack makes the sampler unbounded on standard RAVEN-convention models (±1000 bounds). Kept different by design between implementations — a solver-stack constraint, not a preference; whether MATLAB's own solver path handles it safely wasn't tested here. |
| `loopless_good_reactions` | `True` | Measured more correct than the alternative fixed-threshold (±999) heuristic, which over-excludes reactions that legitimately reach capacity. The proper technique is implemented on the Python side only so far; porting it is a bigger project than a default flip. |
| `n_samples`, `seed`, `min_flux`, `max_attempts` | `1000`, `None`, `False`, `100` | Same on both sides / no evidence favoured a different value. |
| `find_good_reactions.flux_tol` | `1e-9` | Python-only helper (supports the random-objective method's target selection); no cross-implementation comparison. |
| `max_volume_ellipsoid.maxiter` | `150` | Literature: Zhang & Gao 2003. |
| `max_volume_ellipsoid.tol`, `.reg` | `1e-6`, `1e-8` | Same on both sides. |

**Full detail:** [sampling.md](../maintenance/benchmarks/sampling.md) ·
[sampling convergence study](../studies/sampling_convergence_calibration.md) ·
[CHRR/ACHR algorithm reference](flux_sampling_algorithms.md)

---

## FSEOF

RAVEN MATLAB: `FSEOF`. raven-toolbox: `analysis.fseof`.

| Parameter | Default | How determined |
|---|---|---|
| `n_steps` | `10` | Literature: Choi et al. 2010. |
| `max_fraction` | `0.9` | Literature: Choi et al. 2010; measured — values below 0.9 pick up spurious targets on iJO1366. |
| `correlation_threshold` | `0.9` | Literature: Choi et al. 2010; measured — 0.7 adds 4 spurious amplification targets on iJO1366. |
| `flux_eps` | `1e-6` on the Python side; MATLAB has no threshold at all | Measured: an explicit `1e-6` floor avoids classifying accumulated solver noise as a false-positive knockdown target on genome-scale models. MATLAB's target classification is a bare float comparison with no tolerance whatsoever, stricter than any fixed floor — real solver noise there is expected to produce at least as many spurious classifications as measured here, likely more. |

**Full detail:** [fseof.md](../maintenance/benchmarks/fseof.md)

---

## INIT / ftINIT

RAVEN MATLAB: `runINIT`, `scoreComplexModel`, `getINITModel`, `ftINIT`, `prepINITModel`,
`ftINITInternalAlg`, `getINITSteps`. raven-toolbox: `init.run_init`,
`init.score_reactions_from_genes`, `init.gene_scores_from_expression`,
`init.get_init_model`, `init.ftinit`, `init.prep_init_model`, `init.run_ftinit`,
`init.get_init_steps`, `init.classify_reactions`.

| Parameter | Default | How determined |
|---|---|---|
| `prod_weight` | `0.5` | Literature: Ågren et al. 2012 (the original INIT paper). |
| `allow_excretion` | `False` (all three entry points) | Measured zero effect at the default `prod_weight` across all three entry points. |
| `mip_gap` | `None` on the Python side (solver default); `0.0004` on the MATLAB side | Measured on genome-scale Human-GEM: the solver's own gap (~1e-4 on Gurobi) is already at least as tight as the measured-good value, so Python's `None` needs no change; genome-scale users get concrete numbers (`0.01` for the full ftINIT pipeline) documented rather than a single hardcoded default, since the right value depends on single-step vs. full-pipeline use. |
| `time_limit` | `None` (uncapped) on the Python side; `5000 ms` on the MATLAB side | The same genome-scale study found the MATLAB value far too tight (real solves took 42–901s+) but also a real >75-minute uncapped runaway case on degraded input on the Python side — `None` was kept as the safer default there, with the measured working range (120–600s/step) documented for users pushing hard or noisy input. |
| `series` | `'1+1'` | Literature: Gustafsson et al. 2023 (the ftINIT paper). |
| `force_on` | `0.1` | RAVEN's original value; measured near-insensitive across a 0.02–0.5 range. |
| `big_m` | `100.0` | RAVEN's original value; confirmed as an intentional LP-relaxation tightener (not a flux cap) — see the linked study for why. |
| `eps` | `1.0` | Same on both sides. |
| `gene_scores_from_expression.factor`, `.max_score`, `.min_score` | `5.0`, `10.0`, `-5.0` | Literature: Wang et al. 2012. |
| `score_reactions_from_genes.*`, `classify_reactions.*` | various | Same on both sides / standard practice. |

**Full detail:** [init.md](../maintenance/benchmarks/init.md) ·
[INIT parameter calibration study](../studies/init_param_calibration.md) ·
[INIT solver benchmark](../studies/init_solver_benchmark.md)

---

## Gap-filling

RAVEN MATLAB: `fillGaps` (connectivity mode) — ported as `gapfilling.connect_blocked_reactions`.
`gapfilling.fill_gaps_fast_lp` and `gapfilling.fill_gaps_kumar_milp` are literature-based
methods added to the Python implementation, with no MATLAB RAVEN counterpart.

| Parameter | Default | How determined |
|---|---|---|
| `fill_gaps_fast_lp.epsilon` | `0.0001` | Literature: Thiele et al. 2014 (the fastGapFill paper). |
| `connect_blocked_reactions.eps` | `1.0` | Measured on a synthetic supply-limited model; reliable noise filter for RAVEN-convention (±1000-bound) models — lower it for tightly-constrained exchanges. |
| `connect_blocked_reactions.penalty`, `.allow_net_production` | `1.0`, `False` | Measured / standard practice for reconstruction (vs. curation) use. |
| `fill_gaps_kumar_milp.weights` | `(1.0, 2.0)` | Literature: Kumar et al. 2007; measured — confirms reversal is preferred over adding a new reaction at the intended 2:1 ratio. |
| `fill_gaps_kumar_milp.big_m` | `1000.0` | Measured: matches RAVEN's ±1000 bound convention; a smaller value artificially caps reversed-reaction flux. |

**Full detail:** [gapfilling.md](../maintenance/benchmarks/gapfilling.md)

---

## Homology-based reconstruction

RAVEN MATLAB: `getModelFromHomology`, `getBlast`, `getDiamond`. raven-toolbox:
`reconstruction.homology.get_model_from_homology`, `.run_blast`, `.run_diamond`.

| Parameter | Default | How determined |
|---|---|---|
| `evalue` (both aligners) | `1e-4` | Matches the MATLAB value; measured to have no effect on the final model either way, since `get_model_from_homology`'s much stricter filters dominate downstream. |
| `threads` | `max(1, cpu_count-1)` (Python); dynamic "all cores" (MATLAB) | Measured deterministic across thread counts — pure performance change, ~1.9–4× speedup. Not treated as a "value" that needs unifying, since both sides already mean "use the available cores". |
| `max_evalue` | `1e-30` | Measured **inert** from `1e-4` to `1e-50` — identity and alignment-length have already excluded whatever a looser e-value would admit. |
| `min_align_len` | `100` on the Python side; MATLAB still uses `200` | Measured against independent KEGG and OMA ortholog references across a 4-organism relatedness series: recovers 3–4 points of recall over `200` for ≤0.6 points of precision cost. A MATLAB back-port is proposed but not yet done. |
| `min_identity` | `40` | Same study — the binding filter and the measured optimum, weighting a wrongly-transferred reaction as worse than a missing one. |
| `bidirectional` | `True` | Same on both sides. |

An earlier arm of the same study tried scoring against curated GEMs (hanpo-GEM,
rhto-GEM) built by this same reconstruction method, and was retired once it showed
the obvious circularity directly: agreement with the curated model peaked exactly
at that model's own build settings and collapsed one step past them.

**Full detail:** [reconstruction_homology.md](../maintenance/benchmarks/reconstruction_homology.md) ·
[homology cut-off calibration study](../studies/homology_cutoff_calibration.md)

---

## KEGG-based reconstruction

RAVEN MATLAB: `getKEGGModelForOrganism` (the KO-assignment sub-step's exact internal
name hasn't been independently confirmed). raven-toolbox: `reconstruction.kegg.assign_kos`,
`.run_hmmsearch`, `.build_ko_hmm`, and the `get_kegg_model_*` family.

| Parameter | Default | How determined |
|---|---|---|
| `assign_kos.cutoff` | `1e-30` on the Python side; MATLAB still uses `1e-50` | Measured against real KEGG gene→KO annotations across 4 organisms of varying study depth: the tighter value sits inside the tail of *real* matches (~20 orders of magnitude above the noise floor), discarding genuine hits for no gain. MATLAB's port is still pending. |
| `assign_kos.min_score_ratio_g` | `0.9` on the Python side; MATLAB still uses `0.8` | Same study — the real precision lever (0.8→0.95 lifts precision ~0.07–0.10 for ~0.02 recall). MATLAB's port is still pending. |
| `assign_kos.min_score_ratio_ko` | `0.3` | Same study — measured empirically **inert** (varying it changes precision/recall by ≤0.02); kept at the original value since it does nothing either way. |
| `build_ko_hmm.seq_identity` | `0.9` | Literature: CD-HIT's own recommendation for protein clustering. |
| `run_hmmsearch.threads`, `build_ko_hmm.threads` | `max(1, cpu_count-1)` | Performance-only; HMMER is deterministic (to within insignificant float noise) across thread counts. |
| `get_kegg_model_*.keep_spontaneous`, `.keep_undefined_stoich`, `.keep_incomplete`, `.keep_general` | `True`, `True`, `True`, `False` | Same on both sides / standard draft-reconstruction practice. |

**Full detail:** [reconstruction_kegg.md](../maintenance/benchmarks/reconstruction_kegg.md) ·
[KEGG HMM cut-off calibration study](../studies/kegg_hmm_cutoff_calibration.md)

---

## Sub-cellular localisation

RAVEN MATLAB: `predictLocalization`. raven-toolbox: `localization.predict_localization`.
**These solve the problem with genuinely different algorithms**: MATLAB uses
simulated annealing (a stochastic heuristic, no optimality guarantee); the Python
side solves a deterministic MILP. A few parameters below are therefore not really
a "which value is right" question — the same number means something different in
a heuristic's search budget than in a solver's optimality gap.

| Parameter | Default | How determined |
|---|---|---|
| `default_compartment` | `'c'` on the Python side; MATLAB requires the argument every call | Convenience default — cytosol is correct for the vast majority of reactions. Not a value disagreement, since MATLAB simply has no default to compare against. |
| `transport_cost` | `0.5` | Not independently benchmarked; same on both sides. Calibration note: for real (non-integer-scale) predictor scores, this typically needs dialing down — see the linked study. |
| `time_limit` (Python) / `maxTime` (MATLAB) | `None` (uncapped) on the Python side; MATLAB's `maxTime` defaults to `15` minutes | Not a value to unify: MATLAB's number is a simulated-annealing search budget (more time generally means a better heuristic answer, not a proof of optimality), while Python's is a MILP solver cutoff (returns a proven-bounded incumbent). Python's `None` was validated at ~2.5 minutes on yeast-GEM, the primary development-scale model. Whether the MILP needs a cap at all for harder cases hasn't been stress-tested. |
| `mip_gap` | `None` | Python-side-only — MATLAB's heuristic has no analogous optimality-gap concept. |

**Full detail:** [localization.md](../maintenance/benchmarks/localization.md) ·
[yeast-GEM localisation benchmark](../studies/yeast_localization_benchmark.md)

---

## Model manipulation

RAVEN MATLAB: `removeGenes`, `mergeModels`, `addRxns`, `addTransport`, `simplifyModel`
(reversibility-constraining and duplicate-removal modes). raven-toolbox:
`manipulation.remove_genes`, `.merge_models`, `.add_reactions_from_equations`,
`.add_transport_reactions`, `.constrain_reversible_reactions`, `.find_duplicate_reactions`.

| Parameter | Default | How determined |
|---|---|---|
| `remove_genes.blocked_reactions` | `'remove'` on the Python side; MATLAB's equivalent policy keeps the reaction | Measured on e_coli_core: keeping the reaction gives false-positive growth after an essential gene is removed (it survives with an empty gene rule). `'remove'` is correct for gene-essentiality / metabolic-engineering use; the keep-style policy remains available on the Python side too, for annotation-curation workflows. MATLAB's port of this fix is still pending. |
| `constrain_reversible_reactions.eps` | `1e-9` | Matches Gurobi's own default primal feasibility tolerance. |
| `merge_models.match_by` | `'name'` (Python) / `'metNames'` (MATLAB) | Same semantic field (metabolite display name), different schema field name — not a value disagreement to resolve. |
| `remove_genes.remove_orphans`, `add_reactions_from_equations.mets_by`, `find_duplicate_reactions.ignore_direction`, `add_transport_reactions.reversible`/`.only_to_existing` | various | Python-side additions/refinements; sensible defaults, no direct MATLAB equivalent to compare against. |

**Full detail:** [manipulation.md](../maintenance/benchmarks/manipulation.md)

---

## Tasks

RAVEN MATLAB: `checkTasks`. raven-toolbox: `tasks.check_tasks`, `tasks.find_task_essential_reactions`.

| Parameter | Default | How determined |
|---|---|---|
| `close_boundaries` (both functions) | `True` | Matches RAVEN's implied behaviour. |
| `find_task_essential_reactions.tol` | `1e-8` | Python-side-only parameter. |

**Full detail:** [tasks.md](../maintenance/benchmarks/tasks.md)

---

## Still open

- **Genome-scale ACHR sampling convergence** has no cheap fix (see the sampling section
  above and its linked study in full).
- Four values are measured and live on the Python side only, with MATLAB pending: the
  KEGG KO-assignment cut-off and score ratio in `getKEGGModelForOrganism`, homology's
  `min_align_len` in `getModelFromHomology`, the `removeGenes` gene-deletion policy, and
  FSEOF's noise floor (which needs adding to MATLAB rather than changing, since none
  exists there today). Tracked in
  [the master parameter index](../maintenance/benchmarks/index.md#changes-needed-in-matlab-raven-for-parity).
