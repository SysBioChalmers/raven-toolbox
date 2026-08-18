# Changelog

Milestones in the raven-toolbox port. For function-level status see
[docs/raven_migration.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/reference/migration.md); for open work see
[docs/todo.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/reference/todo.md).

## Unreleased

* **`assign_compartments` gap-fill: a reliable flux-based fill instead of cobra's MILP.** The
  universal-DB gap-fill in `localization.certify` no longer calls `cobra.flux_analysis.gapfill`, whose
  indicator MILP, at genome scale, fails to find a valid fill in the **majority** of cases *even when the
  exact reaction that restores growth is present in the universal* (it returns an incumbent its own
  validation rejects, and raises rather than offering a fill). The replacement adds the candidates on a
  working copy, holds biomass at the floor, runs pFBA, and returns the flux-carrying additions — **sorted**,
  so the result does not depend on which co-optimal vertex the solver picked. A plain LP cannot have the
  MILP's failure mode: a returned set is a real flux solution that reaches the floor. On single-reaction
  knockout-recovery (remove an essential reaction, fill from a universal that contains it), recall is
  **60/60 — the exact reaction each time — vs cobra's 45 %**; on realistic incomplete drafts (12% of
  internal reactions dropped, many simultaneous gaps) it restores growth **12/12 vs cobra's 0/12**. The universal must share the draft's
  metabolite namespace: candidates are matched by id, as cobra's gapfill required, and a mismatch now
  **warns** instead of returning `[]` silently. The set is flux-parsimonious (pFBA) but not guaranteed
  reaction-count-minimal; the caller re-certifies with a real FBA, so no false certificate is possible.
* **Deterministic, score-aligned reaction placement in `assign_compartments`.** The placement master
  maximises a gene-localisation objective that never mentions the per-reaction placement variable, so
  each reaction's compartment was a free co-optimum: the pinned solver returned it reproducibly but
  *arbitrarily*, and reaction agreement with curated yeast-GEM was only **52.8 %** (an earlier un-pinned
  build happened to land ~72 %). A lexicographic second pass now fixes the gene layout to the primary
  optimum, then places each reaction in the compartment its own enzymes are predicted to occupy — the
  summed DeepLoc score of the reaction's genes, with a small `default_compartment` prior so genes-free
  and score-tied reactions fall there deterministically. Yeast reaction agreement rises to **72.5 %**
  (1408/1943) and now rests on the localisation evidence rather than a solver tie-break; **gene agreement
  is unchanged** (88.7 %, 716/807 — the gene layout is untouched); coherent placement fragments fewer
  metabolites, so it adds **fewer transports** (1001 → 967); growth and blocked-fraction are unchanged;
  the second solve is warm-started and adds ~4 s. Reproducible run-to-run (identical placement across
  independent runs).
* **`diff_models` compares grRules as logic, not text.** The GPR check now DNF-expands each rule (via the
  existing `manipulation.gpr_to_dnf`), sorts the genes within each isozyme clause and sorts the clauses, so
  operand order no longer registers as a difference: `a and b` == `b and a` and `a or b` == `b or a`. The
  previous heuristic only lowercased and collapsed whitespace, so it flagged logically identical rules that
  differed only in operand order. This brings `diff_models` in line with MATLAB RAVEN's `diffModels`
  ([RAVEN #686](https://github.com/SysBioChalmers/RAVEN/pull/686)); a rule cobra cannot parse falls back to
  the old string comparison, so malformed rules are still compared rather than silently equated.
* **Fix `load_delta_g_csv` recording the ΔG side-car tables' "missing" sentinel as a real measurement.** The
  side-car tables encode "no valid ΔG" as `10000000`, and the loader — written for exactly these files, down to its
  `Var1`/`Var2` defaults — stamped it verbatim, presenting a physically impossible 10⁷ kJ/mol as a
  measurement on **777 of yeast-GEM's 4102 reactions (19.5%)**. yeast-GEM's own `checkrxnDirection.m` gates
  on the same value (`if ~isequal(seed_rxnInfo{...},'10000000.0') %check if database contains valid deltaG
  value`). The sentinel is now treated as missing, as NaN already was, recognised whichever dtype the CSV
  round-trip produces; the new keyword-only `missing_value` (default `DELTA_G_MISSING`) tunes or
  disables it. Real ΔG coverage of yeast-GEM is 78.2%, not the 97.1% the loader previously implied.
* **Wire the confidence facets together.** `confidence.annotate_confidence(model, proposal=..., scores=...)`
  runs every applicable scorer in one call and returns `{facet: reactions_scored}` — `equation` and
  `gene_association` need only the model, `localization` runs only when a proposal and its scores are given
  (skipped, not failed, otherwise). `curation_priority` now drops a placement a curator has settled with
  `mark_curated` from the review queue (new `include_curated=False`), so a settled reaction stops
  resurfacing; `include_curated=True` keeps it. The no-SBO-terms warning now names its remedy,
  `raven_toolbox.annotation.add_sbo_terms(model)`.

## 0.3.0 — 2026-07-16

Compartment localisation and per-reaction confidence tracking, new gap-filling and flux-sampling
algorithms, ftINIT parity with RAVEN, and KEGG artefact hosting moved to the dedicated `raven-data`
repository.

* **`export_for_git` can pin the MATLAB `.mat` variable name.** New keyword-only `varname`, forwarded
  to `save_matlab_model`, so a repository that expects a specific struct name (Human-GEM expects
  `humanGEM`) no longer has to write the `.mat` itself. Default `None` keeps cobra's own fallback (the
  model id), so existing callers are unaffected.
* **ftINIT now matches RAVEN's model preparation, solver parameters, and gap schedule.** A line-by-line
  review against RAVEN's `ftINIT` / `prepINITModel` found the formulation faithfully ported, with the
  real divergences concentrated in model simplification and solver/gap handling. `prep_init_model` runs
  both RAVEN simplifications in `prepINITModel` order and broadens the exchange mask to RAVEN
  `getExchangeRxns`' one-sided rule (any reaction with no products or no substrates) — previously only
  the topological dead-end pass ran, leaving a ~10 % larger model. `ftinit` adopts RAVEN `optimizeProb`
  solver parameters (`Threads=1` for a deterministic MILP incumbent, `Presolve=2`, `1e-9` feasibility /
  optimality / integrality tolerances), nudges tiny reaction scores off zero, forces permanent essential
  reactions at `min(0.99*|carried flux|, force_on)` while tracking carried flux across steps, and
  escalates the MIP gap per step (RAVEN `MILPParams` / `AbsMIPGaps`). New `manipulation.simplify_model`
  mirrors RAVEN's `simplifyModel` boolean-flag interface, adding `remove_zero_interval_reactions`
  (`deleteZeroInterval`) and `remove_no_flux_reactions` (FVA-blocked removal, `deleteMinMax`).
  Task-essential discovery follows RAVEN `checkTasks`, and task gap-fill follows `ftINITFillGaps`. On
  Human-GEM / DLD1 the prep model now matches RAVEN's sizes (`ref_model` 10240 against RAVEN's ~10198;
  merged `min_model` 6959 against ~6917 — previously 11532 / 8252). `allow_excretion` deliberately keeps
  `S*v >= 0`, matching the flag's name and classic INIT rather than RAVEN's `csense 'L'`; it is unused in
  the default `'1+1'` schedule and documented inline.
* **Structural confidence facets: `equation` and `gene_association`.** `confidence` gains
  `score_equation_confidence` (mass & charge balance, formula completeness) and
  `score_gene_association_confidence` (GPR presence + literature corroboration), plus the public
  `equation_exempt` / `gene_association_exempt` predicates and `facet_summary`. Two rules govern every
  score, both forced by `overall = min(facets)`: a facet that does not apply to a reaction is **not
  written at all** (an exchange reaction is imbalanced by construction, so writing `1.0` would make 469
  never-checked yeast-GEM reactions indistinguishable from the 3617 verified ones), and `0.0` always means
  *evidence contradicts the model* — never *evidence is missing*, which keeps `overall == 0.0` a usable
  filter. Exemptions are SBO-driven, never name-driven: a `\bgrowth\b` regex would silence the chemistry
  check on "non-growth associated maintenance reaction". Validated on yeast-GEM, where the gene rubric
  recovers the model's own curator-assigned `Confidence Level` (99.9% / 91.5% / 95.3% per band) without
  ever reading it. `mark_curated` now takes a `facet`, and `clear_confidence` is exported
  ([design](docs/studies/confidence_tracking.md)).
* **Fix `get_elemental_balance` crashing on an unparseable formula.** A parenthesised polymer such as
  `(C5H8)n` (glycogen, starch) makes cobra's `Metabolite.elements` return `None`, which
  `check_mass_balance()` turns into a `ValueError` — so the shipped helper raised part-way through a model
  instead of reporting it. A present-but-uninterpretable formula is now `unknown`, which is what the
  function's own contract already promised.
* **Fix `score_localization_confidence` vetoing reactions it could not measure.** A reaction whose genes are
  absent from the score table was written `0.0`/`connectivity`, which under `overall = min(facets)` vetoed
  the whole reaction on the strength of a *missing input* rather than of evidence. It now abstains; a gene
  that is scored but scores zero at the assigned compartment still earns a real `0.0`.
* **Fix gap-fill materialisation landing on the wrong compartment metabolite.**
  `apply_assignment`'s `_add_universal_reaction` matched a universal candidate's metabolite ids
  verbatim against the target model instead of resolving them through the same base-id/compartment
  lookup `_move_reaction` already uses for relocated reactions. Wherever a draft's non-default-
  compartment species only exist because a relocated reaction created them (any id scheme other
  than the universal database's own), the gap-fill reaction silently materialised as a disconnected
  island: the assignment's own solved objective already accounted for the correct shared node, only
  the applied model was wrong. Fixed with a regression test
  (`test_gapfill_reuses_relocated_compartment_metabolite`).
* **Evidence-aware transport scoring (first increment).** New `localization.transport_evidence` turns
  per-gene transporter evidence into the per-metabolite `transport_cost` mapping the assignment MILPs
  already accept, so a transport is cheap when a transporter gene supports it (right substrate, right
  membrane) and pays the full prior otherwise: `evidence_aware_transport_cost`
  (`cost = base·(1−evidence)`), `annotate_transporters` (bring-your-own annotation table), and
  `TransporterAnnotation`. Carrier-general and organism-agnostic; the `hmmsearch` (Pfam) / `diamond`
  (TCDB) annotation back-ends are the next increment
  ([plan](docs/reference/transport_evidence_scoring.md)).
* **Consolidated `assign_compartments` into raven-toolbox.** The functionality-constrained
  compartment-assignment MILP — biomass/growth floor + flux gating + optional gap-fill + sound
  reaction-level multi-localisation — moves from the standalone `edkerk/assignCompartments` repo into
  `localization/` (as `assign_compartments`/`apply_assignment`/`AssignmentProposal`), coexisting with
  the score-driven `predict_localization` ([design](docs/reference/multi_localization_design.md)). The
  MATLAB port is tracked in [MATLAB back-ports](docs/reference/matlab_raven_backports.md).
* **CarveFungi assignment head-to-head on its own MILP.** Ran CarveFungi's *own* `minmax_reduction`
  carve-MILP (CPLEX, unmodified) with our transport-minimising term swapped into its objective, on its
  real universal-DB candidate set + DeepLoc-injected scores
  ([study](docs/studies/carvefungi_milp_benchmark.md), `scripts/run_carvefungi_cplex.py`;
  `scripts/benchmark_carvefungi_milp.py` is a Gurobi re-implementation for the formulation study).
  Adversarially verified (which caught and fixed a compartment-id parsing bug). Findings: adding our
  transport cost yields ~1.6× fewer inter-compartment transports per reaction (41% fewer) at no
  detectable assignment-accuracy cost (86.0% vs 85.5% recall; 93% identical placements). The carve's
  MILP formulation is hard — neither CPLEX nor a tighter Gurobi port proves optimality — so these are
  deterministic, time-budget-stable near-optimal incumbents, reported with their gaps. Also surfaced
  that CarveFungi's *shipped* yeast localisation file is inert (RefSeq vs ORF id mismatch).
* **Head-to-head vs RAVEN predictLocalization + CarveFungi positioning.** Benchmarked the
  deterministic compartment-assignment MILP against RAVEN's stochastic `predictLocalization` on
  identical yeast-GEM + DeepLoc inputs ([study](docs/studies/predictlocalization_comparison.md),
  `scripts/compare_predictlocalization.py` + `scripts/run_predictlocalization.m`): on the common gene
  set the MILP is ~7 pp more accurate (83.9% vs 76.8%), deterministic (vs 35% of genes flipping
  between SA runs), and faster (90 s vs a multi-minute budget). Also a source-level
  [analysis of CarveFungi](docs/studies/carvefungi_analysis.md), the contemporary
  carve-a-universal-model method, showing how our transport-minimising assignment differs.
* **Curation triage for localisation.** Added `triage_localization` — an optional companion to
  compartment assignment that ranks the genes/reactions whose localisation is shakiest (low DeepLoc
  confidence, borderline top-two margin, multi-source disagreement, no evidence, low-trust
  compartment, multi-localised), each with a plain-English reason, so a curator knows where to look.
  Returns a `ReviewReport`. `load_deeploc` gained `keep_raw_confidence=True` and `LocalizationScores`
  a `raw_confidence` field (per-gene normalisation otherwise discards the confidence the triage needs).
* **Finetuned localisation hyperparameters on the slow yeast run.** Refreshed the triage
  `DEEPLOC_COMPARTMENT_TRUST` table from the slow (ProtT5) data (mitochondrion 0.67 → 0.86, `mm` now
  trustworthy via the validated split, Golgi 0.23 → 0.01) and re-validated the `min_confidence` gate
  (0.7 → 88.3% corroboration, 80% kept) and `membrane_threshold` (0.50 is inside the optimal plateau)
  in a new [finetuning study](docs/studies/localization_finetuning.md)
  (`scripts/finetune_localization_yeast.py`).
* **Cross-species DeepLoc benchmark.** Generalised the predictor benchmark to any curated model
  (`scripts/benchmark_deeploc.py --species {yeast,aracore,icre1355}`, a per-species compartment
  config) and added independent non-yeast eukaryotes. DeepLoc 2.1 generalises across kingdoms — the
  chloroplast is recovered in both *Arabidopsis*
  ([AraCore](docs/studies/deeploc_aracore_benchmark.md), **80.3%** overall, plastid 89.9%) and the
  green alga *Chlamydomonas* ([iCre1355](docs/studies/deeploc_icre1355_benchmark.md), plastid 78%,
  though algal cytosol/mito are poor on an auto-generated model). A gene-level
  [Human-GEM control](docs/studies/deeploc_humangem_benchmark.md)
  (`scripts/benchmark_deeploc_humangem.py`) reaches 84.7% but, crucially, **excludes the 439 (15%)
  gene compartments Human-GEM sourced from DeepLoc2** (which score 93.8% — DeepLoc grading itself).
  The yeast run was refreshed to DeepLoc's slow (ProtT5) model, lifting organelle-collapsed accuracy
  54.6% → 64.6% and mitochondrial-membrane recall 47% → 86%.
* **Optional raw DeepLoc probabilities.** `load_deeploc` gained `normalise=False` to keep DeepLoc's
  calibrated probabilities instead of rescaling each gene's best compartment to 1.0. A whole-model
  yeast-GEM benchmark ([study](docs/studies/deeploc_normalisation_benchmark.md)) finds normalisation
  is **accuracy-neutral** for compartment assignment (raw does not rescue the contested or
  high-confidence calls); the only reproducible difference is that raw assigns fewer genes to
  multiple compartments — a re-scaling of the existing `transport_cost`/`multi_compartment_penalty`
  knobs, not new signal. So normalisation stays the **default** and `normalise=False` is an opt-in
  for callers wanting the calibrated magnitudes (e.g. the `triage_localization` confidence signal).
* **Fuse and tune localisation evidence.** Added `combine_scores` (weighted-sum consensus of several
  `LocalizationScores`, so agreement across DeepLoc / UniProt / COMPARTMENTS is reinforced), and gave
  `load_deeploc` / `load_mulocdeep` a `min_confidence=` gate (drop unreliable low-confidence genes)
  plus, for `load_deeploc`, `membrane_split={"m":"mm"}` (route mitochondrion to its membrane
  sub-compartment using the transmembrane signal — mito only; ER is not separable). Motivated and
  validated by the [DeepLoc 2.1 yeast-GEM benchmark](docs/studies/deeploc_yeast_benchmark.md).
* **Prepare sequence-predictor input.** Added `prepare_deeploc_input` (plus `fetch_protein_sequences`
  and `write_fasta`) to write a DeepLoc-2.1-ready protein FASTA for a model's genes — sequences
  fetched from UniProtKB, headers set to the gene ids so the predictor output lines up with the model
  and `load_deeploc`. DeepLoc 2.1 has no batch API; the FASTA is chunked at the web server's
  500-sequence limit, and genes without a reviewed sequence are reported. Script:
  `scripts/prepare_deeploc_yeast.py`.
* **Localisation loaders modernised.** Added `load_mulocdeep` (MULocDeep wide tables),
  `load_compartments` (the COMPARTMENTS evidence database), `load_uniprot` (curated UniProtKB
  `Subcellular location` exports) and `fetch_uniprot_localization` (the same via the UniProt REST
  API by organism id), plus `DEFAULT_COMPARTMENT_MAP` to rename predictor labels to
  model compartment ids and collapse synonyms. `load_deeploc` gained a `compartment_map`
  argument. **Removed `load_wolfpsort`** — modern multi-label predictors, the COMPARTMENTS
  database and UniProt supersede the single-label WoLF PSORT caller.
* **Flux sampling: CHRR and ACHR, unified under `random_sampling`.** `random_sampling(model,
  method=...)` is the single entry point and dispatches `"achr"` (the new default), `"chrr"`, and
  `"random_objective"` (the historical Bordel et al. 2010 vertex method). **CHRR** — Coordinate
  Hit-and-Run with Rounding (Haraldsdóttir et al. 2017) — does nullspace reduction, maximum-volume
  ellipsoid rounding, then coordinate hit-and-run; it is the recommended sampler for enzyme-constrained
  (ecModel + proteomics) and flux-measured models, whose feasible set is a thin, ill-conditioned slab
  that defeats unrounded chains. cobrapy ships ACHR but no CHRR, so the ACHR path wraps
  `cobra.sampling.ACHRSampler` while CHRR is a genuine new implementation. Also exports
  `max_volume_ellipsoid` (Zhang & Gao 2003 primal-dual interior-point MVE solver), validated against
  analytic cases (box → unit ball, scaled/sheared box, triangle → Steiner inellipse). All methods return
  a unified `FluxSamplingResult`; `RandomSamplingResult` is kept as an alias. Reference:
  `docs/reference/flux_sampling_algorithms.md`. **Breaking:** `random_sampling`'s default changed from
  the random-objective vertex method to ACHR — pass `method="random_objective"` for the previous
  behaviour.
* **Gap-filling algorithms: LP, MILP, and topological.** Three strategies in `raven_toolbox.gapfilling`
  complementing `connect_blocked_reactions`: `fill_gaps_fast_lp` (LP-relaxation connectivity
  gap-filling, fastGapFill / Thiele et al. 2014, with `variant="swift"` for the SWIFTCORE single-LP
  form, Tefagh & Boyd 2020 — no MILP solve); `fill_gaps_kumar_milp` (Kumar et al. 2007 global
  growth-floor MILP, adding directionality-reversal repair on top of database-reaction addition); and
  `analyse_topology` (Meneco-inspired BFS metabolite-producibility scope, reporting unreachable
  metabolites and pruning candidate reactions with no solver call). These complement
  `cobra.flux_analysis.gapfill`, which does objective-based gap-filling without reversal repair.
  References: `docs/reference/gap_filling_algorithms.md`, `docs/reference/cobra_raven_comparison.md`.
* **KEGG artefacts hosted in `raven-data`.** All KEGG artefact URLs move from `raven-toolbox` releases
  to the dedicated `raven-data` repository (`raven-data/releases/download/kegg118/`), keeping the
  toolbox release lean. Adds `scripts/publish_to_raven_data.py` to upload build artefacts to a
  `raven-data` release, and updates the `_DATA_REGISTRY` URLs in `data.py` to match. Maintenance docs
  (`maintaining_binaries.md`, `maintaining_kegg_data.md`, `data_manifest.md`) describe the 3-step
  publish workflow; `docs/reference/matlab_raven_backports.md` is slimmed to the deliberate-omission
  section now that every backport item is complete.
* **kegg118 artefact set.** Regenerates `data/manifest.json` for kegg118 (versions, SHA256 checksums,
  byte sizes for the core bundle, taxonomy, and the prokaryote/eukaryote HMM libraries) and syncs the
  in-code `_DATA_REGISTRY`, which had been left at kegg116 and is the fallback when
  `$RAVEN_PYTHON_MANIFEST` is unset. Fixes the generated release-asset URLs, which pointed at the
  singular `…/release/download/…` path that GitHub does not serve — every download would have 404'd —
  and `scripts/make_registry_snippet.py` now takes a release `--tag` and builds the
  `…/releases/download/<tag>` prefix itself so the typo cannot recur. The HMM libraries are named
  `kegg<version>_prokaryotes` / `_eukaryotes` consistently across the manifest, the Python resolver, and
  MATLAB RAVEN.
* **Binary provisioning: sets, fetch CLI, auto-fetch toggle, native-Windows HMMER.** Provisioning is
  decoupled from pip (extras can only pull PyPI wheels, and downloading binaries during `pip install`
  is a known anti-pattern) into: **binary sets** — `runtime` (blast, diamond, hmmsearch) vs `build`
  (hmmbuild, mafft, cd-hit); an explicit **fetch CLI**, `raven-toolbox-binaries --set runtime|build|all`
  (and `--list`), which is OS-aware, SHA256-verified, skips tools already on PATH, and reports per-tool
  `present/downloaded/unavailable/error`; and the unchanged **lazy first-use download**, now disableable
  via `RAVEN_PYTHON_AUTOFETCH=0` for air-gapped or conda-managed setups. Adds a `windows-x86_64` entry
  on the `hmmer` bundle (HMMER 3.3.2 repackaged from RAVEN 2.10.5) so the KEGG HMM *query* runs on
  native Windows without WSL — searching 3.4-built libraries with 3.3.2 is safe because the toolbox
  ships ASCII `.hmm` libraries, the `HMMER3/f` format is unchanged 3.1→3.4, and 3.4 introduced no
  protein-scoring change. HMM *building* stays WSL/conda-only (MAFFT and CD-HIT have no Windows
  builds). `ensure_binary` now resolves `<name>.exe` on Windows, and sibling DLLs extract next to it.
* **Resumable KEGG artefact build.** `scripts/build_kegg_artefacts.py` skips each stage whose output
  already exists (parsed tables, taxonomy, per-domain HMM library, core bundle), so a build that dies
  partway continues on re-run of the same command instead of restarting the multi-hour HMM step; a new
  `--force` rebuilds from scratch. `build_ko_fastas` keeps already-written `<KO>.fa` files and writes a
  `.ko_fastas_complete` marker so a finished run fast-paths without re-scanning `genes.pep`.
  `organism_gene_ko` is loaded lazily, so a fully-published build re-runs as a clean no-op.
* **Progress reporting for the KEGG build.** An opt-in `progress` flag (tqdm) on the maintainer-side
  KEGG path, so long-running steps report progress instead of appearing to hang: per-file byte bars for
  download and the multi-GB proteome gunzip, a byte bar over the dominant `organism_gene_ko` streaming
  pass in `parse_kegg_dump`, and an "N of M KOs" counter over the HMM build loop. `progress` is
  independent of `verbose`; with both on, per-KO log lines route through `tqdm.write` so they do not
  corrupt the bar.
* **Excel export of `model.ec`.** `export_to_excel` writes a populated `model.ec` (`EcData`) to two
  further sheets, mirroring RAVEN's `exportToExcelFormat`: **ENZYMES** (`ID`, `GENE`, `MW`, `SEQUENCE`,
  `CONC`) and **ENZRXNS** (`ID`, `KCAT`, `SOURCE`, `NOTE`, `EC-NUMBER`, `ENZYMES`), where the `ENZYMES`
  column encodes subunit stoichiometry from `ec.rxn_enz_mat` as `enzyme:count` pairs. Written only when
  the model carries a populated `model.ec`; plain cobra models are unchanged. Export-only — YAML remains
  the round-trippable ecModel format.
* **CI runs on macOS and Windows.** The `test` job gains a lean matrix — ubuntu × Python 3.11/3.12/3.13
  plus macos × 3.12 and windows × 3.12 (5 jobs) — so OS-specific bugs in path handling, `subprocess`,
  and file I/O are caught here rather than in downstream consumers. `lint` / `mypy` / `docs` stay
  ubuntu-only.
* **Docstring rendering fixes.** `set_gam`'s trailing "Returns the (mutated) model for chaining." sat
  inside the `Parameters` block, which griffe read as a phantom parameter named `Returns`; it is now a
  proper `Returns` section. `add_sbo_terms` used 4-space-indented bullets under plain-text headings,
  which CommonMark renders as a code block rather than a list; it now uses bold headers with unindented
  bullets. No API changes.

## 0.2.0 — 2026-06-14

Project rename plus KEGG-reconstruction and CI improvements.

* **Project renamed `raven-python` → `raven-toolbox`.** The import package is now
  `raven_toolbox` (was `raven_python`) and the repository moved to
  `SysBioChalmers/raven-toolbox`. Update imports and any `raven-python` git/URL
  references accordingly. (#34)
* **De-novo KEGG query uses `hmmsearch` instead of `hmmscan`.**
  `get_kegg_model_from_sequences` now runs one `hmmsearch` over the concatenated KO
  library (`-Z` set to the profile count, so E-values and KO assignments are identical
  to the previous `hmmscan` path) — the faster, more parallel search direction.
  `ensure_kegg_hmm_library` no longer runs `hmmpress` (just gunzips); the published
  `.hmm.gz` artefact is unchanged. (#32)
* **Domain-mode `get_kegg_model_for_organism_from_artefacts` auto-resolves the
  taxonomy artefact** from the artefact directory, so `"prokaryotes"` /
  `"eukaryotes"` no longer require an explicit `taxonomy=` path. (#31)
* **Test data no longer ships real KEGG records.** The on-disk
  `tests/data/kegg_dump` is replaced by a session fixture that generates a fully
  fictional KEGG-format dump at runtime, so no KEGG-derived data is redistributed. (#33)
* **Removed the visualization stub and the `[visualization]` extra** — an
  unimplemented placeholder. (#30)
* **CI on Node 24** — `actions/checkout@v5`, `actions/setup-python@v6`. (#35)

## 0.1.0 — 2026-06-10

First release with **published, downloadable KEGG artefacts**, plus a cobra-aligned
hardening pass (no behaviour change on well-formed inputs). Highlights:

* **KEGG artefacts published (`kegg116`):** `ensure_kegg_data` /
  `ensure_kegg_hmm_library` fetch version-pinned, SHA256-verified files from the
  GitHub release. Every artefact is **gzip + version-prefixed**
  (`kegg116_<name>.gz`) so MATLAB and Windows read them with the built-in `gunzip`
  (no external tool) — `organism_gene_ko` moved from xz to gzip for this. The core
  model files (reference model + KO/reaction tables) ship as a single
  `kegg116_core.tar.gz` that `ensure_kegg_data` extracts on first use; the HMM
  libraries and `taxonomy` are separate assets. The **HMM
  libraries ship as one gzip concatenated flatfile per domain**
  (`kegg116_<domain>.hmm.gz`); the client decompresses and `hmmpress`-es once on
  first use, cutting the download ~10× versus the pressed index and letting the
  same artefact serve MATLAB RAVEN.
* **Taxonomy + phylogenetic distance:** publish `kegg116_taxonomy.gz` and add
  `reconstruction.kegg.phyl_dist` (with `PhylDist`), a faithful port of RAVEN's
  `getPhylDist` that regenerates the `keggPhylDist` distance matrix from the
  taxonomy file — so GECKO's organism-distance kcat selection needs no MATLAB
  `.mat`. `ensure_kegg_taxonomy` fetches the artefact.
* **Packaging:** `raven_toolbox.__version__` now derives from the installed package
  metadata (`importlib.metadata`) instead of a hard-coded literal that had drifted
  to `0.0.1`; the docs site reported the wrong version. Pinned `ruff==0.15.15` in
  both the `dev` extra and CI so the lint result is reproducible, and fixed two
  lint errors the unpinned ruff had started flagging.
* **Errors aligned to cobra:** solver/feasibility failures in `run_init`,
  `run_ftinit`, `fill_tasks` and `random_sampling` now raise
  `cobra.exceptions.OptimizationError` (already used elsewhere in the package)
  instead of a bare `RuntimeError`.
* **Consistency:** a single `utils.parse.subsystem_to_str` coerces a reaction
  `subsystem` to cobra's canonical `str` everywhere it is rendered/compared
  (`io.excel`, `comparison.compare`, `curation.batch`, `manipulation.add`) — fixes
  a crash on non-string subsystem items and the silent drop of multi-subsystem
  reactions. GPR score-aggregation (`AGGREGATORS` / `resolve_aggregators`) is now
  shared by `init.score` and `init.genes`. Maintainer-side KEGG-download progress
  uses a module logger instead of `print`.
* **Robustness:** path-traversal guard on bundled-ZIP extraction (`binaries.py`,
  matching the tarfile `filter="data"` precedent); `connect_blocked_reactions`
  rejects a non-positive `penalty`; `random_sampling` refuses a NaN-contaminated
  sample matrix; `ec_data` warns on an all-zero reaction↔enzyme coupling; optional
  `verify=` SHA256 re-check on `ensure_data_file` cache hits; reporter p-value
  guarded against non-finite z-scores. Regression tests added for each.

## 0.1.0a1 — 2026-05-30

First alpha release. Covers the functional scope of RAVEN built on cobrapy:
de-novo reconstruction (KEGG / homology), context-specific modeling (tINIT / ftINIT),
metabolic-task validation, connectivity gap-filling, HPA omics ingestion, sub-cellular
localisation, N-model comparison, reporter metabolites, FSEOF, flux sampling, and the
RAVEN-style I/O formats (YAML / SIF / Excel). Validated against MATLAB RAVEN on Human-GEM
(Jaccard 0.975–0.980).

* **Licensing:** released under the **MIT** license (previously GPL-3.0-or-later).
* **Docs:** Sphinx + MyST documentation site (sources under `docs/`).
* Not yet implemented: visualization (`visualization/`), metabolomics-based (f)tINIT scoring,
  and published binary / KEGG-artefact release bundles. See the README and
  [docs/todo.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/reference/todo.md).

The milestone sections below record the incremental development history leading to this release.

## Infrastructure

* **GitHub Actions CI** ([.github/workflows/ci.yml](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/.github/workflows/ci.yml)) —
  ruff + pytest matrix over Python 3.11/3.12/3.13. Tests that require Gurobi
  auto-skip (no Gurobi on free runners); the known HiGHS upstream blocker
  (`hybrid_interface.Configuration` rejects `lp_method='primal'`) is marked
  `xfail(strict=True)` so CI flips red when optlang fixes it.

## Quality sweep — known-issues section F (design-choice divergences)

Closed the five items in section F (the "design choices that differ from RAVEN"
backlog from the original review). Three docstring/comment fixes; two code
fixes with matching MATLAB back-port proposals in IMPROVEMENTS.md (FS4, B2).

* `run_init` docstring spells out the score-0 semantics divergence between
  classic INIT and ftINIT.
* `get_init_model` inaccurate "same regime" comment replaced with an accurate
  description of the conservative pre-filter.
* `fseof` classifier now uses the slope of `|flux|` (`linregress(enforced, |flux|)`)
  instead of first-vs-last endpoints. A track whose endpoints straddle a
  peak/trough no longer ends up mislabelled.
* `reporter_metabolites` docstring documents the one-sided p-value + z-score
  ordering vs RAVEN's two-tailed sort, and points at the up/down split via
  `gene_fold_changes`.
* `get_elemental_balance` now reports `unknown` for empty-stoichiometry
  reactions (previously vacuously `balanced`). Original review attributed the
  bug to `check_model`; the actual code is in `balance.py`.

Two new regression tests (F3 in `test_analysis_fseof.py`, F5 in
`test_utils_balance.py`). [docs/known_issues.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/reference/known_issues.md) now
fully closed (all sections A–F).

## Quality sweep — known-issues sections C / D / E

Closed all the robustness, efficiency, and dead-code items in one pass.

**Robustness (C):**
* `constrain_reversible_reactions` wraps FVA in try/except + NaN check; both
  backend-raised `OptimizationError` and silent-NaN returns now surface as one
  clear `RuntimeError` (the original `abs(NaN) < eps` silently no-op'd).
* `ensure_binary` downloads through `.part` + `os.replace`, matching `data.py` —
  an interrupted download leaves a `.part`, never a half-complete `.zip`.
* `parse_task_list` (.xlsx) checks `wb.sheetnames` before lookup; missing
  `TASKS` sheet now raises a clear `ValueError` instead of a bare `KeyError`.
* `parse_taxonomy` pads with explicit `""` when a depth level is skipped and
  warns once.

**Efficiency (D):**
* `group_linear_reactions` rewritten with a metabolite worklist (re-enqueue
  the mets touched by each merge); same observable result, O(n+m) work per
  pass instead of restarting the full scan after every merge.
* `parse_kegg_reactions` now caches the parsed stoichiometry on each
  `KeggReaction.stoichiometry`; `build_reference_model` reuses it instead of
  re-parsing.

**Dead code (E):**
* Dropped `KeggReaction.modules` and `.rhea` (parsed but never consumed).
* Dropped the vestigial `only_genes_in_models` parameter from `_ortholog_map`.

Six new regression tests; the only one without a test is the `.part` atomic
download (defensive, needs urlopen mocking).

## Quality sweep — known-issues section B

Closed all four "silent misbehaviour" items from [docs/known_issues.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/reference/known_issues.md):
* `merge_models` warns on `formula` / `charge` conflicts when two source models
  share a name[comp] but disagree (used to silently keep the first-seen).
* `add_reactions_from_equations` warns when creating a metabolite in an
  unregistered compartment — both the `mets_by="id"` and `mets_by="name"` paths
  (id-mode used to skip the check entirely, an asymmetry).
* `parse_task_list` warns when continuation data appears before any task ID
  has been seen (used to silently drop the orphan row).
* `export_model_to_sif` warns up front when a custom label map sends two
  distinct ids to the same label (used to silently collapse nodes).
Four new regression tests cover them.

## Quality sweep — known-issues section A

Closed all six "latent edge-case bug" items from [docs/known_issues.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/reference/known_issues.md):
* `add_reactions_from_equations` no longer misparses `"2 oxoglutarate"` (or any
  leading-number metabolite name) — the resolver tries the full token before
  splitting off a coefficient.
* `add_reactions_from_equations` warns when an equation's terms cancel to a
  zero-metabolite reaction.
* `add_reactions_from_model` tracks ids minted within the batch so two source
  metabolites whose ids both collide with the draft don't collapse onto the
  same generated id.
* `add_transport_reactions` warns on duplicate metabolite names in the source
  or target compartment instead of silently dropping all but one.
* `connect_blocked_reactions` membership-guards the FVA result before
  `.at[]` lookup.
* `assign_kos` rejects `cutoff >= 1` up front — would have crashed inside the
  ratio filter at `log(best_evalue) == 0`.
Six new regression tests cover the user-reachable cases.

## Phase 7 — Localization

* **Sub-cellular localisation by MILP.** [`localization.predict_localization`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/localization/predict.py)
  + [`apply_localization`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/localization/predict.py). Deterministic (not simulated
  annealing); caller-passed `reactions_to_relocate` set with everything else pinned;
  incomplete-model tolerant (no silent reaction removal); `apply=False` returns a diff
  preview; multi-compartment by default with primary-free, extras-penalised scoring.
* **Predictor loaders.** [`load_wolfpsort`, `load_deeploc`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/localization/scores.py),
  with the `gene × compartment` DataFrame contract open for any predictor.
* **Compartment helpers** ([`manipulation/compartments.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/compartments.py)):
  `merge_compartments`, `copy_to_compartment` — useful standalone for model curation.
* **Real-data validation on yeast-GEM** ([docs/yeast_localization_benchmark.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/studies/yeast_localization_benchmark.md))
  — accuracy 0.72 → 0.39 on 298 GPR'd reactions as confident predictor mis-scoring rises
  from 0 % to 50 %; perfect on compartments with disjoint gene sets (c/g/lp/p/v/vm), and
  surfaces a `transport_cost` calibration insight for soft-probability score tables.

## Phase 5 — Data integration & analysis

* **Reporter metabolites, FSEOF, random sampling** ([`analysis/`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/analysis/)).
* **HPA omics ingestion** ([`omics.parse_hpa`, `parse_hpa_rna`, `hpa_gene_scores`, `rna_gene_scores`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/omics/hpa.py))
  — pandas-tidy DataFrames replace RAVEN's sparse-matrix layout; scoring adapters reuse the
  existing GPR walk.
* **N-model comparison** ([`comparison.compare_models`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/comparison/compare.py)).
* **Dynamic FBA** is **not ported** — established Python packages cover it (`dfba`,
  `reframed`, `mewpy`).

## Phase 4d — ftINIT

* **ftINIT pipeline** ([`init.ftinit`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/ftinit.py)) — staged MILP, linear merge,
  task-aware gap-filling, gene pruning.
* **Validated against MATLAB RAVEN on Human-GEM.** 5 Hart2015 cell-line models;
  Jaccard 0.973–0.977 (no-task) and 0.978–0.980 (task-constrained). See
  [docs/humangem_validation.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/studies/humangem_validation.md).
* **Parameter calibration & input-robustness study** ([docs/init_param_calibration.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/studies/init_param_calibration.md))
  — `mip_gap=0.01` is the genome-scale full-pipeline sweet spot (~37% faster than 0.001 at
  Jaccard 0.995); pipeline is robust to expression noise (Jaccard 0.92–0.95) but sensitive
  to sparsity (50–70% dropout → Jaccard 0.59–0.71); the task + gap-fill layer keeps the
  essential-task pass-rate at 67–69/69 across the gradient, whereas tINIT-without-it passes
  only 35/69 even on clean data.
* **Cross-solver portability** ([docs/init_solver_benchmark.md](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/docs/studies/init_solver_benchmark.md))
  + [`tests/test_init_solvers.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/tests/test_init_solvers.py): Gurobi and GLPK pass at toy
  scale; only Gurobi is viable at genome scale today (HiGHS hits an upstream optlang
  `clone()` bug; GLPK ignores `configuration.timeout` on MIP).
* **Engineering wins surfaced by the genome-scale work:** `check_tasks` and
  `fill_tasks._feasible` rewritten in-place (~12× each); `optlang.symbolics.add` builds
  in the MILP construction (the O(n²) sympy `sum()` blow-up was the original genome-scale
  blocker); bounded gap-fill MILP; `rescaleModelForINIT` ported.

## Phase 4c — tINIT

* **INIT MILP and the tINIT pipeline** ([`init.run_init`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/init.py),
  [`init.get_init_model`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/init/build.py)). Clean optlang reformulation;
  RNA-seq scoring via `5·ln(level/ref)`-clamped.

## Phase 4b — Gap-filling

* **Connectivity gap-filling** ([`gapfilling.connect_blocked_reactions`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/gapfilling/fill.py))
  — MILP. Targeted (toward objective) mode delegates to `cobra.gapfill`.

## Phase 4a — Metabolic tasks

* **Task list parsing + `check_tasks`** ([`tasks/`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/tasks/)).

## Phase 3 — Reconstruction

* **Homology-based draft** from a template GEM + BLAST/DIAMOND wrappers
  ([`reconstruction/homology/`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/reconstruction/homology/)) — with structured
  improvements over RAVEN's `getModelFromHomology` (see IMPROVEMENTS H1–H6).
* **KEGG five-step pipeline** ([`reconstruction/kegg/`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/reconstruction/kegg/)):
  dump → parser → HMM library builder → species model → HMM-query draft.
* **MetaCyc reconstruction** **not ported** (and flagged for removal from MATLAB RAVEN —
  see IMPROVEMENTS R-MetaCyc).

## Phase 2 — I/O

* **YAML** aligned to cobra's `!!omap` writer + RAVEN-only fields preserved into `.notes`,
  plus geckopy `ec-*` for enzyme-constrained models
  ([`io/yaml.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/io/yaml.py)).
* **SIF**, **Excel export**, and **Standard-GEM `model/<fmt>/…` git layout**
  ([`io/`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/io/)). Excel import intentionally excluded.

## Phase 1 — Foundation

* **GPR / balance / validation / parsing helpers** ([`utils/`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/utils/)) —
  cobra-absent bits only; the rest are cheatsheeted.
* **Manipulation ergonomic layer** ([`manipulation/`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/)) —
  add/change/remove/transport/transfer/merge/simplify/variance + adopted transforms.
* **External-binary resolver** ([`binaries.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/binaries.py)) — version-pinned
  release-ZIP registry, SHA256-verified cache.

## Phase 0 — Scaffold

* Project structure, packaging, pytest skeleton, license alignment with MATLAB RAVEN
  (GPL-3.0-or-later).
