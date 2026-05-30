# Known issues — deferred findings from the full-codebase review

Low-priority issues found during the 2026-05-26 full critical review (see PROGRESS
commits `db5a1fa`, `2b85830` for the bugs that *were* fixed). These are deliberately
not fixed yet: each is an edge case, robustness gap, efficiency concern, dead code, or
a documented design choice — none affects correctness on normal, well-formed inputs.
Line numbers are indicative; refer to the named function.

## A. Latent edge-case bugs

All six items in this section were closed in a quality-sweep pass (see CHANGELOG
"Quality sweep" entry); regression tests live alongside each fixed function. Kept
here for traceability of the original review.

- ✅ **`manipulation/add.py` — `add_reactions_from_equations` (name mode):** a
  metabolite whose *name* starts with a number (e.g. `"2 oxoglutarate"`) was
  misparsed — the leading number was taken as a coefficient. Fixed by trying the
  full token as a name first and only splitting off a coefficient when the
  remainder names something resolvable. Test:
  `tests/test_manipulation_add.py::test_name_mode_preserves_leading_number_name`.
- ✅ **`manipulation/add.py` — `add_reactions_from_equations`:** an equation whose
  terms net to zero produced a reaction with no metabolites and no warning. Now
  warns. Test: `test_empty_stoichiometry_warns`.
- ✅ **`manipulation/transfer.py` — `add_reactions_from_model` (`_new_met_id`):**
  two source metabolites sharing an `id` but different `name[comp]`, both
  needing minted ids, used to collide. Now tracks ids minted in the batch.
  Test: `tests/test_manipulation_transfer.py::test_intra_batch_id_minting_unique`.
- ✅ **`manipulation/transport.py` — `add_transport_reactions`:** the
  source/target metabolite lookup was keyed by name, silently dropping
  same-name duplicates. Now warns on collision. Test:
  `tests/test_manipulation_transport.py::test_duplicate_name_in_source_compartment_warns`.
- ✅ **`gapfilling/fill.py` — `connect_blocked_reactions`:** the
  `fva.at[r, "maximum"]` access used to crash with `KeyError` if FVA dropped a
  candidate. Now membership-guarded (defensive — the original is unreachable
  with cobra's default FVA, no regression test).
- ✅ **`reconstruction/kegg/query.py` — `assign_kos`:** `cutoff >= 1` would let
  `log(best_evalue) == 0` through and crash inside the ratio filter. Now
  rejected up front with a clear error. Test:
  `tests/test_reconstruction_kegg_query.py::test_cutoff_ge_one_rejected`.

## B. Silent misbehaviour on unusual inputs

All four items in this section were closed in a quality-sweep pass (see CHANGELOG
"Quality sweep — known-issues section B" entry); regression tests live alongside
each fixed function. Kept here for traceability of the original review.

- ✅ **`manipulation/merge.py` — `merge_models`:** name[comp]-matched metabolites
  with differing `formula` or `charge` across models now warn instead of silently
  unifying to the first-seen. Tests
  `tests/test_manipulation_merge.py::test_formula_conflict_warns`,
  `test_charge_conflict_warns`.
- ✅ **`manipulation/add.py`:** both the `mets_by="id"` and `mets_by="name"` paths
  now warn when a new metabolite is created in a compartment that hasn't been
  registered yet (the id-mode path used to skip the check entirely). Tests
  `tests/test_manipulation_add.py::test_id_mode_unknown_compartment_warns`,
  `test_name_comp_unknown_compartment_warns`.
- ✅ **`tasks/tasklist.py` — `parse_task_list`:** continuation rows appearing
  *before* the first ID-bearing row used to be silently dropped. Now warns
  with the file:row pointer. Test
  `tests/test_tasks.py::test_parse_warns_on_data_row_before_first_id`.
- ✅ **`io/sif.py` — `export_model_to_sif`:** when the caller's custom label
  map sends two distinct ids to the same label, the target-side dedup used
  to silently merge the nodes. Now warns up front. Test
  `tests/test_io_sif.py::test_collapsing_label_map_warns`.

## C. Robustness gaps

All four items closed in the same quality sweep (see CHANGELOG); regression
tests live alongside each fixed function.

- ✅ **`manipulation/simplify.py` — `constrain_reversible_reactions`:** the
  FVA call is now wrapped in a try/except + NaN check; both backend-raised
  `OptimizationError` and silent-NaN returns surface as a single clear
  `RuntimeError`. Test
  `tests/test_manipulation_simplify.py::test_constrain_reversible_raises_on_infeasible`.
- ✅ **`binaries.py` — `ensure_binary`:** downloads through a `.part` sibling
  and `os.replace`s into the final name on success, mirroring `data.py`.
  An interrupted download leaves a `.part` (never a half-written `.zip`).
  Defensive — no regression test (needs urlopen mocking).
- ✅ **`tasks/tasklist.py`:** the xlsx reader checks `wb.sheetnames` before
  the `wb["TASKS"]` lookup; a missing sheet now raises a clear `ValueError`
  listing the actual sheets. Test
  `tests/test_tasks.py::test_parse_task_list_xlsx_missing_tasks_sheet`.
- ✅ **`reconstruction/kegg/taxonomy.py`:** depth handling pads with explicit
  `""` placeholders and warns once when a level is skipped (e.g. `####`
  directly under `##`). Test
  `tests/test_reconstruction_kegg_hmm.py::test_parse_taxonomy_handles_skipped_depth`.

## D. Efficiency (correct but slow at scale)

- ✅ **`manipulation/simplify.py` — `group_linear_reactions`:** rewritten with
  a metabolite worklist (re-enqueue the mets touched by each merge) instead of
  the restart-after-every-merge loop. Same observable behaviour, O(n+m) work
  per pass instead of O(n²·m). Test
  `tests/test_manipulation_simplify.py::test_group_linear_merges_long_chain_in_one_pass`.
- ✅ **`reconstruction/kegg/parse.py`:** `parse_kegg_reactions` now caches the
  parsed stoichiometry on each `KeggReaction.stoichiometry`; `build_reference_model`
  reuses it instead of re-parsing. Test
  `tests/test_reconstruction_kegg_parse.py::test_stoichiometry_cached`.

## E. Dead / vestigial code

- ✅ **`reconstruction/kegg/parse.py`:** removed `KeggReaction.modules` and
  `.rhea` (parsed but never consumed by the artefact builders).
- ✅ **`reconstruction/homology/homology.py`:** removed the vestigial
  `only_genes_in_models` parameter from `_ortholog_map` (the actual filtering
  happens earlier in `get_model_from_homology`).

## F. Documented design choices that differ from RAVEN (not bugs)

All five items addressed in the quality-sweep pass (see CHANGELOG). The three
docstring/comment items got documentation fixes; the two correctness items
(F3, F5) got code fixes plus matching RAVEN-back-port proposals in
[IMPROVEMENTS.md](improvements.md) (FS4, B2).

- ✅ **`init/init.py` — `run_init`:** docstring now spells out the score-0
  semantics divergence between classic INIT (score-0 = removable) and ftINIT
  (score-0 = kept by default).
- ✅ **`init/build.py` — `get_init_model`:** the inaccurate "same regime
  run_init will use" comment is replaced — the pre-filter is now correctly
  described as conservative (only reactions blocked under the *most* permissive
  regime are removed, so the strict run_init path never loses a candidate).
- ✅ **`analysis/fseof.py` — `fseof`:** classifier now uses the slope of
  `|flux|` (a real `linregress(enforced, |flux|)` fit) instead of the first vs
  last endpoints. A track whose endpoints straddle a peak/trough no longer
  ends up with the wrong amplify/knockdown label. Matching MATLAB back-port
  proposed as IMPROVEMENTS FS4. Test
  `tests/test_analysis_fseof.py::test_amplify_label_uses_abs_slope_not_endpoint_difference`.
- ✅ **`analysis/reporter.py` — `reporter_metabolites`:** docstring now
  documents the one-sided ``1 - Φ(z)`` `p_value` and the `z_score`-descending
  sort, explicitly contrasted with RAVEN's two-tailed `allPValues` ordering.
  Internally consistent; same information available via the up/down `gene_fold_changes`
  partition.
- ✅ **`utils/balance.py` — `get_elemental_balance`:** empty-stoichiometry
  reactions now report `unknown` instead of `balanced`. (Original review
  attributed the bug to `utils/validate.py::check_model`; the actual code
  lives in `balance.py`.) Matching MATLAB back-port proposed as
  IMPROVEMENTS B2. Test
  `tests/test_utils_balance.py::test_empty_reaction_is_unknown`.
