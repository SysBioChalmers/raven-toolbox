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

- **`init/init.py` — `run_init`:** a reaction with score *exactly* 0 is removable
  (gets a binary with 0 reward), unlike ftINIT's "score 0 ⇒ left in the model". This
  matches classic INIT semantics but is undocumented in `run_init`; worth a one-line
  note to avoid cross-variant confusion.
- **`init/build.py` — `get_init_model`:** the dead-end pre-filter uses
  `open_exchanges=True` even when `allow_excretion=False`; harmless (INIT drops the
  reaction anyway) but the "same regime run_init will use" comment is inaccurate for
  the strict path.
- **`analysis/fseof.py` — `fseof`:** the amplify/knockdown/knockout label is decided
  from the first vs last enforced-flux endpoints, not the (already-computed) regression
  slope; a correlated-but-non-monotone track can be mislabelled. The `|r|` gate limits
  this.
- **`analysis/reporter.py` — `reporter_metabolites`:** reports a one-sided ("up")
  enrichment p-value and sorts by `z_score`, which differs from RAVEN's p-value
  ordering. Internally consistent; flagged for parity awareness.
- **`utils/validate.py` — `check_model`:** a reaction with no metabolites is reported
  `balanced` (vacuously); arguably should be `unknown`.
