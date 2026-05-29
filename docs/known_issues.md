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

- **`manipulation/merge.py` — `merge_models` (`match_by="name"`):** metabolites with
  the same `name[comp]` but differing `formula`/`charge` across models silently unify
  to the first-seen (no warning). *Fix:* warn on a formula/charge conflict.
- **`manipulation/add.py`:** with `allow_new_mets=True`, an unknown/typo'd explicit
  compartment silently creates a metabolite in a brand-new compartment (the
  `mets_by="id"` path never validates compartment at all — asymmetric).
- **`tasks/tasklist.py` — `parse_task_list`:** continuation rows appearing *before* the
  first ID-bearing row are silently dropped (`current is None`). Undocumented.
- **`io/sif.py` — `export_model_to_sif`:** node dedup is label-based, so two distinct
  metabolites/reactions mapped to the same custom label collapse into one edge target.
  Only with custom label maps.

## C. Robustness gaps

- **`manipulation/simplify.py` — `constrain_reversible_reactions`:** uses
  `fraction_of_optimum=0` FVA; on an infeasible model the NaN ranges make the
  `abs(lo) < eps` comparisons silently no-op rather than erroring.
- **`binaries.py` — `ensure_binary`:** an interrupted download leaves a stale
  `_download.zip` in the cache (self-healing on retry, since the next call overwrites,
  but it lingers on error). `data.py` does this better with a `.part` temp + atomic
  `replace`. Also, a cached binary's integrity isn't verified before reuse.
- **`tasks/tasklist.py`:** the xlsx reader assumes a sheet literally named `TASKS`;
  any other name raises a bare `KeyError`. (The whole xlsx path is also untested.)
- **`reconstruction/kegg/taxonomy.py`:** depth handling assumes no skipped levels
  (a `####` directly under a `##`); real KEGG taxonomy is well-formed, and the domain
  classification consumed downstream stays correct, so robustness-only.

## D. Efficiency (correct but slow at scale)

- **`manipulation/simplify.py` — `group_linear_reactions`:** restarts the full
  scan after *every* merge (O(n²·m) on large models). (Note: `init/merge.py`'s
  `merge_linear` is the ftINIT-grade implementation; this `simplifyModel` variant is
  the lossy gene-dropping one.) *Fix:* don't `break` on each merge, or maintain live
  incidence as `merge_linear` discusses.
- **`reconstruction/kegg/parse.py`:** `_parse_equation` runs once per reaction in
  `parse_kegg_reactions` and again in `build_reference_model` — a full redundant parse.
  Maintainer-only, run-once, so low impact.

## E. Dead / vestigial code

- **`reconstruction/kegg/parse.py`:** `KeggReaction.modules` (and `rhea`) are parsed
  and stored but never written to any artefact — dead data collection.
- **`reconstruction/homology/homology.py`:** `only_genes_in_models` is threaded into
  `_ortholog_map`'s signature but never referenced there (the filtering happens
  earlier); vestigial.

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
