# Known issues — deferred findings from the full-codebase review

Low-priority issues found during the 2026-05-26 full critical review (see PROGRESS
commits `db5a1fa`, `2b85830` for the bugs that *were* fixed). These are deliberately
not fixed yet: each is an edge case, robustness gap, efficiency concern, dead code, or
a documented design choice — none affects correctness on normal, well-formed inputs.
Line numbers are indicative; refer to the named function.

## A. Latent edge-case bugs (rare crashes / silent misbehaviour)

- **`manipulation/add.py` — `add_reactions_from_equations` (name mode):** a metabolite
  whose *name* starts with a number (e.g. `"2 oxoglutarate"`) is misparsed — the
  leading number is taken as a stoichiometric coefficient, silently corrupting the
  equation. RAVEN's parser has the same fragility. *Fix:* in `mets_by="name"` mode,
  only treat a leading token as a coefficient if the remainder still resolves to a
  known metabolite, or require explicit numeric+space only when the rest is a met.
- **`manipulation/add.py` — `add_reactions_from_equations`:** an equation whose terms
  net to zero produces a reaction with no metabolites and no warning. *Fix:* warn or
  skip empty reactions.
- **`manipulation/transfer.py` — `add_reactions_from_model` (`_new_met_id`):** two
  source metabolites sharing an `id` but different `name[comp]`, neither already in the
  draft, both pass the "not in model" check and get assigned the same new id →
  `add_metabolites` collision/crash. *Fix:* track already-assigned new ids within the
  batch and route intra-batch collisions through `_new_met_id` too.
- **`manipulation/transport.py` — `add_transport_reactions`:** the source-metabolite
  lookup is keyed by name (`{m.name: m}`), so two metabolites sharing a name in the
  source compartment silently collapse (one is dropped from transport). *Fix:* group
  by name → list, or key by id.
- **`gapfilling/fill.py` — `connect_blocked_reactions`:** `fva.at[r, "maximum"]`
  assumes every candidate appears in the FVA index; a `KeyError` results if FVA drops
  one. *Fix:* `.get`/membership guard.
- **`reconstruction/kegg/query.py` — `assign_kos`:** divides by `log(best_evalue)`,
  which is `0` when the best E-value in a group is exactly `1.0` → `ZeroDivisionError`.
  Guarded in practice (default `cutoff=1e-30` excludes `evalue==1`), but reachable when
  a caller passes `cutoff >= 1`. *Fix:* clamp the cutoff `< 1` or special-case
  `log_best == 0`.

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
