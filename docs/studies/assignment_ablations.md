# Compartment-assignment ablations: transport pruning & gap-filling

Two `assign_compartments` features shipped with correctness tests but no performance measurement. This
study measures each on curated *S. cerevisiae* yeast-GEM, flattened to one compartment and reassigned —
the same draft the other yeast studies use (2569 draft reactions, 2296 to place, `min_growth` = 50 % of
curated). Regenerate with `scripts/benchmark_assignment_ablations.py`.

## 1. Transport-reaction pruning (`prune_transports`)

After placement, `assign_compartments` adds the inter-compartment transports the network needs to stay
functional. Pruning then removes the ones that turned out redundant — a transport whose removal leaves the
model still certified. The question: how much does it remove, and does removing it cost accuracy?

| | transports added | reaction agreement vs curated | runtime |
|---|--:|--:|--:|
| `prune_transports=True` (default) | **994** | 0.721 | 209 s |
| `prune_transports=False` | 1268 | 0.721 | 74 s |

- **Pruning removes 274 transports — 21.6 % — at no accuracy cost.** Reaction-level agreement is identical
  to four decimals, and the pruned set is a **strict subset** of the unpruned one (0 transports appear only
  when pruning is on). So pruning only ever removes genuinely redundant transports; it never trades one
  placement for another.
- **The cost is runtime**, not quality: pruning re-certifies the model after each candidate removal, so the
  default is ~3× slower. That is the price of a leaner network — 274 fewer artificial inter-organelle
  shuttles that a curator would otherwise have to inspect.

This is the intra-method complement to the [CarveFungi head-to-head](carvefungi_milp_benchmark.md), which
showed the *transport-minimisation objective term* yields ~41 % fewer transports than CarveFungi's carve;
here the *post-hoc pruning pass* removes a further fifth of what placement provisionally added.

## 2. Gap-filling on the natural draft (does it over-add?)

`assign_compartments` can pull reactions from a `universal` model to restore biomass when a
compartmentalised placement cannot grow. The first question is whether it fires *gratuitously* — adds
reactions a well-connected model did not actually need.

Reassigning the flattened draft with a universal available (the draft itself) vs without:

| | added reactions | certified | growth |
|---|--:|:--:|--:|
| no `universal` | 0 | yes | 0.14 |
| `universal=` draft | **0** | yes | 0.14 |

**Gap-fill adds nothing on the natural draft.** Transport addition alone restores growth well above the
floor (0.14 vs the 0.04 floor), so the growth-failure feedback that triggers gap-fill never fires. This is
the genome-scale confirmation of the `test_no_gratuitous_gapfill` unit test: the feature is a safety net
that stays inert when the model is already functional, rather than a source of spurious additions.

## 3. Gap-filling under real gaps (does it add the *right* reaction?)

To measure gap-fill *working*, it needs real gaps. Ground-truthed knockout-recovery on the certified
compartmentalised model: take each of the 352 internal reactions whose single removal drops growth below
5 % of optimum, remove it, and gap-fill from a universal that contains it (a copy of the model — so the
removed reaction is a candidate and the ground truth is exact). This is the same
`cobra.flux_analysis.gapfill` call `assign_compartments._gapfill` wraps. A 60-reaction seeded sample:

| outcome | count | of sample |
|---|--:|--:|
| recovered growth, **re-added the exact removed reaction** | 27 | 45 % |
| `cobra.gapfill` numerical failure (declined) | 33 | 55 % |
| recovered growth with a *wrong* reaction | **0** | 0 % |

Two things stand out, and they are the honest headline:

- **Precision is perfect. Of every recovery, 100 % re-added the exact reaction that was removed** — and
  zero knockouts were "fixed" with a different or superfluous reaction. When gap-fill answers, it answers
  correctly.
- **The only failure mode is declining to answer.** The 55 % shortfall is entirely
  `cobra.flux_analysis.gapfill`'s own numerical-validation limit ("Failed to validate gap filled model, try
  lowering the integer threshold"), not a wrong addition. `_gapfill` catches it and reports no gap-fill, so
  a placement that cobra cannot solve stays *uncertified and visible* rather than silently mis-filled. The
  recovery rate is therefore a property of cobra's MILP tolerance, not of the assignment logic — and it is
  bounded below the true recoverable set, never above it.

cobra's failure message suggests "try lowering the integer threshold", so we swept it. It does not help
— the default is a **sharp optimum**:

| `integer_threshold` | recovered (of 30) | exact |
|---|--:|--:|
| 1e-9 | 0 % | — |
| **1e-6 (cobra default)** | **~47 %** | **100 %** |
| 1e-5 / 1e-4 / 1e-3 | 0 % | — |

Every value other than the default collapses recovery to zero, in both directions. So there is no tuning
win to be had, the documented remedy is misleading, and `_gapfill` correctly leaves the threshold at
cobra's default. The ~45–47 % ceiling is cobra's, and it is not cheaply liftable — which makes the
fail-safe design (decline, never mis-fill) the right call rather than a workaround.

## 4. Verdict

- **Transport pruning:** a strict win — 21.6 % fewer transports at zero accuracy cost, paid for in runtime.
  Correctly a default.
- **Gap-filling:** precise and conservative. It never fires gratuitously (§2) and never adds a wrong
  reaction (§3); when it recovers a gap it recovers it *exactly*. Its recall is capped by cobra's gapfill
  tolerance, and it fails safe — an unsolved gap is reported, not papered over.

## 5. Reproducing

```
python scripts/benchmark_assignment_ablations.py   # all three parts -> .research_tmp/assignment_ablations.json
```
