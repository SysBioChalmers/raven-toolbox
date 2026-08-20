# ftINIT extraction determinism — `strict_gap` and `canonical`

What the two opt-in determinism flags on `ftinit()` actually buy on a genome-scale model,
and — the reason this page exists — what they *do not* buy. Companion to the flags'
API documentation in `raven_toolbox.init.ftinit`.

## The problem

The ftINIT extraction MILP is highly degenerate: many different reaction subsets reach the
same score optimum, and the solver returns whichever one its branch-and-cut happens to
land on. raven-toolbox already matches RAVEN's reproducibility levers — `Threads=1` (the
dominant one; multi-threaded Gurobi picks among equal optima non-deterministically) and a
fixed `Seed` — but those make the choice **reproducible**, not **unique**. The chosen
optimum is a property of the solver build, so a different Gurobi version, platform, or
seed lands on a different one. Downstream that surfaces as gene-essentiality predictions
that drift for reasons unrelated to the model or the data.

The two flags attack that directly:

* **`strict_gap`** — one solve per step proven to a fixed *absolute* gap (0.05, below the
  0.1 reaction-score granularity), replacing RAVEN's loose relative-gap escalation, whose
  final near-zero-objective run otherwise accepts an arbitrary within-gap incumbent.
* **`canonical`** — a lexicographic phase 2 over the removable (binary) reactions: hold the
  score objective at its optimum, minimise the count of kept reactions, then their summed
  id rank. Applied to both each extraction step and the task gap-fill.

Both default to off, so the shipped behaviour is unchanged.

## Method

* **Model / data.** Human-GEM v2.0.0, Hart2015 DLD1 cell line, the same prep and scoring
  as the [Human-GEM validation](humangem_validation.md) run (series `1+1`, tasks on).
* **Solver.** Gurobi 13, `Threads=1`, single machine.
* **Determinism probe.** The extraction seed (`ftinit._EXTRACT_SEED`, RAVEN's 1234) is
  changed to 7 and the model rebuilt. Everything else — inputs, code, solver build,
  machine — is held fixed, so any difference between the two models is *purely* the
  solver's tie-break among equal optima. This is a lower bound on cross-version drift: a
  Gurobi upgrade perturbs the search far more than a seed change does.
* **Metrics.** (1) *reaction seed-swing* — the size of the symmetric difference between the
  two kept-reaction sets; (2) *essential-gene seed-swing* — the symmetric difference
  between the two single-gene-deletion essential sets. The second is the metric that
  matters for users; the first is the one the flags act on.

Each cell below is one build per seed (two builds per row). This is an n=1 probe on one
cell line with one seed pair, not a replicated benchmark — see *Limitations*.

## Results

| Configuration | Reaction seed-swing | Essential-gene seed-swing | Build time |
|---|---:|---:|---:|
| Baseline (no flags) | 125 | **5** | ~12 min |
| `canonical` only | 64 | not measured | ~25 min |
| `canonical` + `strict_gap` | **34** | 19 | ~34–80 min |

Out of ~7766 kept reactions and ~281 essential genes. The five baseline essentiality flips
are PC, ACADVL, DUT, ALDH3A2 and RDH5. Essential-gene counts: 281 / 278 for the two
baseline seeds, 283 / 296 for the two `canonical`+`strict_gap` seeds.

## Findings

**The flags do what they claim at the reaction level.** The seed-swing drops 125 → 64 with
`canonical` alone and → 34 with both, a 3.7× reduction. The remaining 34 reactions are the
degeneracy the lexicographic phase 2 cannot see: it canonicalises the *removable-reaction
choice*, while a residual flux-distribution degeneracy (which feeds the next step's small
`ess_force` clamp) is left to the solver.

**But the reaction metric is a poor proxy for what users care about.** All 34 residual
reactions are GPR-less — mostly transport. They move without touching any gene, so driving
that number down does not, by itself, stabilise anything gene-level.

**And gene-essentiality determinism got *worse*, not better: 5 → 19 flips.** Two effects
compound:

1. `canonical` minimises the kept set, so it produces **sparser** models (283 / 296
   essential genes vs the baseline's 281 / 278). A sparser network has fewer alternative
   routes, so more genes sit on a unique path and become essential — and more of them are
   sensitive to whichever residual alternative pathway survives. The transport-level
   degeneracy that the reaction metric dismisses as "GPR-less" propagates into enzymatic
   genes' essentiality through connectivity.
2. `strict_gap` **cannot prove the genome-scale optimum in reasonable time.** One build in
   this probe accumulated ~67 min of per-step timeouts and fell back to arbitrary
   incumbents — exactly the behaviour the flag exists to eliminate, and it adds noise of
   its own. The 34–80 min build-time range in the table is that fallback showing up as
   variance.

The runtime cost is real regardless: 2× for `canonical`, 3–7× for both.

## Conclusions

* These are a correct, tested, opt-in tool for a **more parsimonious and more reproducible
  extraction**. They pin *which* of many equally data-consistent optima is returned; they
  do not make the model biologically more accurate, and they reduce the MILP's fragility
  rather than its error.
* They are **not** a fix for gene-essentiality determinism, and on this benchmark they made
  it worse. Do not reach for them to stabilise a downstream essentiality analysis.
* **For run-to-run essentiality reproducibility the reliable lever is pinning the solver
  stack** — raven-toolbox commit plus `gurobipy` version — which gives bit-identical
  results without any of the runtime cost.
* Reasonable use: producing one or a few stable, parsimonious model artifacts where the
  extra build time is affordable and the parsimony is wanted for its own sake. For most
  work the baseline is fine.
* Worth testing if the flags are pursued further: `canonical` **without** `strict_gap`
  (untested at gene level here) isolates parsimony from the timeout noise, and is the
  cheaper of the two. Fixing the essentiality drift likely needs the residual transport
  degeneracy canonicalised too, not a tighter primary gap.

## Limitations

One cell line, one seed pair, one build per configuration, one machine — enough to
establish the direction and rough magnitude of each effect, not to put error bars on any
number. The `canonical`-only gene-essentiality figure was never measured, so the
attribution of the 5 → 19 regression between the two flags is inferred from the mechanism,
not measured. Cross-Gurobi-version drift — the failure mode the flags actually target — was
not measured at all; the seed probe stands in for it.

## Reproducing

Requires Gurobi and the Human-GEM work directory from the
[Human-GEM validation](humangem_validation.md) run.

```python
import cobra
from raven_toolbox.init import ftinit
from raven_toolbox.init import ftinit as ftinit_mod

def build(seed, **flags):
    ftinit_mod._EXTRACT_SEED = seed          # determinism probe; not a public parameter
    model = ftinit(prep, rxn_scores=scores, time_limit=900, **flags)
    return {r.id for r in model.reactions}, model

for label, flags in [("baseline", {}),
                     ("canonical", {"canonical": True}),
                     ("canonical+strict", {"canonical": True, "strict_gap": True})]:
    a, ma = build(1234, **flags)
    b, mb = build(7, **flags)
    print(label, "reaction swing:", len(a ^ b))
    ea = {g.id for g in cobra.flux_analysis.find_essential_genes(ma)}
    eb = {g.id for g in cobra.flux_analysis.find_essential_genes(mb)}
    print(label, "essential-gene swing:", len(ea ^ eb), f"({len(ea)} vs {len(eb)})")
```

Budget several hours: the `canonical`+`strict_gap` pair alone is 1–3 h of solve time.
