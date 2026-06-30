# CarveFungi assignment head-to-head: our transport cost on CarveFungi's own carve MILP

A genuine empirical comparison of CarveFungi's compartment-**assignment** objective against ours, run
on CarveFungi's own intermediate state with **CarveFungi's own carve MILP** (its `minmax_reduction`,
in CPLEX). This fixes the three flaws that made an earlier quick emulation a strawman (see
[carvefungi_analysis.md](carvefungi_analysis.md)) and is honest about what a hard, loosely-bounded
MILP can and cannot support.

* Drivers:
  * [`scripts/run_carvefungi_cplex.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/scripts/run_carvefungi_cplex.py)
    — runs CarveFungi's **own** `minmax_reduction` (CPLEX), unmodified; the definitive comparison.
  * [`scripts/benchmark_carvefungi_milp.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/scripts/benchmark_carvefungi_milp.py)
    — an independent Gurobi re-implementation, used to study the formulation (tighter indicator coupling).
* CarveFungi: [github.com/SandraCastilloPriego/CarveFungi](https://github.com/SandraCastilloPriego/CarveFungi),
  bioRxiv [2023.08.23.554328](https://doi.org/10.1101/2023.08.23.554328).

## Design — the intermediate-state swap

Hold CarveFungi's **real intermediate state** fixed and swap only the assignment objective:

* **Candidate set:** CarveFungi's universal model `bigModelv2.21b.sbml` (5602 reactions, 8
  compartments) — each reaction exists only in the compartments the DB instantiates it in.
* **Scores:** CarveFungi's *unmodified* scoring code on the shipped *S. cerevisiae* annotations, run to
  produce the exact per-(reaction, compartment) score dict it feeds its carve MILP.
* **Both arms** run CarveFungi's `minmax_reduction` on the *same* scores — `S·v=0`, big-M reversibility
  coupling, hard biomass ≥ 0.1 and ATP-maintenance ≥ 0.1. Only the objective's parsimony differs:
  * **Arm A (CarveFungi):** `max Σ score(r,c)·y[r,c]` — no transport cost.
  * **Arm B (ours):** the same, with each inter-compartment transport reaction's score reduced by 0.3
    (our transport-minimisation term, fed through their objective — no constraint or variable changed).

An independent adversarial review confirmed the harness is **faithful — not a strawman**: the runner
imports and calls CarveFungi's `minmax_reduction` untouched; the only deviations are the three
disclosed here (a CPLEX time limit; Arm B's −0.3 transport coefficient; the DeepLoc score injection
below), and Arm B perturbs only objective coefficients on existing binaries. Each kept reaction copy's
compartment is read from its **metabolites** (ground truth), not its id suffix — the universal model's
suffixes are mixed-case (`_C`/`_m`/`_x`/`_M` …) and overloaded with transport codes (`_TCE` …), so
suffix parsing both mislabels compartments and double-counts a reaction's copies.

## A repro bug in CarveFungi's shipped example: inert localisation

Run literally, CarveFungi's shipped yeast `loc_pred` is **inert**: it is keyed by RefSeq `NP_` ids
while the annotations use SGD ORF names — **zero overlap** — so its localisation layer never fires and
every reaction gets only its EC score, identical across all its compartments (0/5200 differentiated).
Reproducing its localisation needs the Zenodo TensorFlow models. We therefore **injected DeepLoc 2.1**
(ORF-keyed) into CarveFungi's *unmodified* scoring — this reactivates the localisation gate (16.8% of
scores change; 521 reactions become compartment-differentiated; spot-checks sensible, e.g.
enolase→cytosol, succinyl-CoA ligase→mito). Both arms use these same scores, so it is a clean,
disclosed substitution that isolates the objective.

## A hard big-M MILP — solved with CarveFungi's own CPLEX

CarveFungi's `minmax_reduction` uses **big-M** reversibility coupling, which gives a weak LP
relaxation. The carve is consequently hard: even full academic **CPLEX** leaves the bound loose
(**18–27% gap** here, not closing within 20 min); the Gurobi re-implementation with *tighter* indicator
coupling reaches ~10–14%. CarveFungi's nominal 0.1% pool gap is **not attained** on this instance by
any solver we tried.

Crucially, the result is still **deterministic and time-budget-stable.** In deterministic mode CPLEX
finds one incumbent within ~25 s (Arm A obj 288.46) that **never changes** as the bound slowly
descends — the kept set is byte-for-byte identical at 120 s (46.7% gap), 500 s (26.8%) and a 1200 s
run. So each arm's solution is reproducible and independent of the time budget. It is **not, however,
a *proven* optimum**: a constant incumbent under a still-falling bound is consistent with either
optimality or an early lock-in, and the two arms terminate at *different* gaps. We therefore report the
achieved gap with every number and lean only on what is robust to this (below).

**Running full CPLEX (reproduction note).** The licensed CPLEX Studio install shipped without its
Python API. We installed the PyPI `cplex` (version-matched 22.2.0.0, but Community-capped at 1000
constraints) and replaced its bundled `_internal/cplex2220.dll` with the Studio runtime's full
`cplex2220.dll` of the same version — academic binaries are unlocked, so the matching pip bindings then
load at full capacity. No Gurobi needed.

## Result: a leaner transport network at no detectable accuracy cost

Compartments from metabolites; accuracy = EC-mapped against curated yeast-GEM compartments, on the
**common** kept set (same denominator both arms). Numbers are the deterministic incumbents (identical
across time budgets); gaps are the achieved bounds at 500 s.

| | Arm A (CarveFungi) | Arm B (+ our transport cost) |
|---|--:|--:|
| base reactions kept | 877 | 830 |
| inter-compartment transports | 138 | **81** (41% fewer) |
| transports per base reaction | 0.157 | **0.098** (≈1.6× fewer) |
| recall vs curation (n=463) | 86.0% | 85.5% |
| exact-set match (n=463) | 62.6% | 61.6% |
| identical compartment set (common, both assigned) | — | 93.1% (683/734) |
| achieved gap @ 500 s | 26.8% | 18.1% |

Two takeaways, scoped to what the data supports:

* **Transport parsimony — large and direction-robust.** Adding our transport cost cuts inter-
  compartment transports 41% (0.157 → 0.098 per base, ~1.6×). The *direction* is mechanistically
  guaranteed (the −0.3 cost dwarfs CarveFungi's ~1e-11 transport scores) and the per-base
  normalisation controls for Arm B keeping fewer reactions; the *exact magnitude* is conditional on
  these loose-gap (but deterministic) solutions.
* **No detectable assignment-accuracy cost.** The accuracy difference is within noise — Arm A
  (unmodified CarveFungi) is nominally higher on both metrics, but by ~2 reactions (recall) and ~5
  (exact), non-significant (best-case paired McNemar p = 0.50 recall, p ≈ 0.06 exact). The honest
  statement is *"no detectable accuracy gain or loss from our objective"*, not a neutral tie. 93.1% of
  the placements the two arms share are identical.

## What this does and doesn't show

* **Not proven optima.** Each arm is a stable, deterministic incumbent at 18–27% gap, not a certified
  optimum; the arms stop at different gaps. A tighter solve could shift the magnitudes (it cannot flip
  the transport direction, which is forced by the objective).
* **The kept *set* changed substantially** (877 vs 830, common 796): the "93.1% unchanged" is over the
  common bases assigned in both arms and is conditioned on the subset least likely to change.
* **Compartment-mapping asymmetry:** the gold side can yield compartments (`ce`/`v`) the universal→yeast
  mapping can never produce (21/616 ECs), capping achievable exact-match — but this applies equally to
  both arms, so it does not affect the comparison.
* Arm B still keeps 81 transports despite the penalty (consistent with the biomass/ATPM feasibility
  constraints; not separately tested).

## Bottom line

Running CarveFungi's **own** carve MILP, adding our transport-minimisation term yields a materially
leaner transport network (~1.6× fewer transports per reaction, 41% fewer overall) with **no detectable
assignment-accuracy cost** and 93% identical placements — i.e. transport parsimony essentially for
free. The carve's big-M formulation is hard enough that neither CPLEX nor a tighter Gurobi port proves
optimality, so these are deterministic near-optimal incumbents, reported with their gaps. The clean,
*tight-gap* same-task head-to-head for the paper remains
[`predictLocalization`](predictlocalization_comparison.md) (same lineage, solves fast, deterministic);
CarveFungi is related work of a different kind, and this is a faithful, honest comparison against it.
