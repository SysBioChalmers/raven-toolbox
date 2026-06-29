# CarveFungi assignment head-to-head: faithful harness, tractability, and the robust finding

A genuine empirical comparison of CarveFungi's compartment-**assignment** objective against ours, on
CarveFungi's own intermediate state. This fixes the three flaws that made an earlier quick emulation a
strawman (see [carvefungi_analysis.md](carvefungi_analysis.md)) — and is honest about what the result
can and cannot support.

* Driver: [`scripts/benchmark_carvefungi_milp.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/scripts/benchmark_carvefungi_milp.py)
* CarveFungi: [github.com/SandraCastilloPriego/CarveFungi](https://github.com/SandraCastilloPriego/CarveFungi),
  bioRxiv [2023.08.23.554328](https://doi.org/10.1101/2023.08.23.554328).

## Design — the intermediate-state swap

Hold CarveFungi's **real intermediate state** fixed and swap only the assignment objective:

* **Candidate set:** CarveFungi's universal model `bigModelv2.21b.sbml` — each reaction exists only in
  the compartments the DB instantiates it in (mean 2.16, max 4; mostly cytosol/mito/peroxisome/ER).
* **Scores:** CarveFungi's *unmodified* scoring code (`EggNogScoring`) on the shipped *S. cerevisiae*
  EggNOG annotations, run to produce the exact per-(reaction, compartment) score dict it feeds its
  carve MILP.
* **Both arms** enforce `S·v=0` + ε-flux indicator coupling (a reaction is on only if it carries ≥ε
  flux) + hard biomass ≥ 0.1 and ATP-maintenance ≥ 0.1, on the *same* scores. Only the objective's
  parsimony terms differ:
  * **Arm A (CarveFungi):** `max Σ score(r,c)·y[r,c]` — no transport cost, no multi-localisation penalty.
  * **Arm B (ours):** the same, plus a per-transport cost (our transport-minimisation term).

A Gurobi port of CarveFungi's `minmax_reduction` runs both arms (CPLEX is unavailable in this Python
3.14 env). An independent adversarial review confirmed the port is **faithful — not a strawman**: it
respects the candidate set, the ε-flux connectivity, CarveFungi's score signs (−3 unannotated default,
near-zero transports), exchange/uptake handling, and biomass/ATPM, and feeds both arms identical
scores so the objective is genuinely isolated.

## A repro bug in CarveFungi's shipped example: inert localisation

Run literally, CarveFungi's shipped yeast `loc_pred` is **inert**: it is keyed by RefSeq `NP_` ids
while the EggNOG annotations use SGD ORF names — **zero overlap** — so its localisation layer never
fires and every reaction gets only its EC score, identical across all its compartments (0/5200
differentiated). Reproducing its localisation needs the Zenodo TensorFlow models. We therefore
**injected DeepLoc 2.1** (ORF-keyed) into CarveFungi's *unmodified* `scoring_compartments` — this
reactivates the localisation gate (16.8% of scores change; 521 reactions become compartment-
differentiated; spot-checks sensible, e.g. enolase→cytosol, succinyl-CoA ligase→mito). Both arms use
these same scores, so it is a clean, disclosed substitution that isolates the objective.

## What we can conclude — and what we can't

**The carve is a hard MILP.** With Gurobi's tight indicator coupling, Arm A reaches only ~11%
optimality gap in 600 s and Arm B ~14% (the multi-localisation-penalty variant is far worse, ~170%+,
essentially unsolved). CarveFungi itself runs CPLEX to a 0.1% pool gap. At ~10–14% gaps the
**agreement-with-curation comparison is gap-sensitive and unstable across runs** (it flips with the
solve and with the multi-loc parameter), so **we do not report it as a finding.** A definitive
accuracy head-to-head needs a tight-gap solve — i.e. CPLEX (below).

**The one gap-robust finding — transport parsimony.** Our transport cost cuts the inter-compartment
transport reactions sharply, and this survives the loose gap because the −0.3 cost dwarfs CarveFungi's
~1e-11 transport scores:

| | Arm A (CarveFungi) | Arm B (+ our transport cost) |
|---|--:|--:|
| transports kept | 401 | 75 |
| transports per reaction kept | 0.24 | **0.053** (≈4.5× fewer) |

The per-reaction rate (not just the total) drops ~4.5×, so it is a real objective effect, not merely
a consequence of Arm B keeping fewer reactions. Mechanistically this is our transport-minimisation
doing exactly what it is designed to: a less cross-compartment-shuttling network. Whether that comes
at any cost to assignment accuracy is the gap-sensitive question CPLEX must settle.

## Getting the definitive numbers (CPLEX)

The trustworthy accuracy head-to-head needs CarveFungi's MILP solved to a tight gap, which means
running its own `minmax_reduction` (CPLEX) in a **Python 3.11/3.12** env where `import cplex` works
(CPLEX 22.2 has no Python 3.14 binding; this repo's env is 3.14). The recipe, all reproducible:

1. Produce the DeepLoc-injected score dict with CarveFungi's unmodified scoring (as above).
2. Run `minmax_reduction(universal, scores)` → **Arm A**.
3. Run it again with each transport reaction's score reduced by the transport cost → **Arm B**.
4. Score both on the **common** kept set against EC-mapped curated yeast-GEM compartments
   (`gold_reference`/`ec_eval` in the driver), to a tight, equal gap.

`scripts/benchmark_carvefungi_milp.py` implements the comparison + the EC-mapped gold reference; only
the solver needs to be CPLEX for the accuracy claim to be trustworthy.

## Bottom line

The harness is faithful and the methods' mechanistic difference is real, but a *tight-gap* empirical
accuracy comparison against CarveFungi requires CPLEX and is left as a reproducible recipe. The clean,
same-task empirical head-to-head for the paper remains
[`predictLocalization`](predictlocalization_comparison.md) (same lineage, solves fast, deterministic);
CarveFungi is related work of a different kind, and against its own pipeline our transport-minimisation
robustly yields a ~4.5×-more-parsimonious transport network.
