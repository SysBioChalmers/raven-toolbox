# DeepLoc normalisation benchmark (whole yeast-GEM)

Does keeping DeepLoc's **raw** probabilities (`load_deeploc(normalise=False)`) instead of rescaling
each gene's best compartment to 1.0 change a compartment assignment's agreement with curated
yeast-GEM? Run on the **entire** model. Short answer: **no — it is accuracy-neutral**, so
normalisation stays the default and `normalise=False` is an opt-in.

* Driver: [`scripts/benchmark_deeploc_normalisation.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/scripts/benchmark_deeploc_normalisation.py)
* Model: `yeast-GEM.xml` — 4102 reactions, 1143 genes, 14 compartments (perfect DeepLoc-gene
  coverage: the input FASTA was prepared from this model).
* Scores: the three committed yeast DeepLoc 2.1 CSVs (`data/deeploc/yeast-GEM_deeploc_00{1,2,3}.csv`),
  loaded twice from the *same* files — `normalise=True` (top→1.0, the default) vs `normalise=False`
  (raw probabilities) — and swept over `transport_cost`. `multi_compartment_penalty=0.5`,
  `mip_gap=0.01`, Gurobi.

## Method

The established flatten benchmark (see the
[yeast localization benchmark](yeast_localization_benchmark.md)), on the full model: flatten
yeast-GEM to one compartment with
[`merge_compartments`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/manipulation/compartments.py)
so the MILP cannot lean on metabolite topology, then ask
[`predict_localization`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/src/raven_toolbox/localization/predict.py)
to re-place every non-boundary single-compartment GPR'd reaction (**truth set: 2207 reactions**)
from the DeepLoc score table.

The MILP objective is `max Σ score[g,c]·y[g,c] − transport_cost·transports −
multi_compartment_penalty·extra_compartments`. The score enters **linearly**, so normalisation makes
every gene vote with weight 1.0, while raw lets a confident gene outvote a shaky one and — because
`transport_cost` and the penalty are **absolute constants** — raises the bar a reaction must clear to
leave the default compartment or a gene to occupy a second one. Normalisation is a monotone per-row
rescale, so **both arms agree on each gene's preferred (argmax) compartment**; they can only diverge
on placement.

Two properties bound what this benchmark can see, and both were confirmed by an independent
methodology review:

1. **DeepLoc addresses only 9 of 14 compartments.** The membrane sub-compartments `erm/gm/mm/vm` and
   the lipid particle `lp` have no DeepLoc column, so the 902 reactions truly there score 0 for
   **both** arms (they cancel in any cross-arm comparison). `addressable` accuracy restricts to the
   1305 reactions whose curated compartment is DeepLoc-reachable.
2. **The `c` majority + a transport-free default give a high floor.** `c` is 698/2207 of the truth
   set and is transport-free, so a trivial *leave everything in `c`* baseline already scores
   **overall 0.316 / addressable 0.535**. Read the numbers below against that, not against zero.

### Reproducibility floor (read this before the tables)

`predict_localization`'s objective has **no term in the reaction-placement variable** — when a gene
sits in several compartments the MILP does not pin *which* of them its reaction lands in. That tie is
broken by the solver (sensitive to variable/set ordering and `mip_gap`), so:

* **Gene-level** quantities (multi-localisation counts) are **exactly reproducible**.
* **Reaction-level** `overall`/`addressable` are stable to ~1pp run-to-run.
* **`macro`** (the unweighted mean of per-compartment accuracies) **is not** — it equal-weights
  compartments with n = 7–13, so a handful of tie-broken reactions swing it by ~0.1.

Two independent runs of the *identical* configuration (`transport_cost=0.01`) make this concrete:

| run | arm | overall | addressable | macro | multi-loc genes |
|---|---|---:|---:|---:|---:|
| 1 | normalised | 0.396 | 0.670 | **0.362** | 617 |
| 2 | normalised | 0.402 | 0.680 | **0.475** | 617 |
| 1 | raw | 0.399 | 0.675 | 0.472 | 392 |
| 2 | raw | 0.399 | 0.674 | 0.470 | 392 |

`multi-loc genes` is identical to the unit; `overall`/`addressable` move ≤1pp; normalised `macro`
swings **0.113** between runs — far larger than any cross-arm difference below. **Treat sub-1pp
accuracy deltas and all `macro` differences as within noise.**

## Accuracy vs. transport_cost

`moved` = reactions placed outside the default `c`. `*` = best operating point with
`transport_cost > 0` (by `addressable`); `(deg.)` marks `transport_cost = 0`, which scores highest
but is degenerate (no transport term ⇒ multi-compartment-gene reactions are an arbitrary tie-break),
so the sections below use the matched, well-conditioned `transport_cost = 0.01`.

| arm | transport_cost | overall | addressable | macro | moved | multi-loc genes |
|---|---:|---:|---:|---:|---:|---:|
| normalised (deg.) | 0.0 | 0.418 | 0.707 | 0.474 | 1271 | 610 |
| normalised \* | 0.01 | 0.402 | 0.680 | 0.475 | 756 | 617 |
| normalised | 0.02 | 0.391 | 0.661 | 0.339 | 584 | 625 |
| normalised | 0.05 | 0.377 | 0.637 | 0.308 | 221 | 651 |
| normalised | 0.1 | 0.372 | 0.628 | 0.264 | 147 | 674 |
| raw (deg.) | 0.0 | 0.420 | 0.711 | 0.500 | 1345 | 376 |
| raw \* | 0.01 | 0.399 | 0.674 | 0.470 | 828 | 392 |
| raw | 0.02 | 0.389 | 0.657 | 0.448 | 671 | 401 |
| raw | 0.05 | 0.388 | 0.657 | 0.434 | 278 | 434 |
| raw | 0.1 | 0.374 | 0.633 | 0.279 | 153 | 476 |

## Head-to-head at matched transport_cost = 0.01

Only the scores differ (same model, truth set, knobs, operating point).

| metric | normalised | raw | delta (raw − norm) |
|---|---:|---:|---:|
| overall | 0.402 | 0.399 | −0.004 |
| addressable | 0.680 | 0.674 | −0.006 |
| macro | 0.475 | 0.470 | −0.004 |
| multi-loc genes | 617 | 392 | −225 |

Against the *leave-everything-in-`c`* baseline (overall 0.316 / addressable 0.535), both arms add
the same ~9 / ~14 pp. Every accuracy delta is within the reproducibility floor — **the arms are
indistinguishable on agreement.** The only difference that exceeds noise is structural and at the
gene level: **raw assigns far fewer genes to multiple compartments** (392 vs 617), because its
calibrated magnitudes (~0.5–0.9) clear the absolute `multi_compartment_penalty` less often than the
normalised 1.0. The same effect is reachable by *raising* `multi_compartment_penalty` on normalised
scores — i.e. raw is a re-scaling of the existing knobs, not new information.

## Where the arms actually differ

**Contested subset.** Only 157 of 2207 reactions are placed differently at `transport_cost = 0.01`;
the other 2050 are identical by construction (shared argmax). On the contested ones **both arms place
poorly and raw is no better**: normalised 29/157 = 0.185, raw 21/157 = 0.134. (This subset is itself
tie-break-sensitive; the takeaway is "both bad", not the exact split.)

**Confidence-stratified (single-gene reactions, binned by the gene's raw DeepLoc top probability).**
If raw probabilities helped, the gain would concentrate in the high-confidence bins. It does not —
raw is flat-to-slightly-worse everywhere, including the high-confidence tail:

| raw top prob | n | normalised | raw | delta |
|---|---:|---:|---:|---:|
| [0.0, 0.5] | 23 | 0.217 | 0.217 | +0.000 |
| [0.5, 0.7] | 633 | 0.248 | 0.235 | −0.013 |
| [0.7, 0.9] | 792 | 0.451 | 0.451 | +0.000 |
| [0.9, 1.0] | 240 | 0.438 | 0.404 | −0.033 |

**Per-compartment (matched `transport_cost = 0.01`).** Differences sit inside the noise floor; the
big addressable compartments (`m` n=219, `p` n=115) favour normalised or tie, and the tiny ones
(`e`, `g`, `v`) are equal here — the `g` swing that inflated `macro` in an earlier run did not
reproduce.

| compartment | n | normalised | raw |
|---|---:|---:|---:|
| c | 698 | 0.926 | 0.920 |
| ce | 110 | 0.000 | 0.000 |
| e | 9 | 0.778 | 0.778 |
| er | 90 | 0.611 | 0.622 |
| erm (unreachable) | 348 | 0.000 | 0.000 |
| g | 13 | 1.000 | 1.000 |
| gm (unreachable) | 52 | 0.000 | 0.000 |
| lp (unreachable) | 146 | 0.000 | 0.000 |
| m | 219 | 0.594 | 0.594 |
| mm (unreachable) | 296 | 0.000 | 0.000 |
| n | 44 | 0.068 | 0.068 |
| p | 115 | 0.296 | 0.252 |
| v | 7 | 0.000 | 0.000 |
| vm (unreachable) | 60 | 0.000 | 0.000 |

## Conclusion

* **Per-gene normalisation is accuracy-neutral for assignment.** On the whole model, normalised and
  raw agree with curated yeast-GEM to within the reproducibility floor on every metric. This is
  partly structural — the arms share each gene's argmax, so 2050/2207 reactions are placed
  identically and they can only diverge on 157 contested reactions, where both do badly.
* **Confidence-gating does not help.** Raw does not rescue the high-confidence calls (it is slightly
  *worse* on the [0.9, 1.0] bin) nor the contested ones. The hypothesis that keeping calibrated
  magnitudes would let confident genes be placed better is **not supported** here — the most robust
  finding in the study (large n, conservative direction).
* **The one reproducible difference is a knob re-scaling.** Raw's fewer multi-localisations (392 vs
  617) follow mechanically from comparing calibrated probabilities against the *absolute*
  `multi_compartment_penalty`; the same is obtainable by tuning that penalty (and `transport_cost`)
  on normalised scores — the calibration lesson already in the
  [yeast localization benchmark](yeast_localization_benchmark.md).
* **Scope.** This isolates the *scoring* effect: topology is flattened and there are no
  flux/functionality constraints, so the absolute accuracies are not expected real-world performance,
  only a fair within-benchmark comparison of the two score scalings.

**Decision:** keep per-gene normalisation the **default** (`parseScores` convention, comparable
multi-source scales); expose `load_deeploc(normalise=False)` as an **opt-in** for callers who want
DeepLoc's calibrated magnitudes — e.g. the confidence signal
[`triage_localization`](../guide/localization.md) consumes, or to avoid hand-tuning the transport
scale. Flipping the default is not justified by this benchmark.

## Reproducing

```bash
python scripts/benchmark_deeploc_normalisation.py \
    --yeast-gem /path/to/yeast-GEM/model/yeast-GEM.xml \
    --doc /tmp/normbench.md      # scratch — this study page is curated by hand
```

The script regenerates every table above (the head-to-head, contested, confidence and
per-compartment sections at the matched `--operating-point`, default 0.01). This page is curated, so
point `--doc` at a scratch file rather than overwriting it. A single run cannot reproduce the
cross-run reproducibility table either — that compares two independent runs to expose the `macro`
tie-break noise; numbers in the other tables are from one representative run and vary within the
floor described above.
