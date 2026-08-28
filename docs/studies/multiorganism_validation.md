# Multi-organism validation (P1): certified assignment across four kingdoms

The generality test for `assign_compartments` — does the method work beyond yeast, at genome
scale, and on the **chloroplast** (a compartment CarveFungi's fungal 4-category scheme cannot
represent). Answer: **yes on all counts.** Four organisms spanning fungi, mammals, plants, and green
algae, all placed from DeepLoc-2 evidence and FBA-certified, with **no per-organism method changes** —
only a compartment-label map and a metabolite base-key per organism.

*Scripts:* `scripts/benchmark_certified_yeast.py` (yeast, with the CarveFungi head-to-head) and
`scripts/benchmark_certified_multiorg.py` (the other three). Results in
`.research_tmp/certified_{yeast,multiorg}.json`.

## Results

| organism | kingdom | model | reactions | compartments | chloroplast | wall | certified | reaction agr | gene agr | transports | blocked |
|---|---|---|--:|--:|:--:|--:|:--:|--:|--:|--:|--:|
| *S. cerevisiae* | fungal | yeast-GEM | 2569 | 14 | no | ~140 s | **yes** | 72.5 % | **88.7 %** | 967 | 29.7 % |
| *H. sapiens* | mammalian | Human-GEM | 12877 | 9 | no | 190–464 s | **yes** | 56.9 % | 79.4 % | 3608 | 10.7 % |
| *A. thaliana* | plant | AraCore v2.1 | 585 | 6 | **yes** | 35 s | **yes** | **82.1 %** | **93.7 %** | 200 | 0.0 % |
| *C. reinhardtii* | algal | iCre1355 | 2394 | 11 | **yes** | 93 s | **yes** | 54.1 % | 70.5 % | 1083 | 32.9 % |

All rows are measured on the **deterministic score-aligned placement** (each reaction placed by its own
enzymes' DeepLoc scores; see [yeast_validation.md](yeast_validation.md)), not the earlier arbitrary
co-optimum. **Gene** agreement is unchanged from the arbitrary-placement runs — the tie-break fixes only
the reaction placement, leaving the gene layout untouched — which is why every gene-agreement figure
above matches its pre-tie-break value. Reaction agreement rose or held (Human-GEM 52.8 % → 56.9 % the
largest move; AraCore and iCre1355 were already near their ceilings), and coherent placement needs fewer
transports across the board. The one non-obvious shift is iCre1355's blocked fraction (19.0 % → 32.9 %):
tighter, less transport-heavy placement leaves more reactions dead-ended, though the model still
certifies. Placement is reproducible run-to-run.

Agreement is scored against each model's **own** compartment annotation (a *circular* truth — the model
being re-placed is also the reference). This measures self-consistency and functionality, not accuracy;
the independent-ground-truth work (P2: HPA for human, SUBA5 for Arabidopsis, the Chlamydomonas
Chloroplast Protein Atlas, UniProt experimental across all) replaces it for the accuracy claim.

## What P1 establishes

1. **It scales past 4 compartments at genome size.** Human-GEM — 12877 reactions, **9 compartments**,
   5221 reactions re-placed — certifies (materialised biomass above the floor by a real FBA) in
   ~3–8 min. iCre1355 reaches **11 compartments** (the richest organelle set here) in ~1.5 min. The
   flux-free master + FBA-certification design does not degrade as compartments and reactions grow;
   this was the open scaling risk, now closed.
2. **It handles the chloroplast** — the generality hook CarveFungi structurally cannot represent.
   DeepLoc's "Plastid" label maps to the chloroplast code (`h` in both plant models), and both AraCore
   (plant) and iCre1355 (alga) place and certify with a chloroplast, AraCore at **82.1 % reaction /
   93.7 % gene** self-agreement.
3. **It is organism- and predictor-agnostic in practice, not just in principle.** The same function,
   the same DeepLoc-2 evidence, no algorithm changes — only a `{DeepLoc-label -> compartment-code}` map
   and a metabolite base-key (name-based for yeast/Human-GEM; bracket-stripped `Glc[c]->Glc` for
   AraCore; suffix-stripped `10fthf_c->10fthf` for iCre1355) — carry it from a 585-reaction plant core
   model to a 12877-reaction human GEM.
4. **Every model stays functional.** All four materialise a growing model (blocked-reaction fractions
   0–33 %, mostly the models' intrinsic dead-ends), the property that distinguishes this method.

## Honest caveats

- **Reaction agreement drops on the larger models** (52–53 % for Human-GEM and iCre1355 vs 72–82 % for
  yeast/AraCore). Part is the circular truth (more compartments = more ways to disagree at the reaction
  level even when the gene call is right — gene agreement stays 70–94 %), part is DeepLoc coverage on
  non-model-organism proteomes, and part may be genuine over-/under-placement to be examined once
  independent truth is in. Do not present these as accuracy until P2.
- **Runtime varies** (Human-GEM 190 s vs 464 s across two runs) with machine load and the master's
  non-determinism; report a range, not a point.
- **Organelle-genome-encoded proteins** (chloroplast/mito-encoded, e.g. iCre1355 Chre*Cp*) have no
  DeepLoc score and bias organelle recall — quantify in P2.
- **Model coverage differs**: AraCore is a curated *core* model (6 central-metabolism compartments, no
  ER/Golgi/vacuole), so its high agreement is on an easier, smaller problem; Human-GEM/iCre1355 are the
  genome-scale stress tests.

## Next (P2)

Replace the circular per-model truth with the independent experimental resources already sourced (HPA
Enhanced+Supported; SUBA5; Chloroplast Protein Atlas; UniProt ECO:0000269), add the RAVEN
`predictLocalization` and DeepLoc-argmax baselines, and apply McNemar/CI per organism — the accuracy
claim across kingdoms rests on that, not on the self-agreement numbers above.
