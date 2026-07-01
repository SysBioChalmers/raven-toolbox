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
  * [`scripts/analyse_carvefungi_transports.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/scripts/analyse_carvefungi_transports.py)
    — the transport-fidelity investigation below (curated match, functional impact, connectivity).
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

"Free" here means *in compartment-assignment accuracy*. Whether the leaner transport network is also
*biologically better* is a separate question — investigated next.

## Is less better? A transport-fidelity investigation

Arm B drops 66 of Arm A's 138 transports (72 are shared, 9 are B-only), concentrated on
cytosol↔mito (27), cytosol↔extracellular (20, mostly sugar/polyol export), cytosol↔ER (7) and
cytosol↔peroxisome (7). The cargo includes textbook shuttles — citrate, (S)-malate, oxaloacetate,
2-oxoglutarate (the malate–aspartate and citrate shuttles), plus trehalose/fructose/mannose export.
Three checks ask whether dropping them is an improvement. Driver:
[`scripts/analyse_carvefungi_transports.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/scripts/analyse_carvefungi_transports.py).

**1. The reduction is *not selective.*** Match each carved transport (metabolite + compartment pair)
against the curated, literature-backed yeast-GEM transportome. If the cut were "smart", dropped
transports would match curation *less* than kept ones. They don't:

| transport set | matches a curated yeast-GEM transport |
|---|--:|
| shared (both arms keep) | 39% |
| **dropped by Arm B** | **42%** |
| kept by Arm A (all) | 41% |

Arm B sheds real and spurious transports at the same rate (classifiable n = 46/59/105; these are
qualitative rates, not a significance test). This is expected: the carve has **no transporter-level
evidence** (transport scores are ~1e-11), so a blanket −0.3 penalty just removes whatever the network
can do without while keeping biomass feasible.

**2. The dropped transports are functionally load-bearing.** Map the dropped transports to their
curated counterparts (31 yeast-GEM reactions) and delete them from yeast-GEM — a *proper* model with
GPRs and validated growth (baseline 0.081). Growth collapses to **0**, and **5 are individually
essential**: 2-oxoadipate/2-oxoglutarate, 2-dehydropantoate (a CoA precursor), NADP⁺, NADPH, and
serine transport. The *shared* (kept-by-both) transports are equally load-bearing — their 24 curated
counterparts also collapse growth when removed, with 3 individually essential — so essential transports
are spread across kept and dropped alike. That is precisely the **indiscriminate** point: Arm B's cut
is not concentrated on the dispensable ones. This measures importance in *curated* biology — Arm B
itself stays feasible (biomass ≥ 0.1 by construction; the carve routes around the cuts), so the finding
is that its network *diverges* from curated yeast, not that it fails to grow. (Gene essentiality cannot
be tested on the carved models themselves — the universal DB has 0 GPRs — which is why this uses
yeast-GEM.)

**3. Connectivity barely changes — the carve re-routes through exchanges.** Structurally (internal
network, exchanges excluded; the carve guarantees *flux*-connectivity *with* exchanges by its ε-flux
coupling, so this exposes the latent gaps exchanges otherwise hide):

| | Arm A | Arm B |
|---|--:|--:|
| connected components | 2 (giant + a 3-rxn island) | 2 (giant + the same island) |
| dead-end metabolites | 186 (17.3%) | 194 (18.6%) |

Dropping transports does **not** fragment the network into isolated sub-networks. It does strand
modestly more metabolites: 19 are mass-balanced in A but dead-end in B (vs 8 the other way), localised
to cytosol (8), mito (7), peroxisome (3), ER (1) — exactly the dropped-transport cargo (2-dehydro-
pantoate, formate, butyrate, peroxisomal citrate/ammonium…). These then lean on boundary exchanges
(secrete/import) rather than internal transport to stay balanced — feasible, but less biologically
self-contained.

**Verdict: less is more *parsimonious*, not more *correct*.** The reduction is indiscriminate, removes
functionally essential curated transports, and modestly raises dead-ends — though it does not break
global connectivity, because the carve re-routes via the (artificial) environment. The root cause is
that the transport penalty is a blanket prior applied **without transporter-level evidence.**

## Toward evidence-aware transport scoring

The fix follows directly: make the transport cost *evidence-aware*, so the reduction becomes selective
— penalise transports with **no** transporter support while retaining those with sequence-level
evidence. This mirrors how the localisation module already scores *metabolic* reactions by gene
localisation; it simply extends the same predictor-agnostic, sequence-derived evidence to transport.

**Evidence sources — every carrier, every membrane.** All four are sequence-, HMM-, or
orthology-derived, so they cover *any* transporter family across *any* membrane (the dropped transports
span c↔mito 27, c↔extracellular 20, c↔ER 7, c↔peroxisome 7 — the scoring must not privilege one
membrane):

* **Transporter family (Pfam / hmmer) — the backbone.** One `hmmscan` against the transporter clans
  flags carrier genes of all families: the mitochondrial carrier family (MCF, `PF00153`/SLC25 — the
  c↔mito carriers Arm B dropped), major facilitator (MFS), ABC, amino-acid/sugar permeases,
  aquaporins, P-type ATPases, and so on. Family identity also gives a coarse substrate class.
* **Transporter classification (TCDB, via DIAMOND).** A `diamond blastp` against TCDB assigns a TC
  number → substrate class **and** mechanism (uni/sym/antiport): the substrate-specific gold standard.
* **Compartment placement (DeepLoc, already in-pipeline).** *Which* membrane a carrier sits on follows
  from the gene's predicted **compartment** — the reliable organelle outputs (trust 0.78–0.88), *not*
  the noisy membrane-*type* output (`mm` ≈ 0.86 but `erm`/`gm`/`vm` ≈ 0). A carrier-family gene
  predicted in compartment *X* supports transports across *X*'s boundary; this generalises to every
  compartment.
* **Orthology (already available).** EggNOG/KEGG orthogroups flag transporter orthologs with
  substrate/direction.

**Scoring.** Replace the constant `transport_cost` with a per-transport cost. For a candidate transport
*t* moving metabolite *m* across membrane *M* = {c₁,c₂}:

```
evidence(t)       = max over genes g of  conf_transporter(g) · compartment_match(g, M) · substrate_match(g, m)
transport_cost(t) = base_cost · (1 − evidence(t))      # supported → cheap; unsupported → full prior
```

`conf_transporter` from the Pfam/TCDB hit strength, `compartment_match` from the DeepLoc compartment,
`substrate_match` from the TC/family substrate class vs *m*'s class. It drops straight into the
assignment MILP objective (the transport term becomes per-reaction), symmetric with the existing
per-gene localisation scoring, and recovers today's constant −0.3 when `evidence = 0`.

**Organism-agnostic by design.** None of the evidence is a species-specific transporter table — they
are universal HMMs (Pfam), a cross-organism sequence DB (TCDB), cross-species orthogroups (EggNOG/KEGG),
and a eukaryote-wide predictor (DeepLoc). The only per-organism input is the **proteome FASTA**; the
compartment set comes from the target model (the module already maps cross-kingdom compartments, e.g.
plastid for plants). So the same pipeline runs unchanged on any eukaryote — consistent with the
module's cross-kingdom DeepLoc validation (yeast, *Arabidopsis*, *Chlamydomonas*, human). yeast-GEM is
only the *benchmark* here, not a dependency.

**Phased implementation** (by evidence maturity — all-carrier from the start, not membrane-by-membrane):

1. **Family scan** (Pfam/hmmer) over all transporter clans → per-gene "is a carrier" + coarse substrate
   + compartment placement from DeepLoc. Covers every membrane immediately.
2. **TCDB** (DIAMOND) → TC-number substrate specificity + mechanism → substrate-matched scoring.
3. **Consensus/refinement** — combine family + TCDB + orthology + DeepLoc, add transport
   directionality, resolve conflicts.

Mitochondria are merely the cleanest *validation* exemplar (the one membrane where DeepLoc's
membrane-type output independently corroborates, and whose curated essential carriers — malate/2-OG/
citrate — are textbook); the scoring itself is membrane- and organism-agnostic.

**Validation.** Reuse this study's benchmark — curated transport precision/recall plus the functional
(essentiality) test — before vs after. Success criteria: the kept-transport curated-match rate rises
*above* the dropped rate (the cut becomes selective), the 5 individually-essential transports are
retained, and the gains reproduce on a **non-fungal** model (e.g. AraCore) to confirm
organism-agnosticism.

**Design caveats.** Substrate matching (metabolite → substrate class) is the hard part; start coarse
(sugars / amino acids / organic acids / ions / nucleotides / lipids). Absence of evidence ≠ absence of
a transporter, and annotation completeness varies by organism — so keep a *mild, tunable* prior on
unsupported transports; never hard-forbid.

## Result: evidence-aware scoring makes the cut selective

The evidence-aware cost above is now implemented (`raven_toolbox.localization.transport_evidence`; see
[the reference](../reference/transport_evidence_scoring.md)) and scored against **this carve**:
`analyse_carvefungi_transports.py` annotates the yeast proteome (364 transporter genes via
`hmmsearch` + `diamond`), builds the per-metabolite `transport_cost`, and scores every approach on the
*same* candidate set (Arm A's 138 carved transports) against the curated yeast-GEM transportome (43
curated, 9 individually essential):

| approach | kept | curated replicated | essential kept | spurious kept |
|---|--:|--:|--:|--:|
| CarveFungi (no transport penalty) | 138 | 43/43 | 9/9 | 95 |
| CarveFungi (blanket −0.3) | 72 | 18/43 | **3/9** | 54 |
| ours: coarse | 98 | 41/43 | **9/9** | 57 |
| ours: + ChEBI | 99 | 41/43 | **9/9** | 58 |
| ours: + ChEBI + sibling 0.5 | 107 | **43/43** | **9/9** | 64 |

The blanket penalty is leanest but **indiscriminate** — it drops 6 of 9 essential carriers (incl.
2-oxoglutarate, 2-dehydropantoate, NADP⁺/NADPH, serine) and more than half the curated transports.
Every evidence-aware variant retains **all 9 essential** and 41–43 of 43 curated at moderate parsimony
(98–107 kept). "ours" isolates the transport-cost objective (keep the network-needed transports plus
those the evidence supports); it is not a re-solve of the carve MILP. The non-fungal (AraCore)
reproduction remains the outstanding organism-agnosticism check.

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
assignment-accuracy cost** and 93% identical placements. But "leaner" is not automatically "better":
the cut is **indiscriminate** — it drops curated, functionally essential transports (≥5 individually
essential in yeast-GEM) at the same rate as spurious ones, because the carve has no transporter-level
evidence. The actionable conclusion is the **evidence-aware transport scoring** proposed above:
penalise only *unsupported* transports. The carve's big-M formulation is hard enough that neither CPLEX
nor a tighter Gurobi port proves optimality, so these are deterministic near-optimal incumbents,
reported with their gaps. The clean, *tight-gap* same-task head-to-head for the paper remains
[`predictLocalization`](predictlocalization_comparison.md) (same lineage, solves fast, deterministic);
CarveFungi is related work of a different kind, and this is a faithful, honest comparison against it.
