# CarveFungi: how it works, and how our compartment-assignment MILP differs

[CarveFungi](https://github.com/SandraCastilloPriego/CarveFungi) (Castillo-Priego et al., bioRxiv
[2023.08.23.554328](https://doi.org/10.1101/2023.08.23.554328)) is the only recent method besides
RAVEN's `predictLocalization` that does network-aware, predictor-driven, optimization-based
compartmentalisation of eukaryotic metabolic models. This note records how it actually works (from a
source read of the repo) and the mechanistic differences from our assignment MILP — for the
methods-comparison section of a paper.

## What it is

An **end-to-end fungal GEM reconstruction** pipeline: EggNOG annotation → a 20-model deep-learning
protein-localisation ensemble → per-(reaction, compartment) scores → a **carveMe-derived MILP** (IBM
CPLEX) that *carves* a fungal universal model into a functional, compartmentalised model. It is not a
stand-alone compartment assigner — localisation is one signal inside a reconstruction-and-carving
optimisation.

## The MILP (the core to compare against)

Solved with CPLEX over a *single, already-compartmentalised universal model* (`bigModelv2.21b.sbml`,
8 fungal compartments: c, m, x/peroxisome, r/ER, l/lipid particle, n, g, e). Each (reaction,
compartment) pair is a distinct universal reaction with its own binary.

* **Variables:** continuous flux per universal reaction (objective coeff 0); binary direction
  indicators `yf_r`/`yr_r` per reaction (objective coeff = the reaction's score); binary uptake
  indicators for exchanges.
* **Objective:** `max Σ_r score(r)·(yf_r + yr_r) + Σ uptake_score·y_e` — a score-weighted **reaction
  selection**. Compartmentalisation is implicit in *which compartment-specific copy* is kept.
* **Constraints:** steady-state `S·v = 0`; big-M/ε indicator coupling so a reaction is "on" iff it
  carries ≥ ε flux (the connectivity guarantee — no dead reactions); **hard biomass ≥ min_growth and
  ATP-maintenance ≥ min_atpm** (every output model grows by construction).
* **No transport-cost term and no multi-localisation penalty.** Transport reactions are ordinary
  compartmental reactions given a near-zero score (1e-11), so the optimiser keeps a transporter only
  when connectivity forces it — it does not *reward* transport, but it does not *penalise* it either,
  and it never counts transported metabolites. The same enzymatic reaction may be switched on in
  several compartments simultaneously with **zero** penalty and no coupling between a gene's copies.

## Scoring (localisation → reaction score)

Two-stage. (A) EggNOG EC evidence → a per-reaction confidence in `[0, 3]` (unannotated reactions get
−3; transport/spontaneous get ~0). (B) a localisation **gate**: the deep predictor emits 4 classes
(E=ER/secretory, M=mito, P=peroxisome, O=other) mapped onto the 8 compartments with per-compartment
cutoffs (peroxisome 0.10, mito/cytosol 0.45, ER 0.20). The fused per-(reaction, compartment) score is
roughly `(SCALE − ec_score) · (P(class) − cutoff)`: the **sign** comes from whether the protein clears
the compartment's cutoff (below ⇒ negative ⇒ that compartment is penalised), and the **magnitude** is
*larger when EC evidence is weak* (localisation only decides placement when annotation is ambiguous;
strong EC evidence dominates). A safeguard reverts to the plain EC score if every compartmental copy
would go negative, so a confident enzyme is never deleted outright.

## Determinism

**Ensemble by design, not a single optimum.** It uses CPLEX's solution pool (`relgap = 0.001`) to
enumerate near-optimal solutions, keeps the top 5, runs the whole carve **twice** (constrained vs.
open medium), and crosses the two sets into one ensemble model whose reaction *names* encode in how
many members each reaction appears. The tutorial instructs the user to "create few models and select
the one with bigger objective." There is no single canonical output.

## The one mechanistic difference that matters

**CarveFungi maximises annotation-and-localisation-weighted reaction *selection* from a universal
network, subject to mass balance + forced biomass (a "carve a universal model" objective). Our method
minimises a localisation-disagreement + transport-cost + multi-localisation-penalty *placement*
objective over a fixed reaction set (an "assign with transport minimisation" objective).** CarveFungi
decides *whether a reaction exists and grows*; ours decides *where existing reactions live, as
parsimoniously as the transport budget allows*. Placement parsimony is simply not in CarveFungi's
objective — only network functionality is.

## What we should borrow

* **Functionality as a first-class MILP constraint.** CarveFungi makes biomass + ATPM hard
  constraints with full `S·v=0` + ε-flux indicator coupling *inside* the selection MILP, guaranteeing
  every kept reaction carries flux. Our ε-flux/biomass mode (`assign_compartments`) is optional;
  CarveFungi is a good template for a "guaranteed-functional" assignment mode.
* **Confidence-weighted localisation.** Down-weighting the localisation signal when annotation
  confidence is high (and up-weighting it when ambiguous) is a principled idea our flat score scale
  lacks.
* **Per-compartment cutoffs** (predictor calibration differs by organelle — cf. our own per-
  compartment trust table) and the **"don't delete a confident enzyme"** safeguard.
* **Solution-pool enumeration** as an optional *diagnostic* of alternative-optima sensitivity (while
  keeping our single deterministic optimum as the headline output).

## Where our approach is clearly more general

* **Predictor-agnostic.** CarveFungi hard-wires its 4-class ensemble and the 8→4 compartment map;
  swapping in DeepLoc/UniProt/experimental data means rewriting the scorer. Ours takes any
  gene × compartment score table.
* **A real transport-cost objective and a multi-localisation penalty.** CarveFungi has neither, so it
  cannot discourage biologically implausible compartment fragmentation or parsimoniously limit
  multi-localisation; ours does both explicitly.
* **Determinism.** One reproducible global optimum vs. an intrinsic ensemble — better for curation and
  reproducibility.
* **Organism generality.** CarveFungi is fungal-locked (curated fungal universal DB, hand-coded fungal
  biomass, genus-specific tweaks). Ours works from any curated model — demonstrated across fungus,
  plant, alga and human.
* **De-circularised evaluation.** We score against curation with the localisation-derived information
  removed (e.g. dropping Human-GEM's DeepLoc2-sourced compartments); CarveFungi reports nothing
  comparable and even hard-codes a few forced reactions.

## Where CarveFungi is genuinely stronger / different (state honestly)

* It is **end-to-end reconstruction** — it selects which reactions exist at all from a universal DB and
  compartmentalises in one pass, with **implicit gap-filling** (it adds low-score reactions when
  functionality requires). Ours assigns compartments to an existing reaction set; it does not carve or
  gap-fill.
* It **bundles a novel sequence+structure localisation predictor** (a self-contained capability),
  where we depend on an external predictor.
* **Guaranteed functionality by default** (hard biomass+ATPM), a stronger default than our optional
  ε-flux mode.

## Can we benchmark the assignment head-to-head? (why it needs the real MILP)

We tried the cheap route — emulate CarveFungi's assignment as "place each reaction in every
compartment whose DeepLoc score clears a cutoff" and compare to ours on yeast-GEM — and an
adversarial source-checked review (3 independent reviewers) **rejected it as a strawman**. The
emulation makes CarveFungi look like it over-multi-localises (≈4.5 compartments/reaction), but that is
an artefact of three mechanisms it omits, each verified against the code:

1. **Candidate set.** CarveFungi has a binary only for (reaction, compartment) pairs that *exist in its
   universal model*. `bigModelv2.21b.sbml` instantiates each enzymatic reaction in a mean of **2.16
   compartments (max 4)**, drawn almost entirely from {cytosol, mitochondrion, peroxisome, ER}. A
   reaction simply cannot be placed where no copy exists, so the emulation's 4.5 comps/rxn over the
   full 9-compartment yeast-GEM is structurally impossible in CarveFungi, and 5 of those 9
   compartments hold ~zero enzymatic copies.
2. **Connectivity.** The ε-flux indicator coupling + hard biomass/ATPM under `S·v=0`
   (`CarveMeFuncPool.py`) keeps a copy only if it carries ≥ ε flux in a connected, growing network;
   transporters are themselves scored reactions (~1e-11) that must be co-selected, so transport is not
   free. Above-cutoff copies that cannot be wired up are pruned — the emulation kept them all.
3. **Score sign.** The default reaction score is −1.0 and a below-cutoff copy gets a *negative* score
   (`(SCALE − ec_score)·(prob − cutoff)`), so `max Σ score·y` actively switches copies **off**.
   CarveFungi's objective already pushes toward fewer, higher-evidence compartments — it is not the
   parsimony-indifferent objective the emulation assumed.

The lesson: a fair head-to-head on the *assignment* requires running CarveFungi's **actual** MILP on
its real candidate set, **not** a thresholding shortcut on yeast-GEM. We did exactly that — running
CarveFungi's own `minmax_reduction` (CPLEX, unmodified) with our transport term swapped into its
objective — in [the CarveFungi head-to-head study](carvefungi_milp_benchmark.md): our transport cost
gives ~1.6× fewer transports per reaction (41% fewer) at no detectable assignment-accuracy cost. It
still carries a gold-reference caveat (mapping universal-DB reactions to curated yeast-GEM compartments
is only partial via EC/KEGG), and the big-M carve does not prove optimality, so the numbers are
deterministic, time-budget-stable near-optimal incumbents rather than certified optima.

## Bottom line for the paper

Cite CarveFungi as the contemporary network-aware compartmentalisation method, but position it as
*different in kind*: a carve-for-functionality reconstruction with a bundled fungal predictor and an
ensemble output, versus our predictor-agnostic, transport-minimising, deterministic, multi-
localisation-sound **assignment** that generalises across kingdoms. Borrow its functionality-coupling
and confidence-weighting ideas. The clean head-to-head benchmark belongs with `predictLocalization`
(same task, same lineage — see [the comparison study](predictlocalization_comparison.md)); CarveFungi
is related work, not a same-objective baseline, and a fair empirical comparison needs its real MILP
(above), not a cheap emulation.

Source read: `bin/CarveMeFuncPool.py` (MILP), `bin/EggNogScoring.py` (scoring),
`bin/compartmentPrediction/predict.py` (predictor), `bin/CreateModelEggNogPool.py` (pipeline) at
[github.com/SandraCastilloPriego/CarveFungi](https://github.com/SandraCastilloPriego/CarveFungi).
