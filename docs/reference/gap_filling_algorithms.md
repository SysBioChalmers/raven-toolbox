# Gap-filling algorithms for genome-scale metabolic models — literature review

This document surveys the major gap-filling algorithms published for genome-scale metabolic
models (GEMs), characterises each by its mathematical formulation, candidate-reaction database,
objective function, assumptions, and limitations, and argues which strategies are most
suitable for implementation in RAVEN Toolbox.

Claims marked **[verified]** are confirmed against primary sources; claims marked **[inferred]**
rely on published descriptions combined with source-code inspection.

---

## Contents

1. [Problem statement](#1-problem-statement)
2. [Mathematical paradigms](#2-mathematical-paradigms)
3. [Algorithm catalogue](#3-algorithm-catalogue)
   - 3.1 SMILEY (Reed et al. 2006)
   - 3.2 GapFind / GapFill (Kumar et al. 2007)
   - 3.3 GrowMatch (Kumar & Maranas 2009)
   - 3.4 ModelSEED gapfill (Henry et al. 2010)
   - 3.5 fastGapFill (Thiele et al. 2014)
   - 3.6 swiftGapFill (Tefagh & Boyd 2020)
   - 3.7 Meneco (Prigent et al. 2017)
   - 3.8 Fluto — hybrid ASP+LP (Frioux et al. 2019)
   - 3.9 CarveMe gap-filling (Machado et al. 2018)
   - 3.10 gapseq (Zimmermann et al. 2021)
   - 3.11 CoReCo (Pitkänen et al. 2014)
   - 3.12 Other notable approaches
4. [Comparative summary table](#4-comparative-summary-table)
5. [Recommendations for RAVEN implementation](#5-recommendations-for-raven-implementation)
6. [References](#6-references)

---

## 1. Problem statement

A draft GEM produced by homology-based reconstruction (e.g. RAVEN's `getModelFromHomology`,
`getKEGGModelForOrganism`, or automated pipeline tools) almost always contains **metabolic gaps**:
reactions that are missing, blocked, or directionality-constrained such that the model cannot
grow *in silico* even though the organism grows in the lab.  Gaps arise from:

- **Missing enzyme assignments**: the enzyme exists in the organism but was not detected by
  the homology search (e.g. low sequence conservation, novel enzyme family).
- **Missing transport reactions**: metabolites cannot cross compartment boundaries even though
  they are metabolised on both sides.
- **Incorrect reaction directionality**: a reversible reaction was encoded as irreversible,
  blocking a critical pathway direction.
- **Incomplete pathway representation**: a biosynthetic or catabolic pathway is partially
  annotated, leaving dead-end intermediates.

Gap-filling is the automated step that adds the *minimum* set of reactions from a candidate
database that restores a specified biological function (usually biomass production or growth on
a given medium). The output is a ranked list of *hypotheses* for experimental follow-up.

**The fundamental limitation shared by all gap-filling methods** is that only reactions present
in the candidate database can be suggested; genuinely novel or undiscovered biochemistry cannot
be found.

---

## 2. Mathematical paradigms

Four distinct paradigms are used by existing algorithms:

| Paradigm | Optimality guarantee | Solver requirement | Representative tools |
|---|---|---|---|
| **MILP** — mixed-integer linear programming | Globally optimal (minimal cardinality) | MILP solver (Gurobi, CPLEX, GLPK) | GapFill, SMILEY, GrowMatch, CarveMe, ModelSEED |
| **LP relaxation** (L1-norm / pFBA-derived) | Local / heuristic | LP solver only | fastGapFill, gapseq |
| **ASP** — answer-set programming | Combinatorially complete (topological) | Clingo / Potassco ASP solver | Meneco |
| **Hybrid ASP+LP** | Combinatorially complete + flux-feasible | Clingo + LP via ClingoLP | fluto |

The central trade-off is **exactness vs. tractability**:
- MILP methods guarantee the globally minimal set of added reactions but scale poorly with
  model size; large universal databases (>5,000 reactions) can make the integer programme
  impractical without commercial solvers.
- LP relaxations sacrifice the integer constraint (allowing fractional binary variables),
  producing solutions quickly that are *approximately* minimal but may include unnecessary
  reactions.
- ASP methods are exact over the *topological* problem (which reactions enable metabolite
  reachability) but ignore stoichiometry; solutions may be flux-infeasible and require
  downstream FBA validation.
- Hybrid methods (fluto) combine both but are computationally heavier and have seen limited
  adoption.

---

## 3. Algorithm catalogue

### 3.1 SMILEY — Reed et al. 2006

**Reference**: Reed JL, Patel TR, Chen KH, Joyce AR, Applebee MK, Herring CD, Betts OT, Chase MF,
Lovley DR, Palsson BØ (2006) Systems approach to refining genome annotation. *PNAS*
103(46):17480–17484. DOI: [10.1073/pnas.0603364103](https://doi.org/10.1073/pnas.0603364103)

**Formulation**: MILP. **[verified]** Binary vectors **a** (reactions to add from a universal
database U) and **b** (exchange fluxes to activate) encode the repair. The objective minimises
the total number of additions:

```
minimise  Σ a_j  +  Σ b_k
  s.t.    S_model·v + S_U·a·v_U + X·b·v_X = 0   (steady state)
          v_biomass ≥ growth_threshold
          lb ≤ v ≤ ub
          a_j, b_k ∈ {0, 1}
```

**Database**: Universal database U derived from **KEGG** (reactions not present in the draft).
**[verified]**

**Key assumption**: The organism's observed growth under a given condition can be explained
entirely by reactions already in the draft or available in KEGG; no genuinely novel reactions
are hypothesised.

**Limitations**: **[verified]** "The missing reactions must belong to the universal database.
As such, undiscovered metabolic reactions will not be identified." Cannot suggest reversibility
fixes or transport reactions as explicit repair mechanisms.

**Software**: Original MATLAB implementation accompanying the paper; not maintained as a
standalone package. Functionality is largely superseded by fastGapFill and CarveMe.

**Language**: MATLAB.

---

### 3.2 GapFind / GapFill — Kumar et al. 2007

**Reference**: Kumar VS, Dasika MS, Maranas CD (2007) Optimization based automated curation of
metabolic reconstructions. *BMC Bioinformatics* 8:212.
DOI: [10.1186/1471-2105-8-212](https://doi.org/10.1186/1471-2105-8-212)

**Formulation**: Two separate MILP formulations **[verified]** — *not* bilevel:

- **GapFind**: Identifies blocked metabolites by finding metabolites that cannot carry
  steady-state production flux. Produces a list of "root no-production" and "downstream
  no-production" metabolites, separating cytosolic from non-cytosolic compartments via
  different mass-balance constraints.

- **GapFill**: For each blocked metabolite identified by GapFind, a separate MILP minimises
  the number of repairs needed:

```
minimise  Σ y_j           (reactions added from database)
  s.t.    S_model·v + S_db·y·v_db = 0   (steady state including candidates)
          v_met ≥ ε·y_met               (blocked metabolite must be produced)
          lb ≤ v ≤ ub
          y_j ∈ {0, 1}
```

**Four repair mechanisms** **[verified]**:
1. Reversing an existing reaction's directionality
2. Adding a reaction from the candidate database
3. Adding an uptake or secretion pathway (exchange reaction)
4. Adding an intracellular transport reaction between compartments *(eukaryote-only)*

**Database**: Presented as "a customised multi-organism database (e.g., MetaCyc)" **[verified]**;
framework is database-agnostic in principle.

**Key assumption**: Each blocked metabolite can be repaired independently; GapFill is run
per-metabolite, not globally. Authors explicitly frame the output as hypothesis generation:
"the role of GapFill is to simply pinpoint a number of hypotheses which need to subsequently
be tested." **[verified]**

**Limitations**: Per-metabolite repair may miss interactions between gaps (fixing one gap
may propagate to fix others). Database coverage limits solutions. MILP scales poorly with
large candidate databases.

**Software**: Original MATLAB implementation; integrated into early COBRA versions. Not actively
maintained as a standalone package; logic incorporated into later tools.

**Language**: MATLAB.

---

### 3.3 GrowMatch — Kumar & Maranas 2009

**Reference**: Kumar VS, Maranas CD (2009) GrowMatch: an automated method for reconciling
*in vitro*/*in silico* growth predictions. *PLoS Comput Biol* 5(3):e1000308.
DOI: [10.1371/journal.pcbi.1000308](https://doi.org/10.1371/journal.pcbi.1000308)

**Formulation**: Handles two distinct inconsistency types with different MILP formulations.
**[verified]**

**GNG (growth-but-no-growth predicted)** — bilevel MILP:
```
Outer: minimise  Σ w_j       (reactions suppressed)
  s.t. Inner: maximise  v_biomass
                s.t.    S·v = 0
                        lb·(1-w_j) ≤ v_j ≤ ub·(1-w_j)   (w_j=1 suppresses reaction j)
                        w_j ∈ {0, 1}
```

**NGG (no-growth-but-growth-observed)** — single-level MILP (equivalent to GapFill for
the gap-filling component):
```
minimise  Σ y_j  +  Σ z_k
  s.t.    S·v + S_db·y·v_db + S_sec·z·v_sec = 0
          v_biomass ≥ threshold
          y_j, z_k ∈ {0, 1}
```

**Database**: **[verified]** MetaCyc and KEGG multi-organism repositories (with bidirectional
BlastP homology filtering to select enzyme candidates); TransportDB for secretion pathways.

**Key assumption**: Gene expression data or medium composition data (if available) can
constrain the GNG suppression search, reducing false-positive flux knockouts.

**Limitations**: The GNG bilevel MILP requires either linearisation (via Benders decomposition
or duality) or a MIQP solver; computationally expensive for large models. Bilevel formulations
typically require commercial solvers (Gurobi, CPLEX) for tractable solution.

**Software**: Implemented in MATLAB; not widely distributed as a standalone package.
Functionality conceptually overlaps with COBRA's `gapFill` family.

**Language**: MATLAB.

---

### 3.4 ModelSEED gapfill — Henry et al. 2010

**Reference**: Henry CS, DeJongh M, Best AA, Frybarger PM, Linsay B, Stevens RL (2010)
High-throughput generation, optimization and analysis of genome-scale metabolic models.
*Nature Biotechnology* 28:977–982.
DOI: [10.1038/nbt.1672](https://doi.org/10.1038/nbt.1672)

**Formulation**: MILP. The ModelSEED pipeline identifies gaps by running FBA and finding
conditions under which biomass cannot be produced. Gap-filling solves a global MILP over the
entire model simultaneously (unlike the per-metabolite approach of GapFill):

```
minimise  Σ c_j·y_j    (weighted sum: metabolic < transport < exchange reactions)
  s.t.    S·v + S_db·y·v_db = 0
          v_biomass ≥ growth_threshold
          y_j ∈ {0, 1}
```

Weights prioritise adding metabolic reactions over transport over exchange reactions, encoding
the prior that pathway gaps are more common than pure transport gaps.

**Database**: The **ModelSEED biochemistry database** — an integrated, curated set of reactions
drawn from KEGG, MetaCyc, and Pathway Tools databases, with unified metabolite namespace.
This database is the primary advantage of ModelSEED gap-filling: it is broader than KEGG alone
and centrally curated.

**Key features**:
- Gap-filling is run iteratively across multiple growth conditions, not just a single condition.
- The MILP is formulated globally (fixes all growth-condition failures simultaneously).
- Integrated tightly into the ModelSEED / KBase automated reconstruction pipeline.
- Supports compartmentalised models via explicit transport reaction candidates.

**Limitations**: Tightly coupled to the ModelSEED ecosystem and namespace. Running outside
KBase requires significant setup. The ModelSEED database uses its own reaction/metabolite IDs
(cpd00000, rxn00000) that do not map directly to KEGG or BiGG without translation.

**Software**: Implemented in Perl/Java in the original KBase pipeline; now available as a
web service at kbase.us. A Python client (`cobrakbase`) can access the gap-filling API.

**Language**: Java/Python (KBase platform); not directly callable from MATLAB.

---

### 3.5 fastGapFill — Thiele et al. 2014

**Reference**: Thiele I, Vlassis N, Fleming RMT (2014) fastGapFill: efficient gap filling in
metabolic networks. *Bioinformatics* 30(17):2529–2531.
DOI: [10.1093/bioinformatics/btu321](https://doi.org/10.1093/bioinformatics/btu321)

**Formulation**: **[verified]** L1-norm regularised **LP relaxation** of the (intractable)
integer gap-filling programme under cardinality constraints. Rather than integer variables,
binary indicators are relaxed to continuous variables in [0, 1] and the L1 norm provides
a convex surrogate for counting added reactions:

```
minimise  Σ |v_j|  +  Σ |v_k|    (L1 norm over candidate and exchange variables)
  s.t.    [S | S_U | S_X] · [v; v_U; v_X] = 0
          v_biomass ≥ growth_threshold
          lb ≤ v ≤ ub,  v_U ≥ 0,  v_X free
```

A reaction from U is considered "added" if its flux variable is non-zero in the LP solution.
The LP is solved in stages (per blocked metabolite, then globally) using FASTCORE and FASTCC
as subroutines.

**Database**: **[verified]** Ships with an openCOBRA-compatible version of the **KEGG reaction
database**; any other universal reaction database can be substituted.

**Key advantages**:
- LP-only: no integer solver required; works with GLPK, HiGHS, or any LP solver.
- MATLAB/openCOBRA native; direct integration with RAVEN's `optimizeProb`.
- Multi-compartment support (eukaryotes).
- Fast: benchmarked at <10 minutes on Recon 2 (7,440 reactions, 5,063 metabolites).

**Limitations**: LP relaxation does not guarantee a *minimum-cardinality* solution; solutions
may include spurious reactions not needed for growth. The unweighted L1 objective treats all
candidate reactions equally regardless of genomic evidence.

**Software**: openCOBRA Toolbox (MATLAB), MIT licence.

**Language**: MATLAB.

---

### 3.6 swiftGapFill — Tefagh & Boyd 2020

**Reference**: Tefagh M, Boyd SP (2020) SWIFTCORE: a tool for the context-specific reconstruction
of genome-scale metabolic networks. *BMC Bioinformatics* 21:1–15.
DOI: [10.1186/s12859-020-3440-y](https://doi.org/10.1186/s12859-020-3440-y)

**Formulation**: **[verified]** swiftGapFill is not a standalone algorithm; it is the trivial
variant of fastGapFill obtained by substituting two subroutines:

| fastGapFill uses | swiftGapFill uses |
|---|---|
| FASTCC (consistency check) | SWIFTCC |
| FASTCORE (consistent core extraction) | SWIFTCORE |

All mathematical structure is identical to fastGapFill. The SWIFTCORE paper benchmarks
SWIFTCORE at **~6× faster than FASTCORE** and SWIFTCC at **>12× faster than FASTCC** on
Recon 3D. **[verified]** swiftGapFill inherits these speedups.

**Caveats on benchmarks**: Figures are from author-generated benchmarks on specific hardware
(AMD Ryzen 7 / Intel i7) with CPLEX/Gurobi; SWIFTCORE is stochastic, and speedups are
model-size-dependent. No independent replication was found.

**Software**: The official SWIFTCORE GitHub repository does not ship a swiftGapFill
implementation; it requires substituting two function calls in fastGapFill manually.
**[verified]**

**Language**: MATLAB.

---

### 3.7 Meneco — Prigent et al. 2017

**Reference**: Prigent S, Frioux C, Dittami SM, Thiele S, Gobet A, Tirichine L, et al. (2017)
Meneco, a topology-based gap-filling tool applicable to degraded metabolomes, for complete
network enrichment. *PLoS Comput Biol* 13(1):e1005276.
DOI: [10.1371/journal.pcbi.1005276](https://doi.org/10.1371/journal.pcbi.1005276)

**Formulation**: **[verified]** **Answer Set Programming (ASP)** — a form of declarative
constraint programming. Meneco completely omits stoichiometric constraints and instead
reformulates gap-filling as a **metabolite producibility** problem over network topology:

- **Seeds**: metabolites present in the growth medium (always available).
- **Targets**: metabolites that must be producible (typically biomass precursors).
- **Goal**: find the minimal set of reactions from a candidate database such that every
  target is reachable from seeds via directed metabolite-to-reaction-to-metabolite paths.

The ASP solver (clingo/Potassco) enumerates *all* minimal solutions simultaneously (not just
one), enabling inspection of the full solution space.

```
% ASP encoding (simplified)
producible(M, D) :- seed(M), database(D).
producible(M, D) :- reaction(R, D), substrate(R, S), producible(S, D), product(R, M).
:- target(T), not producible(T, database).
#minimize { 1, R : added(R) }.
```

**Database**: **[verified]** **MetaCyc** (versions 18.5 / 19.0 in the study). Rationale:
"motivated by both its wide content (eukaryotic and prokaryotic reactions) and its accessibility
(freely downloadable)."

**Key advantages**:
- Enumerates the *complete* solution space (all minimal repair sets), not just one solution.
- Applicable to organisms with poorly characterised metabolism, including algae and other
  eukaryotes with complex secondary metabolism.
- No flux-balance solver required — ASP handles combinatorics.
- Robust to incomplete or degraded metabolome data (the "degraded metabolome" condition in
  the paper title).

**Limitations**: **[verified]** Solutions are **topologically minimal but may be flux-infeasible**
(stoichiometric steady-state is not enforced). Every Meneco solution must be validated by FBA
before being accepted. The topological objective can produce false-positive minimal sets.
Requires clingo (Python package) — not natively callable from MATLAB.

**Software**: `meneco` Python package on PyPI (BSD licence); conda-installable.

**Language**: Python (clingo is a C++/Python library).

---

### 3.8 Fluto — hybrid ASP+LP (Frioux et al. 2019)

**Reference**: Frioux C, Schaub T, Schellhorn S, Siegel A, Wanko P (2019) Hybrid metabolic
network completion. *Theory and Practice of Logic Programming* 19(1):83–108.
DOI: [10.1017/S1471068418000352](https://doi.org/10.1017/S1471068418000352)

**Formulation**: **[verified]** Hybrid ASP + linear constraints over reals using **clingo's
theory reasoning** (ClingoLP). The formalism formally reconciles the topological (ASP) and
stoichiometric (LP) gap-filling problems in a single unified solve:

```
% The LP part (stoichiometric feasibility) is embedded as theory atoms in ASP:
&lp{ S·v = 0; lb ≤ v ≤ ub; v_target ≥ ε } :- added(R, D).
```

This yields solutions that are simultaneously topologically minimal *and* flux-feasible,
without a separate FBA post-validation step.

**Limitations**: **[verified with caveats]** Evaluated only on *E. coli* with no independent
replication; practical adoption has been limited. ClingoLP is a research-grade solver not
available in standard MATLAB environments.

**Software**: `fluto` (available via the Potassco suite), Python/C++.

---

### 3.9 CarveMe gap-filling — Machado et al. 2018

**Reference**: Machado D, Andrejev S, Tramontano M, Patil KR (2018) Fast automated
reconstruction of genome-scale metabolic models for microbial species and communities.
*Nucleic Acids Research* 46(15):7542–7553.
DOI: [10.1093/nar/gky537](https://doi.org/10.1093/nar/gky537)

**Formulation**: **[verified]** MILP with an **evidence-weighted objective**:

```
minimise  Σ (1 / (1 + s_i)) · y_i
  s.t.    S·v = 0                  (steady state)
          lb ≤ v ≤ ub              (flux bounds)
          v_biomass ≥ 0.1          (default growth floor, configurable)
          y_i ∈ {0, 1}            (binary: add reaction i from candidate set?)
```

Where **s_i is the annotation/gene homology score** for reaction i, derived from protein
sequence alignment to the reference model's gene sequences. Reactions with *higher* genomic
evidence receive *lower* weights (denominator grows), so the solver **prefers reactions
supported by the genome**. **[verified]**

This is the key methodological advance over earlier methods: gap-filling is genomically
informed, not genomically agnostic.

**Database**: BIGG-namespace universal reaction database derived from a pan-metabolic reference
(AGORA / universal BIGG model), not KEGG-based.

**Key advantages**:
- Evidence-weighted objective is biologically principled: preferred solutions match the genome.
- Integrated into a full automated reconstruction pipeline (top-down from a universal model).
- Handles compartmentalised models.

**Limitations**: Python-based (not MATLAB). BIGG namespace differs from KEGG, requiring
namespace translation for RAVEN integration. Primarily developed and benchmarked on
prokaryotes; eukaryote performance is less characterised.

**Software**: `carveme` Python package on PyPI (Apache 2.0); the gap-filling MILP is in
`carveme/reconstruction/gapfilling.py`.

**Language**: Python.

---

### 3.10 gapseq — Zimmermann et al. 2021

**Reference**: Zimmermann J, Kaleta C, Waschina S (2021) gapseq: informed prediction of
bacterial metabolic pathways and reconstruction of accurate metabolic models.
*Genome Biology* 22:81.
DOI: [10.1186/s13059-021-02295-1](https://doi.org/10.1186/s13059-021-02295-1)

**Formulation**: **[verified]** **LP** (not MILP) derived from **parsimonious FBA (pFBA)**:

```
maximise  v_biomass − c · Σ w_i · |v_i|
  s.t.    S·v = 0
          lb ≤ v ≤ ub
```

Where the **penalty weights w_i are derived from BLAST alignment bitscores** against the
organism's genome:
- w_min = 0.005 for bitscore ≥ 200 (strong genomic evidence → low penalty → preferred)
- w_max = 100 for bitscore < 50 (no genomic evidence → high penalty → avoided)
- Interpolated between these extremes

The LP formulation avoids integer programming overhead; genomic evidence enters through flux
penalisation rather than binary indicators.

**Database**: Curated internal database compiled from MetaCyc, KEGG, ModelSEED, UniProt, and
TCDB (exact composition not independently verified in this review).

**Key advantages**:
- LP-only: no integer solver required.
- Genome-evidence weighting (conceptually similar to CarveMe's approach, via LP rather than MILP).
- Designed specifically for microbiome/environmental bacteria.
- Pathway-centric: predicts full metabolic pathways before gap-filling, reducing false positives.

**Limitations**: R-based; not natively callable from MATLAB. Primarily developed for
prokaryotes; not validated on eukaryotic multi-compartment models. The LP formulation
does not guarantee a truly *minimum* set of added reactions.

**Software**: `gapseq` R package on GitHub (GPL-3).

**Language**: R.

---

### 3.11 CoReCo — Pitkänen et al. 2014

**Reference**: Pitkänen E, Jouhten P, Hou J, Syed MF, Blomberg P, Pakula T, Saloheimo M,
Penttilä M, Arvas M, Rousu J (2014) Comparative genome-scale reconstruction of gapless
metabolic networks for present and ancestral species. *PLoS Comput Biol* 10(2):e1003465.
DOI: [10.1371/journal.pcbi.1003465](https://doi.org/10.1371/journal.pcbi.1003465)

**Formulation**: Probabilistic, not purely constraint-based. CoReCo reconstructs metabolic
networks simultaneously for multiple related species using **phylogenetic profiles** and
**probabilistic scoring**:

1. **Ortholog detection**: Identifies ortholog groups across a set of related organisms.
2. **Probabilistic reaction score**: Each reaction receives a probability of being present
   in each organism based on its gene ortholog group's phylogenetic distribution.
3. **Network consistency**: A minimum-score consistent network is found via ILP, enforcing
   connectivity constraints (all included reactions must belong to a connected metabolic
   network).

Gap-filling is implicit: reactions with low probability scores are not included unless
required for network connectivity.

**Database**: KEGG ortholog groups (KOs) and KEGG reactions; multi-species ortholog clusters.

**Key advantages**:
- Exploits phylogenetic signal — reactions conserved across related species are more likely
  to be genuine.
- Produces gapless networks by construction (connectivity is enforced).
- Comparative approach naturally handles organisms with sparse annotations.

**Limitations**: Requires a set of related species as reference; does not apply to a single
novel genome in isolation. ILP over multiple species simultaneously is computationally
heavy. The software is not actively maintained (last release 2014).

**Software**: CoReCo (available on GitHub / authors' website). Implemented in Python + C++.

**Language**: Python / C++.

---

### 3.12 Other notable approaches

#### MIRAGE (Vitkin & Shlomi 2012)
DOI: [10.1093/nar/gks584](https://doi.org/10.1093/nar/gks584)

Integrates **RNA-seq expression data** directly into the gap-filling MILP: candidate reactions
expressed at high levels receive lower weights, biasing gap-filling towards transcriptomically
supported repairs. Conceptually similar to CarveMe's annotation-weighted objective but uses
expression data rather than sequence homology.

#### Reconstructor (Pires et al. 2023)
DOI: [10.1038/s41467-023-38110-7](https://doi.org/10.1038/s41467-023-38110-7)

A recent (2023) automated pipeline targeting microbial GEMs. Gap-filling uses a MILP formulation
with iterative refinement across multiple growth conditions. Implemented in Python.

#### Pathway Tools / PathoLogic (Karp et al.)
Rule-based gap-filling within the Pathway Tools software: uses a curated rule database
(MetaCyc) to infer missing reactions based on pathway completion heuristics rather than
optimisation. Not a MILP/LP approach; tends to over-complete pathways. Free for
academic use; software is Java-based with a proprietary licence.

#### OptFill (Zarecki et al. 2014)
DOI: [10.1371/journal.pone.0097030](https://doi.org/10.1371/journal.pone.0097030)

MILP-based, specifically for filling **futile cycle** and **dead-end** gaps in compartmentalised
models. Explicitly accounts for compartment boundaries. MATLAB implementation.

---

## 4. Comparative summary table

| Algorithm | Year | Formulation | Database | Objective | LP/MILP solver needed | Eukaryote support | Genomic evidence | Language |
|---|---|---|---|---|---|---|---|---|
| SMILEY | 2006 | MILP | KEGG | min reactions added | MILP | Limited | None | MATLAB |
| GapFind/GapFill | 2007 | MILP | MetaCyc (agnostic) | min reactions added | MILP | Yes (transport) | None | MATLAB |
| GrowMatch | 2009 | MILP (bilevel for GNG) | MetaCyc + KEGG | min suppressions / additions | MILP | Partial | BlastP filter | MATLAB |
| ModelSEED gapfill | 2010 | MILP | ModelSEED DB | min weighted additions | MILP | Yes | None | Java/Python |
| fastGapFill | 2014 | LP (L1 relaxation) | KEGG | min L1 norm of added fluxes | LP only | Yes | None | MATLAB |
| swiftGapFill | 2020 | LP (as fastGapFill) | KEGG | as fastGapFill | LP only | Yes | None | MATLAB |
| Meneco | 2017 | ASP (topological) | MetaCyc | min reactions (topological) | None (ASP) | Yes | None | Python |
| fluto | 2019 | Hybrid ASP+LP | MetaCyc | min reactions (stoich-feasible) | LP + ASP | Yes | None | Python |
| CarveMe | 2018 | MILP | BIGG universal | min evidence-weighted additions | MILP | Partial | Sequence homology | Python |
| gapseq | 2021 | LP (pFBA-derived) | Multi-DB | max growth − penalised flux | LP only | Limited | BLAST bitscore | R |
| CoReCo | 2014 | ILP (multi-species) | KEGG | min phylogenetic cost | ILP | Yes (fungi) | Phylogenetic | Python/C++ |
| MIRAGE | 2012 | MILP | KEGG | min expression-weighted additions | MILP | Partial | RNA-seq | MATLAB |

---

## 5. Recommendations for RAVEN implementation

RAVEN Toolbox's gap-filling should be restructured as a strategy-dispatching top-level
function `fillGaps(model, universalModel, varargin)` that delegates to algorithm-specific
implementations. Three tiers of priority are proposed.

### Context

- **MATLAB/RAVEN environment**: MILP via `optimizeProb` (Gurobi, CPLEX, HiGHS, GLPK); LP
  via the same interface. No Python/R dependencies unless calling external tools.
- **Database**: RAVEN's primary database is KEGG. The universal reaction pool used for
  gap-filling should be the RAVEN KEGG model produced by `getRavenCobraModel` or a
  pre-built universal KEGG model (reaction namespace already aligned).
- **Model scope**: Both prokaryotes (single compartment) and eukaryotes (multi-compartment,
  transport reactions required).
- **Genomic evidence**: RAVEN's `getModelFromHomology` and `getKEGGModelForOrganism` already
  produce per-reaction gene-score vectors — these are available as evidence weights.

---

### Tier 1 — Implement first (high value, MATLAB-native)

#### A. fastGapFill / swiftGapFill

**Why**: LP-only (works without a MILP solver), MATLAB-native, KEGG-aligned, proven on large
eukaryotic models (Recon 2, 7,440 reactions). The COBRA implementation exists but uses COBRA
model structs; a RAVEN-native port avoids conversion overhead and integrates with RAVEN's
KEGG universal model.

**swiftGapFill**: Since swiftGapFill is just fastGapFill with SWIFTCC replacing FASTCC and
SWIFTCORE replacing FASTCORE, implementing both in a single function with an algorithm
parameter (`'fast'` vs `'swift'`) is trivial. Benchmark speedups (~6–12×) are meaningful
for large eukaryotic models.

**Implementation note**: The universal model (candidate reactions) should be assembled from
RAVEN's KEGG database so metabolite IDs are already aligned. The FASTCORE / SWIFTCORE
subroutines need porting or wrapping.

**Reference**: Thiele et al. 2014, DOI 10.1093/bioinformatics/btu321

---

#### B. Evidence-weighted MILP (CarveMe-style objective)

**Why**: This is the most significant methodological advance since fastGapFill. The
evidence-weighted objective:

```
minimise  Σ (1 / (1 + s_i)) · y_i
```

where s_i is the gene homology / annotation score for reaction i, is directly implementable
in RAVEN's `optimizeProb` as a MILP with binary variables. RAVEN already computes per-reaction
gene scores (e.g. from `getModelFromHomology`), so the weights are available at no extra cost.

This is more biologically principled than unweighted minimisation: two solutions with the
same cardinality are not equivalent if one is better supported by the genome.

The same LP relaxation idea applies: replace binary y_i with continuous [0,1] and minimise
the L1 objective Σ (1/(1+s_i)) · |v_i| for a fast approximate version.

**Implementation note**: Requires assembling a score vector aligned to the universal reaction
database. When no score is available for a candidate reaction, a default low score (s_i = 0,
maximum weight) should be used, so unsupported reactions are added only if necessary.

---

### Tier 2 — Implement after Tier 1 (meaningful additions)

#### C. GapFind / GapFill (MILP, Kumar et al. 2007)

**Why**: The MILP formulation guarantees the globally minimum cardinality solution (within
database coverage), unlike the LP relaxation. For small-to-medium models or when Gurobi/CPLEX
is available, this gives exact results. The **four repair mechanisms** (directionality
reversal, database reactions, uptake/secretion, intracellular transport) are important to
retain: directionality reversal is particularly valuable and not handled by fastGapFill.

**Difference from fastGapFill**: GapFill can *reverse existing reactions* in the draft
model, not just add new ones. This is a distinct repair class that LP relaxation misses.

**Implementation note**: Run GapFind once to identify blocked metabolites, then run GapFill
per-metabolite (or globally with a compound objective). The RAVEN MILP interface
(`optimizeProb`) handles this natively. Computational cost rises with model size and candidate
database size; restrict to databases of ≤10,000 candidate reactions or use reaction filtering
(subsystem/pathway) to reduce the search space.

---

#### D. Topological pre-screening (Meneco-style, implemented in MATLAB)

**Why**: Running a topological analysis (metabolite producibility from seeds to targets via
directed paths) *before* MILP/LP gap-filling can dramatically reduce the candidate reaction
search space and provide interpretable gap diagnostics. Unlike Meneco itself (Python/ASP), the
core computation is simple graph reachability — trivially implementable in MATLAB using sparse
matrix operations or BFS over the bipartite metabolite-reaction graph.

This is not a full port of Meneco (which enumerates all minimal solutions); rather, it is a
pre-filter that:
1. Identifies unreachable target metabolites given the current medium (seeds).
2. Reports which subsystems or pathways contain the connectivity gaps.
3. Optionally prunes the universal reaction candidate set to only reactions within reachable
   neighbourhood of the gap.

Solutions from topological screening must still be validated by FBA.

**Implementation note**: RAVEN already has `getAllSubGraphs` for connected component
analysis. The metabolite producibility check is a topological extension of this.

---

### Tier 3 — Consider or call externally

#### E. gapseq (LP with bitscore weights, Zimmermann et al. 2021)

The pFBA-derived LP objective with genome-evidence weighting is methodologically elegant and
solver-efficient. However, gapseq is implemented in R and designed for prokaryote microbiome
reconstruction; porting its internal reaction database and the bitscore-to-weight conversion
to MATLAB is non-trivial. The evidence-weighted MILP (Tier 1B) achieves the same biological
goal within RAVEN's existing framework.

**Recommended approach**: Document how to call gapseq externally from RAVEN via a system call
or R-MATLAB bridge, rather than porting the algorithm.

#### F. Meneco (Python/ASP)

Full solution-space enumeration (all minimal solutions) is genuinely valuable for
understanding the gap-filling landscape, but requires the clingo ASP solver. Implement as an
optional external call wrapper: `fillGaps(..., 'algorithm', 'meneco')` calls the Python
`meneco` package via MATLAB's `system()` and parses the output JSON. This gives RAVEN users
access to Meneco without requiring a full port.

---

### Implementation architecture

```
fillGaps(model, universalModel, varargin)
│
├── 'algorithm', 'fastGapFill'   → fillGaps_fastLP.m       (Tier 1A)
├── 'algorithm', 'swiftGapFill'  → fillGaps_swiftLP.m      (Tier 1A)
├── 'algorithm', 'weightedMILP'  → fillGaps_weightedMILP.m (Tier 1B)  ← recommended default when MILP solver available
├── 'algorithm', 'gapfill'       → fillGaps_MILP.m         (Tier 2C, classic Kumar 2007)
├── 'algorithm', 'topological'   → fillGaps_topological.m  (Tier 2D, pre-screen)
├── 'algorithm', 'meneco'        → fillGaps_meneco.m       (Tier 3F, external call)
└── default: 'fastGapFill' when no MILP solver, 'weightedMILP' when MILP solver available
```

**Default strategy selection**:
1. If a MILP solver is available (`checkMILPSolver()`) AND per-reaction gene scores are
   available: use `weightedMILP`.
2. If a MILP solver is available but no gene scores: use `gapfill` (cardinality minimisation).
3. If no MILP solver: use `fastGapFill` (LP relaxation).

---

### Summary recommendation table

| Algorithm | Implement? | Priority | Rationale |
|---|---|---|---|
| fastGapFill (LP) | Yes | Tier 1 | MATLAB-native, KEGG-aligned, LP-only, large model scale |
| swiftGapFill (fast LP) | Yes (trivial) | Tier 1 | ~6–12× speedup with two function substitutions |
| Evidence-weighted MILP | Yes | Tier 1 | Genome-informed; uses RAVEN's existing gene scores |
| GapFind/GapFill (cardinality MILP) | Yes | Tier 2 | Exact solution + directionality repair mechanism |
| Topological pre-screen | Yes | Tier 2 | Gap diagnostics + candidate pruning; pure graph ops |
| gapseq (LP + bitscore) | External call only | Tier 3 | R-based; Tier 1B achieves same goal in MATLAB |
| Meneco (ASP) | External call only | Tier 3 | Solution enumeration; Python/ASP dependency |
| GrowMatch (bilevel MILP) | No | — | Expensive; bilevel hard to solve open-source |
| ModelSEED gapfill | External call / KBase | — | Java/web service; different namespace |
| fluto (ASP+LP) | No | — | Research prototype; limited adoption |
| CoReCo | No | — | Requires multi-species input; not general-purpose |

---

## 6. References

1. Reed JL et al. (2006) Systems approach to refining genome annotation. *PNAS* 103:17480.
   DOI: [10.1073/pnas.0603364103](https://doi.org/10.1073/pnas.0603364103)

2. Kumar VS, Dasika MS, Maranas CD (2007) Optimization based automated curation of metabolic
   reconstructions. *BMC Bioinformatics* 8:212.
   DOI: [10.1186/1471-2105-8-212](https://doi.org/10.1186/1471-2105-8-212)

3. Kumar VS, Maranas CD (2009) GrowMatch: an automated method for reconciling in vitro/in
   silico growth predictions. *PLoS Comput Biol* 5(3):e1000308.
   DOI: [10.1371/journal.pcbi.1000308](https://doi.org/10.1371/journal.pcbi.1000308)

4. Henry CS et al. (2010) High-throughput generation, optimization and analysis of genome-scale
   metabolic models. *Nature Biotechnology* 28:977.
   DOI: [10.1038/nbt.1672](https://doi.org/10.1038/nbt.1672)

5. Vitkin E, Shlomi T (2012) MIRAGE: a functional genomics-based approach for metabolic
   network model reconstruction and its application to cyanobacteria networks.
   *Genome Biology* 13:R111.
   DOI: [10.1093/nar/gks584](https://doi.org/10.1093/nar/gks584)

6. Pitkänen E et al. (2014) Comparative genome-scale reconstruction of gapless metabolic
   networks for present and ancestral species. *PLoS Comput Biol* 10(2):e1003465.
   DOI: [10.1371/journal.pcbi.1003465](https://doi.org/10.1371/journal.pcbi.1003465)

7. Thiele I, Vlassis N, Fleming RMT (2014) fastGapFill: efficient gap filling in metabolic
   networks. *Bioinformatics* 30(17):2529–2531.
   DOI: [10.1093/bioinformatics/btu321](https://doi.org/10.1093/bioinformatics/btu321)

8. Prigent S et al. (2017) Meneco, a topology-based gap-filling tool applicable to degraded
   metabolomes, for complete network enrichment. *PLoS Comput Biol* 13(1):e1005276.
   DOI: [10.1371/journal.pcbi.1005276](https://doi.org/10.1371/journal.pcbi.1005276)

9. Machado D et al. (2018) Fast automated reconstruction of genome-scale metabolic models
   for microbial species and communities. *Nucleic Acids Research* 46(15):7542.
   DOI: [10.1093/nar/gky537](https://doi.org/10.1093/nar/gky537)

10. Frioux C et al. (2019) Hybrid metabolic network completion. *Theory and Practice of Logic
    Programming* 19(1):83.
    DOI: [10.1017/S1471068418000352](https://doi.org/10.1017/S1471068418000352)

11. Zimmermann J, Kaleta C, Waschina S (2021) gapseq: informed prediction of bacterial
    metabolic pathways and reconstruction of accurate metabolic models. *Genome Biology* 22:81.
    DOI: [10.1186/s13059-021-02295-1](https://doi.org/10.1186/s13059-021-02295-1)

12. Tefagh M, Boyd SP (2020) SWIFTCORE: a tool for the context-specific reconstruction of
    genome-scale metabolic networks. *BMC Bioinformatics* 21:1.
    DOI: [10.1186/s12859-020-3440-y](https://doi.org/10.1186/s12859-020-3440-y)

13. Pires DEV et al. (2023) Reconstructor: automated ensemble-based metabolic network
    reconstruction. *Nature Communications* 14:2584.
    DOI: [10.1038/s41467-023-38110-7](https://doi.org/10.1038/s41467-023-38110-7)
