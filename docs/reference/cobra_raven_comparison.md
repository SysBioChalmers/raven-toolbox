# COBRA Toolbox vs RAVEN Toolbox — Functionality Comparison

This document compares the COBRA Toolbox (MATLAB) and RAVEN Toolbox (MATLAB) feature-by-feature.
The primary aim is to identify what COBRA has that RAVEN lacks, with particular attention to
genome-scale model reconstruction — RAVEN's core purpose — but covering all functional areas.

**Sources**

- COBRA Toolbox: v3.x documentation, module index at opencobra.github.io, and the COBRA v3.0 paper (arXiv:1710.04038)
- RAVEN Toolbox: RAVEN v2.x source tree (`analysis/`, `reconstruction/`, `io/`, `external/`, `core/`)

---

## Contents

1. [Executive summary](#1-executive-summary)
2. [Genome-scale model reconstruction](#2-genome-scale-model-reconstruction)
3. [Model editing and manipulation](#3-model-editing-and-manipulation)
4. [Flux balance analysis and LP/MILP solving](#4-flux-balance-analysis-and-lpmilp-solving)
5. [Flux variability and enumerating alternatives](#5-flux-variability-and-enumerating-alternatives)
6. [Knockout and deletion analysis](#6-knockout-and-deletion-analysis)
7. [Strain design and metabolic engineering](#7-strain-design-and-metabolic-engineering)
8. [Omics-driven context-specific modelling](#8-omics-driven-context-specific-modelling)
9. [Thermodynamics](#9-thermodynamics)
10. [Flux sampling](#10-flux-sampling)
11. [Multi-species and microbiome modelling](#11-multi-species-and-microbiome-modelling)
12. [Network topology and elementary flux modes](#12-network-topology-and-elementary-flux-modes)
13. [Visualization and reporting](#13-visualization-and-reporting)
14. [Metabolomics integration](#14-metabolomics-integration)
15. [Model I/O and database connectivity](#15-model-io-and-database-connectivity)
16. [Summary gap table](#16-summary-gap-table)
17. [RAVEN strengths not in COBRA](#17-raven-strengths-not-in-cobra)

---

## 1. Executive summary

RAVEN and COBRA serve overlapping but distinct purposes.

**RAVEN's core strength** is automated genome-scale model reconstruction: pulling protein sequences
from a genome, running HMM-based homology search against curated databases, scoring gene-reaction
associations, and extracting tissue-specific models via tINIT/ftINIT from transcriptomic data.
Its analysis toolbox covers the essential FBA operations but is narrower than COBRA's.

**COBRA's core strength** is constraint-based analysis breadth: thermodynamic integration,
advanced strain design (OptKnock, OptForce, gMCS), multiple loopless and parsimonious FBA variants,
Markov-chain flux sampling, microbiome community models, and whole-body physiology models.
Its reconstruction pipeline (DEMETER, rBioNet) targets microbial and human databases and
is less focused on de-novo genome-based reconstruction.

For a lab doing **human tissue or yeast GEM reconstruction** from genomic or transcriptomic data,
RAVEN covers the critical path well. For labs doing **thermodynamically constrained analysis**,
**strain design**, **microbiome modelling**, or **multi-omics integration** beyond transcriptomics,
COBRA has capabilities that RAVEN lacks entirely.

---

## 2. Genome-scale model reconstruction

### RAVEN

| Function | Description |
|---|---|
| `getKEGGModelForOrganism` | Full KEGG-based draft reconstruction from organism taxonomy ID; pulls protein sequences, runs DIAMOND BLASTP against KEGG protein database, scores gene-reaction associations |
| `getModelFromHomology` | Template-based reconstruction via bidirectional HMM/BLAST homology; the backbone of high-quality GEM reconstruction from reference models |
| `makeRaven` | Assembles a RAVEN model struct from scored homology results |
| `addRxnsGenesMets` | Adds reactions, genes, and metabolites to an existing draft |
| `gapReport` | Reports metabolic gaps (dead ends, disconnected subsystems) |
| `fillGaps` | LP-based heuristic gap-filling against a universal reaction database |
| `runINIT` | INIT (Integrative Network Inference for Tissues) — tissue-specific draft from expression scores |
| `tINIT` | The original INIT implementation |
| `ftINIT` | Fast tINIT with MILP-based extraction; state-of-the-art for human tissue models |
| `getINITModel2` | Wrapper calling tINIT/ftINIT with HPA/GTEx data loading |
| `checkModel` | Validates RAVEN model struct fields |
| `simplifyModel` | Removes empty genes/metabolites, normalises identifiers |

### COBRA

| Function | Description |
|---|---|
| `fastGapFill` | MILP-based gap-filler (Thiele et al. 2014); fills stoichiometric gaps using a universal reaction database (KEGG + VMH) including transport reactions |
| `swiftGapFill` | Alternative MILP gap-filler using SWIFTCORE |
| `prepareFastGapFill` | Builds the S/U/X matrix input for fastGapFill |
| `generateSUXMatrix` | Assembles combined stoichiometric + universal + transport matrices |
| `gapFind` | Identifies blocked metabolites (anything downstream of a gap) |
| `identifyBlockedRxns` | Detects all blocked reactions via FASTCORE |
| `postProcessGapFillSolutions` | Classifies added reactions (metabolic/transport/exchange) |
| `fastcc` | Tests stoichiometric consistency; returns the flux-consistent subset |
| `checkCobraModelUnique` | Checks uniqueness of reaction/metabolite names |
| `checkDuplicateRxn` | Detects and removes duplicate reactions |
| `checkModelPreFBA` | Validates stoichiometric and flux consistency before FBA |
| `detectDeadEnds` | Returns metabolites participating in only one reaction |
| `fastLeakTest` | Tests whether metabolites leak in a closed model |
| `removeDeadEnds` | Removes dead-end metabolites and their reactions |
| DEMETER pipeline | Automated refinement of microbial draft reconstructions from KBase/ModelSEED; includes taxonomy-based GPR refinement and 100+ quality tests |
| rBioNet | GUI-based reconstruction tool with VMH database integration |
| Model borgifier | Submodule for merging and comparing two reconstructions |

### Gaps (COBRA has, RAVEN lacks)

- **MILP gap-filling**: RAVEN's `fillGaps` is an LP heuristic; COBRA's `fastGapFill` / `swiftGapFill`
  use MILP to find the globally minimal set of reactions. This produces more parsimonious fills.
- **Stoichiometric consistency checking**: `fastcc`, `fastLeakTest`, `checkDuplicateRxn`,
  `detectDeadEnds`. RAVEN's `checkModel` validates struct fields but does not test flux consistency.
- **DEMETER / KBase pipeline**: For building or refining microbial community member models
  from KBase/ModelSEED outputs, COBRA has a dedicated automated pipeline.
- **VMH / BiGG-centred reconstruction**: If building models against the VMH or BIGG namespace,
  COBRA's rBioNet and DEMETER pipelines are tailored for it.

### RAVEN strengths not in COBRA

- **KEGG-centric genome reconstruction**: `getKEGGModelForOrganism` is a mature, maintained
  pipeline that maps a genome directly to KEGG reactions and produces a draft requiring only
  gap-filling. COBRA has no equivalent end-to-end KEGG-based reconstruction.
- **HMM-based homology scoring**: `getModelFromHomology` supports HMMER as the search engine,
  giving finer granularity on gene confidence than BLAST alone.
- **ftINIT quality**: For human tissue-specific models based on HPA/GTEx protein/RNA data,
  ftINIT (RAVEN) is the state-of-the-art; COBRA's equivalent methods (GIMME, iMAT, FASTCORE)
  predate it and are less accurate on human models.

---

## 3. Model editing and manipulation

### RAVEN

RAVEN has a comprehensive set of model manipulation utilities:

`addRxns`, `addMets`, `addGenes`, `removeRxns`, `removeMets`, `removeGenes`,
`setParam` (constraint modification), `setExchangeBounds`, `getExchangeRxns`,
`contractModel` (merge identical reactions), `expandModel` (split reversible reactions),
`convertToReversible`, `convertToIrreversible`, `mergeModels`, `combineModels`,
`getModelDiff`, `sortModel`, `copyToRaven` (struct copy), `exportForGit`
(deterministic field ordering for version control).

### COBRA

`addReaction`, `removeRxns`, `addMetabolite`, `changeRxnBounds`, `changeObjective`,
`convertToIrreversible`, `mergeModels`, `findExcRxns`, `extractSubNetwork`,
`subNetwork`, `buildRxnGeneMat`, `tncore` (gene essentiality-driven network reduction),
`pruneModel` (removes blocked reactions below a flux threshold).

### Assessment

Both toolboxes cover the essentials. RAVEN's `exportForGit` and `getModelDiff`
are distinctive for version-controlled model development. COBRA's `tncore` and
`extractSubNetwork` are useful for focused subnetwork analysis.

---

## 4. Flux balance analysis and LP/MILP solving

### RAVEN

| Function | Description |
|---|---|
| `solveLP` | Standard LP: maximise/minimise `c'v` subject to `Sv=0`, `lb≤v≤ub` |
| `optimizeProb` | General LP/MILP solver wrapping GLPK, Gurobi, CPLEX, MOSEK |
| `getMinNrFluxes` | MILP: find flux distribution with minimum number of active reactions |
| `getAllowedBounds` | Flux variability: min/max flux for each reaction |
| `haveFlux` | Identifies reactions that can carry non-zero flux |
| `getEssentialRxns` | Returns reactions whose deletion abolishes objective |
| `getMinimalMedium` | MILP: minimal set of uptake exchanges that allow growth |
| `runDynamicFBA` | Dynamic FBA (Euler method, static optimisation) |

### COBRA

| Function | Description |
|---|---|
| `optimizeCbModel` | Flexible FBA: LP, MILP, QP, zero-norm cardinality |
| `pFBA` | Parsimonious FBA (Lewis et al. 2010): optimise objective, then minimise total flux; classifies gene usage |
| `geometricFBA` | Finds unique central point inside optimal FBA polytope (Smallbone & Simeonidis 2009) |
| `sparseFBA` | Finds minimal set of active reactions achieving objective (L0) |
| `relaxedFBA` | Finds minimal relaxations of bounds to make infeasible model feasible |
| `relaxFBA_cappedL1` | Same via capped-L1 norm approximation |
| `enumerateOptimalSolutions` | Finds multiple alternative optimal flux distributions |
| `minimizeModelFlux` | Minimises total flux (taxicab, zero-norm, or Euclidean) |
| `cycleFreeFlux` | Removes stoichiometrically balanced cycles from FBA solutions |
| `addLoopLawConstraints` | Adds thermodynamic loop-law constraints to LP/MILP |
| `dynamicFBA` | Dynamic FBA with substrate/biomass dynamics (Varma & Palsson 1994) |
| `dynamicRFBA` | Dynamic regulatory FBA; time evolution with Boolean regulatory network |
| `mdFBA` | Metabolite dilution FBA — accounts for growth-rate-dependent metabolite dilution |

### Gaps

- **pFBA**: RAVEN has no parsimonious FBA. Gene classification (essential / pFBA-optimal /
  enzymatically complex / metabolically blocked) requires COBRA.
- **Loopless FBA**: RAVEN has no cycle removal or loop-law constraints. Solutions can contain
  thermodynamically infeasible cycles; COBRA's `cycleFreeFlux` / `addLoopLawConstraints` remove them.
- **Geometric FBA**: No equivalent in RAVEN.
- **Infeasibility relaxation**: RAVEN returns infeasibility; `relaxedFBA` finds the *minimal*
  bound adjustments needed to restore feasibility.
- **Multiple optimal solutions**: RAVEN has no equivalent to `enumerateOptimalSolutions`.
- **Regulatory dynamic FBA**: `dynamicRFBA` integrates a Boolean gene-regulatory network
  over time; RAVEN has no equivalent.

---

## 5. Flux variability and enumerating alternatives

### RAVEN

`getAllowedBounds` — basic FVA for a subset of reactions. Single-threaded.

### COBRA

| Function | Description |
|---|---|
| `fluxVariability` | Full FVA; supports loopless mode, multi-threading, and per-reaction objective constraints |
| `fastFVA` | Accelerated FVA via multi-threaded CPLEX solver |
| `mtFVA` | Multi-threaded FVA |
| `fvaJaccardIndex` | Compares flux ranges between two models using Jaccard index |

### Gaps

- Loopless FVA mode — eliminates cycle-containing solutions from variability ranges.
- Multi-threaded FVA — faster for large models.
- Model-to-model FVA comparison (`fvaJaccardIndex`).

---

## 6. Knockout and deletion analysis

### RAVEN

| Function | Description |
|---|---|
| `findGeneDeletions` | Single gene deletion via FBA, returns growth for each KO |

### COBRA

| Function | Description |
|---|---|
| `singleGeneDeletion` | Single gene KO via FBA, MOMA, or linearMOMA |
| `singleRxnDeletion` | Single reaction KO |
| `doubleGeneDeletion` | All pairwise gene knockout combinations |
| `MOMA` | Quadratic MOMA — minimises flux adjustment between wild-type and KO |
| `linearMOMA` | Linear MOMA variant (LP) |
| `ROOM` | Regulatory on/off minimisation (MILP) — minimises number of altered reactions |
| `linearROOM` | LP relaxation of ROOM |
| `findEpistaticInteractions` | Double deletion analysis for synthetic sick/lethal pairs |
| `essentialRxn4MultipleModels` | Reaction deletion across multiple model variants |
| Fast-SL submodule | Identifies all synthetic lethal gene pairs |

### Gaps

- **Quadratic / ROOM knockout prediction**: After a gene deletion, MOMA predicts the
  most likely adaptive flux redistribution by minimising deviation from wild-type fluxes;
  ROOM minimises the number of reactions that change. Both are more biologically realistic
  than simply re-optimising the objective. RAVEN has neither.
- **Double gene deletion**: RAVEN only supports single-gene knockouts. COBRA's
  `doubleGeneDeletion` runs pairwise combinatorial analysis.
- **Synthetic lethality**: RAVEN has no equivalent to Fast-SL.

---

## 7. Strain design and metabolic engineering

### RAVEN

| Function | Description |
|---|---|
| `FSEOF` | Flux Scanning based on Enforced Objective Flux — ranks reactions by co-flux with target |
| `runProductionEnvelope` | Production envelope for biomass vs. target product |
| `runRobustnessAnalysis` | Objective sensitivity as a reaction's flux varies |
| `runPhenotypePhasePlane` | 2-reaction phenotype phase plane |
| `runSimpleOptKnock` | Simplified single-level OptKnock approximation |

### COBRA

**OptKnock family**

| Function | Description |
|---|---|
| `OptKnock` | Bilevel MILP to find knockout sets that couple biomass growth with target overproduction |
| `analyzeOptKnock` | Validates and analyses OptKnock solutions |
| `createBilevelMILPproblem` | Low-level bilevel MILP builder |
| `optGene` | Evolutionary programming-based knockout strategy |
| `GDLS` | Genetic Design by Local Search — fast neighbourhood search |

**OptForce**

| Function | Description |
|---|---|
| `FVAOptForce` | FVA for wild-type vs. mutant with production constraints |
| `findMustU` / `findMustL` | First-order must-upregulate / must-downregulate sets |
| `findMustUU` / `findMustUL` / `findMustLL` | Second-order must-sets (pairwise) |
| `optForce` | Finds minimal k-intervention sets for targeted overproduction |
| `analyzeOptForceSol` | Evaluates proposed interventions |

**Minimal Cut Sets (drug targets / growth coupling)**

| Function | Description |
|---|---|
| `calculateGeneMCS` | Computes genetic Minimal Cut Sets (gMCSs) using CPLEX warm-start |
| `calculateMCS` | Reaction-level minimal cut sets |
| `checkGeneMCS` | Validates minimality of computed gMCSs |

**OptEnvelope**

| Function | Description |
|---|---|
| `optEnvelope` | Finds minimum knockout set giving equivalent production envelope |
| `minActiveRxns` / `sequentialOEReinserts` / `milpOEReinserts` | Supporting envelope optimisation |

**Other**

| Function | Description |
|---|---|
| `productionEnvelope` | Richer production envelope with flux scanning variants |
| `theoretMaxProd` | Theoretical maximum production for a target |
| `phenotypePhasePlane` | 3D surface + shadow prices for two reactions |

### Gaps

COBRA's strain design toolbox is far broader than RAVEN's:

- **True OptKnock** (bilevel MILP, not heuristic): RAVEN's `runSimpleOptKnock` uses a
  single-level approximation; COBRA has the full bilevel formulation.
- **OptForce**: Identifies reactions that *must* change (up or down) to achieve a production
  target — a complementary perspective to knockout-only methods.
- **GDLS / optGene**: Heuristic design strategies for when the MILP becomes intractable.
- **Genetic MCS**: `calculateGeneMCS` computes Minimal Cut Sets at the gene level — essential
  for identifying combination drug targets or growth-coupled designs.

---

## 8. Omics-driven context-specific modelling

### RAVEN

| Function | Description |
|---|---|
| `tINIT` / `ftINIT` | INIT and fast-INIT (MILP) for tissue-specific models from HPA/GTEx expression data |
| `getINITModel2` | High-level wrapper loading HPA/GTEx data and calling ftINIT |
| `scoreComplexModel` | Scores gene-reaction-metabolite associations from expression data |
| `parseScores` | Converts raw expression matrix to per-reaction scores |

### COBRA

**Algorithm suite (all under `createTissueSpecificModel`)**

| Algorithm | Reference |
|---|---|
| GIMME | Becker & Palsson 2008 — minimises low-expression reactions while maintaining objective |
| iMAT | Zur et al. 2010 — balances inclusion of high- vs. exclusion of low-expression reactions |
| FASTCORE | Vlassis et al. 2014 — extracts minimal model containing a user-defined core |
| SWIFTCORE | Fastest consistency-based extraction |
| MBA | Model Building Algorithm — uses high/medium confidence reaction sets |
| eFlux | Constrains reaction bounds proportional to gene expression |
| MOOMIN | Metabolic model-based optimisation of omics integration |

**Multi-omics**

| Function | Description |
|---|---|
| `XomicsToModel` | Integrates transcriptomics + proteomics + metabolomics + literature into a thermodynamically consistent context-specific model |
| `XomicsToMultipleModels` | Batch generation of model variants with different parameter combinations |
| `modelExtraction` | Sub-step within XomicsToModel |

### Gaps

- **Multiple extraction algorithms**: RAVEN only implements INIT/ftINIT. For benchmarking
  or when HPA/GTEx data are unavailable (e.g., microorganisms), COBRA's suite (GIMME, iMAT,
  FASTCORE, etc.) offers alternatives.
- **Proteomics + metabolomics integration** in model extraction: `XomicsToModel` combines
  transcriptomics, proteomics, metabolomics, and literature curation in one pipeline.
  RAVEN's pipeline is transcriptomics-only.
- **Note**: For human tissue models where HPA/GTEx data are available, ftINIT is the
  current gold standard; COBRA's methods are not inherently superior here.

---

## 9. Thermodynamics

### RAVEN

None. RAVEN has no thermodynamic constraint methods.

### COBRA

**vonBertalanffy / thermoFBA module**

| Function | Description |
|---|---|
| `setupThermoModel` | Estimates standard transformed Gibbs energies (ΔrG'°) and sets directional constraints at in vivo conditions; requires ChemAxon/Open Babel |
| `initVonBertalanffy` | Initialises the module |
| `estimateDfGt0` / `estimateDrGt0` | Gibbs energies of formation/reaction with uncertainty bounds |
| `estimateDG_temp` | ΔG at specified temperature |
| `configureSetupThermoModelInputs` | Sets compartment pH, ionic strength, concentration bounds |
| `cycleFreeFlux` | Cycle-free (loopless) FBA solution |
| `addLoopLawConstraints` | Embeds thermodynamic loop constraints into LP/MILP |
| `consistentPotentials` | Chemical potentials satisfying thermodynamic directionality |
| `checkThermodynamicConsistency` | Validates thermodynamic feasibility |
| `fastSNP` | Generates minimal feasible basis for internal cycles (Saa & Nielsen 2016) |
| `findMinNull` | Minimal nullspace for internal cycles via MILP (Chan et al. 2017) |

### Gaps

RAVEN has no thermodynamic capabilities at all. If thermodynamic driving force analysis,
Max-min Driving Force (MDF), or thermodynamically feasible flux distributions are needed,
COBRA (or the separate matTFA toolbox) is required.

---

## 10. Flux sampling

### RAVEN

| Function | Description |
|---|---|
| `randomSampling` | Random-objective flux sampling — samples diverse flux states by solving FBA with random objective vectors |
| `analyzeSampling` | Correlates per-reaction expression t-scores with flux samples; returns enrichment scores |
| `getFluxZ` | Computes Z-scores comparing two flux sample sets |

### COBRA

| Function | Description |
|---|---|
| `sampleCbModel` | Unified entry point; dispatches to ACHR, CHRR, or RHMC |
| ACHR sampler | Artificial centering hit-and-run — ergodic, good convergence |
| CHRR sampler | Coordinate hit-and-run with rounding — uniform sampling from flux polytope |
| RHMC sampler | Riemannian manifold Hamiltonian Monte Carlo — fastest for high-dimensional models |
| `calcSampleStats` | Mean, SD, mode, median, skewness, kurtosis |
| `calcSampleDifference` | Pairwise differences between two sample sets |
| `compareSampleTraj` | Compares flux histograms across samples |
| `compareTwoSamplesStat` | KS test, rank-sum, chi-square, T-test between samples |
| `identifyCorrelSets` | Finds correlated reaction groups in flux samples |
| `plotHistConv` / `plotSampleHist` / `sampleScatterMatrix` | Visualisation and diagnostics |

### Gaps

- **Markov-chain sampling** (ACHR, CHRR, RHMC): these produce provably uniform samples
  from the feasible flux polytope. RAVEN's random-objective approach is fast but does not
  guarantee uniform coverage, and is known to under-sample extreme flux states.
- **Statistical comparison tools**: COBRA has built-in tests for comparing two sample sets;
  RAVEN's `getFluxZ` is a lightweight scoring function, not a hypothesis test.

---

## 11. Multi-species and microbiome modelling

### RAVEN

None.

### COBRA

**SteadyCom (community metabolic models)**

| Function | Description |
|---|---|
| `SteadyCom` | Finds maximum community growth rate at community steady-state |
| `SteadyComFVA` | FVA on community models across growth rates |
| `SteadyComPOA` | Pairwise analysis of member biomasses and reaction fluxes |

**Microbiome Modeling Toolbox / mgPipe (AGORA2 pipeline)**

| Function | Description |
|---|---|
| `mgPipe` / `initMgPipe` | Integrates metagenomic abundances with AGORA member models |
| `fastSetupCreator` | Builds microbiota community model with lumen/fecal compartments |
| `createPersonalizedModel` | Individual-specific community models from abundance data |
| `buildModelStorage` | Adds internal exchange space and coupling constraints |
| `createPanModels` | Pan-models at species/genus/family level |
| `microbiotaModelSimulator` | Applies diet constraints, runs FVA |
| `adaptVMHDietToAGORA` | Converts VMH diet definitions to AGORA-compatible constraints |
| `translateMetagenome2AGORA` | Maps metagenomic taxonomic IDs to AGORA model IDs |

**Whole-body physiological models (Harvey/Harvetta, Recon3D-based)**

| Function | Description |
|---|---|
| `optimizeWBModel` | LP optimisation on whole-body models |
| `Test4HumanFctExtv5` | Tests ~460 human metabolic functions |
| `checkIEM_WBM` | Simulates inborn errors of metabolism in whole-body models |
| `performIEMAnalysis` | Tests biomarker metabolites in biofluids for genetic disorders |
| `organEssentiality` | Determines organ essentiality |
| `calcOrganFract` / `calculateBMR` | Physiological parameterisation |

### Gaps

RAVEN has no multi-species or microbiome modelling capability. These are substantial,
specialised areas and their absence reflects RAVEN's focus on single-organism reconstruction.

---

## 12. Network topology and elementary flux modes

### RAVEN

| Function | Description |
|---|---|
| `getAllSubGraphs` | Identifies metabolite connectivity subgraphs (connected components) |
| `traceFluxPath` | Traces the carbon/flux path between two reactions via metabolite producers/consumers |
| `compareFluxes` | Diffs two flux distributions; reports turned-on, turned-off, flipped, and changed reactions |
| `followFluxes` | Interactive command-line exploration of a flux distribution |
| `followChanged` | Displays the reactions that changed most between two conditions |

### COBRA

| Function | Description |
|---|---|
| `findBlockedReaction` / `findBlockedIrrRxns` | Reactions with zero flux capacity |
| `findMetabolicJunctions` | Identifies metabolic branchpoints |
| `findCarbonRxns` | Finds reactions above a carbon count threshold |
| `conservedMoieties` | Computes conserved moieties (L-matrix) |
| `reactingMoieties` | Identifies reacting moiety pairs |
| `computeExtremeRays` | Computes extreme rays of the flux cone |
| EFMviz suite | Wrapper for EFMtool (Java): import, filter, yield analysis, enrichment analysis, backbone extraction, SBML export of EFM submodels |
| `surfNet` | Interactive network navigation |

### Gaps

- **Elementary flux modes (EFMs)**: COBRA interfaces with EFMtool to compute and analyse EFMs.
  RAVEN has no EFM support.
- **Conserved moieties**: Useful for detecting modelling errors and for metabolite grouping.
  Not present in RAVEN.
- **Extreme ray computation**: Not present in RAVEN.

### RAVEN strengths

- `traceFluxPath` and `compareFluxes` are practical diagnostic tools not present in COBRA;
  they make it easy to understand where flux flows and what changed between conditions.

---

## 13. Visualization and reporting

### RAVEN

| Function | Description |
|---|---|
| `drawFluxes` | Network diagram with reaction fluxes overlaid |
| `plotFluxes` | Bar/scatter plot of flux distributions |
| `printFluxes` | Formatted text table of a flux solution |
| `followFluxes` | CLI exploration of flux paths |
| Various plotting helpers in analysis/ | Phase planes, robustness plots, production envelopes |

### COBRA

| Function | Description |
|---|---|
| paint4Net | `draw_by_met`, `draw_by_rxn`, `bio_paint4net` — renders metabolic networks |
| Escher integration | Exports model/flux data to Escher JSON for interactive browser maps |
| CellDesigner conversion | Converts models to CellDesigner XML |
| EFMviz | Visualises EFMs as reaction/metabolite/gene networks |
| ReconMap | Interface to VMHDB/ReconMap web cartography |
| `metaboliteMassBalancePlot` | Bar chart of top producers/consumers of a metabolite |
| `upsetPlot` | Upset plots for set comparisons across models |
| `visualizeEpistasis` | Heatmap/network visualisation of double-deletion interactions |
| `plotEssentialRxns` | Heatmap of essential reactions across models |
| `phenotypePhasePlane` | 3D surface + shadow prices |
| `plotSampleHist` / `sampleScatterMatrix` | Sampling diagnostics |

### Gaps

- **Escher export**: COBRA can export directly to Escher, the standard web-based metabolic
  network visualisation. RAVEN does not.
- **EFMviz**: Visualising EFM composition, backbone networks, and enrichment.
- **Upset plots / epistasis heatmaps**: Multi-model comparison and double-deletion visualisation.

---

## 14. Metabolomics integration

### RAVEN

Limited to the production/uptake exchange bounds and ftINIT's optional metabolomics scoring
block (which is not yet fully implemented in all pathways).

### COBRA — metabotools

| Function | Description |
|---|---|
| `prepIntegrationQuant` | Generates uptake/secretion profiles from extracellular metabolomics flux data |
| `defineUptakeSecretionProfiles` | Calculates slope ratios for comparative conditions |
| `setMediumConstraints` | Applies medium composition constraints |
| `setQualitativeConstraints` | Applies detection-limit constraints from MS data |
| `setSemiQuantConstraints` | Relative difference-based constraints between two conditions |
| `setQuantConstraints` | Full quantitative metabolomics data integration as model constraints |
| `setConstraintsOnBiomassReaction` | Applies growth rate constraints from doubling times |
| `calculateLODs` | Converts detection limits (ng/mL → mM) |
| `predictFluxSplits` / `computeFluxSplits` | Pathway contribution analysis for a metabolite |
| `MetaboAnnotator` | Annotates metabolites via structure, InChI, VMH, PubChem, HMDB |
| `conc2Rate` | Converts concentrations to uptake rates from cell parameters |

### Gaps

The metabotools suite is not present in RAVEN. For constraining models with extracellular
metabolomics (uptake/secretion rates from NMR or LC-MS), COBRA provides a complete pipeline.

---

## 15. Model I/O and database connectivity

### RAVEN

| Area | RAVEN |
|---|---|
| SBML | Read/write SBML L2V1 and L3 (via `importModel` / `exportModel`) |
| Excel | Structured Excel format for model editing (`importExcelModel` / `exportToExcelFormat`) |
| JSON | COBRA-compatible JSON export |
| YAML | RAVEN YAML format (human-readable, Git-friendly) |
| KEGG | Full REST API integration: pathway, reaction, compound, enzyme queries |
| UniProt | Gene/protein queries for annotation enrichment |
| HPA / GTEx | Direct integration for tINIT expression data loading |
| Ensembl | Gene ID mapping |
| Git-stable export | `exportForGit` produces deterministic field ordering for clean diffs |

### COBRA

| Area | COBRA |
|---|---|
| SBML | Read/write via libSBML (more complete SBML compliance including FBC v2) |
| Excel | Read/write |
| JSON | BIGG-compatible JSON |
| BiGG | BiGG database REST API; metabolite/reaction annotation lookup |
| VMH | Virtual Metabolic Human database integration |
| AGORA2 | Direct interface to 7302 microbial GEMs |
| Recon3D | Integration with the human metabolic reconstruction |

### Gaps

- **SBML compliance**: COBRA's libSBML wrapper handles FBC v2 and Groups extensions more fully.
- **BiGG / VMH databases**: COBRA is tightly coupled to these namespaces; RAVEN uses KEGG and
  does not have native BiGG/VMH lookup tools.
- **AGORA2 access**: Not in RAVEN.

### RAVEN strengths

- **YAML format** and **`exportForGit`**: COBRA has no equivalent to RAVEN's Git-friendly export.
  The YAML format produces minimal, human-readable diffs — important for collaborative GEM
  development in version control.
- **Excel-to-model pipeline**: RAVEN's Excel model format is mature and widely used in the
  community for manually curated models.
- **KEGG integration**: RAVEN's KEGG functions are more comprehensive than COBRA's.

---

## 16. Summary gap table

Functions and capabilities present in COBRA but absent or minimal in RAVEN, grouped by priority
for a GEM reconstruction lab.

### High relevance for reconstruction labs

| Capability | COBRA functions | RAVEN status |
|---|---|---|
| MILP gap-filling | `fastGapFill`, `swiftGapFill` | Only heuristic LP |
| Stoichiometric consistency | `fastcc`, `fastLeakTest`, `checkDuplicateRxn`, `detectDeadEnds` | Struct validation only |
| Multiple omics extraction algorithms | GIMME, iMAT, FASTCORE, MBA, SWIFTCORE, eFlux | ftINIT only |
| Multi-omics model extraction | `XomicsToModel` | Not present |
| Metabolomics constraints | Full metabotools suite | Not present |
| Double gene deletion | `doubleGeneDeletion` | Single KO only |
| pFBA (gene classification) | `pFBA` | Not present |
| Loopless FBA | `cycleFreeFlux`, `addLoopLawConstraints` | Not present |

### Medium relevance

| Capability | COBRA functions | RAVEN status |
|---|---|---|
| Strain design (bilevel) | `OptKnock`, `optGene`, `GDLS` | Simplified approximation only |
| OptForce (must sets) | `findMustU/L/UU/UL/LL`, `optForce` | Not present |
| Genetic MCS | `calculateGeneMCS`, `calculateMCS` | Not present |
| MOMA / ROOM knockout prediction | `MOMA`, `ROOM`, `linearMOMA`, `linearROOM` | Not present |
| Markov-chain flux sampling | ACHR, CHRR, RHMC | Random-objective only |
| Thermodynamic constraints | vonBertalanffy, thermoFBA | Not present |
| EFM analysis | EFMviz + EFMtool | Not present |
| Conserved moieties | `conservedMoieties` | Not present |

### Specialised / niche

| Capability | COBRA functions | RAVEN status |
|---|---|---|
| Microbiome community models | mgPipe, SteadyCom | Not present |
| Whole-body models | Harvey/Harvetta toolbox | Not present |
| Synthetic lethality | Fast-SL | Not present |
| Thermodynamic model setup | `setupThermoModel` (vonBertalanffy) | Not present |
| Escher / CellDesigner export | Escher JSON, CellDesigner XML | Not present |

---

## 17. RAVEN strengths not in COBRA

These are areas where RAVEN has functionality that COBRA lacks:

| Capability | RAVEN functions | COBRA status |
|---|---|---|
| KEGG-based de-novo reconstruction | `getKEGGModelForOrganism` | No equivalent pipeline |
| HMM-based gene calling | `getModelFromHomology` (HMMER/BLAST) | BLAST only via DEMETER |
| ftINIT (tissue-specific, HPA/GTEx) | `ftINIT`, `getINITModel2` | Older methods only |
| Git-stable model export | `exportForGit` | Not present |
| YAML model format | RAVEN YAML | Not present |
| Flux path tracing | `traceFluxPath` | Not present |
| Flux diff between conditions | `compareFluxes` | Not present |
| Minimal medium (MILP) | `getMinimalMedium` | `findMinimalMedia` (similar) |
| Minimal flux (MILP, per-reaction) | `getMinNrFluxes` | `sparseFBA` (similar) |
| Reporter metabolites | `reporterMetabolites` | Not present |
| FSEOF | `FSEOF` | Not present |
| KEGG database queries | `getKEGGComps`, `getKEGGRxns`, `getRxnsFromKEGG`, etc. | Limited |

---

*Last updated: 2026-06-19. Based on RAVEN v2.x (develop3 branch) and COBRA Toolbox v3.x.*
