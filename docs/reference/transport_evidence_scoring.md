# Evidence-aware transport scoring — design & implementation plan

A plan to replace the **blanket** inter-compartment transport penalty in the localisation assignment
with a **per-transport, evidence-aware** cost, in both **RAVEN** (MATLAB) and **raven-toolbox**
(Python). Carrier-general (any transporter family, any membrane) and organism-agnostic (sequence-only).

## Motivation

The assignment step — RAVEN [`predictLocalization`](https://github.com/SysBioChalmers/RAVEN/blob/main/localization/predictLocalization.m)
and raven-toolbox `predict_localization` — adds a transport reaction whenever a gene is placed away
from a metabolite it must reach, charging a fixed `transportCost`/`transport_cost` per transport.
The [CarveFungi head-to-head study](../studies/carvefungi_milp_benchmark.md) showed this blanket
penalty is **indiscriminate**: it drops curated, functionally essential transporters (malate–aspartate
and citrate shuttles, CoA-precursor and NADPH carriers — 5 individually essential in yeast-GEM) at the
same rate as spurious ones, because the cost ignores whether a real transporter exists. The fix is to
let **transporter-level sequence evidence** modulate the cost: cheap for supported transports, full
prior for unsupported ones.

## Key design decisions

1. **Local binaries, not a web API** (see the dedicated answer below). Reuse the `hmmsearch` and
   `diamond` binaries **both repos already bundle** — no new executable. Only two small reference
   databases are new.
2. **No change to the MILP.** Both solvers already accept a per-metabolite transport cost
   (`predictLocalization`: a `transportCost` vector of length `numel(model.mets)`, negative encourages;
   `predict_localization`: `transport_cost: float | Mapping[str, float]`). The new functionality is an
   **evidence → cost** layer that produces that vector/mapping; the assignment code is untouched.
3. **Predictor- and organism-agnostic.** Membrane placement reuses the **DeepLoc compartment output
   already in the pipeline** (the reliable organelle calls, trust 0.78–0.88 — not the noisy
   membrane-*type* output). Every evidence source is a sequence search, HMM, or orthology lookup with
   no species-specific tables, so the only per-organism input is the proteome FASTA.
4. **A shared intermediate contract.** A per-gene *transporter-annotation table* decouples the
   annotation step (binaries) from the scoring step (MILP). This enables caching, a **bring-your-own-
   annotation** mode (feed a table from any tool, incl. web services), and identical semantics across
   the two implementations.

## Architecture (shared across RAVEN and raven-toolbox)

```
proteome FASTA ─┬─ hmmsearch  vs Pfam transporter HMMs ─┐
                └─ diamond    vs TCDB sequences         ├─► per-gene transporter annotation
DeepLoc output (already in pipeline) ───────────────────┘   {gene: is_carrier·conf, family,
                                                              substrate class, TC mechanism, compartment}
                                                                          │
   per candidate transport t (metabolite m across membrane M={c1,c2}):   ▼
   evidence(t) = max over genes g of  conf_transporter(g)·compartment_match(g,M)·substrate_match(g,m)
   transport_cost(t) = base_cost · (1 − evidence(t))     # supported → cheap; unsupported → full prior
                                                                          │
                                                                          ▼
                       existing MILP (predictLocalization / predict_localization), unchanged
```

* **Stage 1 — Annotate.** `hmmsearch` against a transporter-family HMM set (MCF `PF00153`/SLC25, MFS,
  ABC, amino-acid/sugar permeases, aquaporins, P-type ATPases, …) gives "is a carrier" + coarse
  substrate class; `diamond blastp` against TCDB gives a TC number → substrate class **and** mechanism
  (uni/sym/antiport). Output: the per-gene annotation table.
* **Stage 2 — Place.** Join the annotation with the DeepLoc compartment per gene: a carrier gene
  predicted in compartment *X* supports transports across *X*'s boundary (`compartment_match`).
* **Stage 3 — Score.** Compute `evidence(t)` and `transport_cost(t)` per transport; `conf_transporter`
  from Pfam/TCDB hit strength, `substrate_match` from the TC/family substrate class vs *m*'s class.
* **Stage 4 — Solve.** Feed the per-transport costs to the existing assignment MILP.

## Binaries vs API — the explicit answer

**Primarily local binaries, reusing what both repos already provision; an API is only an optional
fallback.**

* **Local (default).** Annotation is proteome-scale (thousands of proteins). The two tools needed —
  `hmmsearch` (HMMER) and `diamond` — are **already bundled and invoked** in both repos:
  * raven-toolbox: `src/raven_toolbox/binaries.py` provisions `diamond` + `hmmsearch` (raven-data
    release ZIPs), used by `reconstruction/kegg/hmm.py` and `reconstruction/homology/blast.py`.
  * RAVEN: `downloadRavenBinaries.m` bundles DIAMOND + HMMER; `getDiamond.m` and
    `getKEGGModelForOrganism.m` already shell out to them.
  So **no new binary is required** — only two new reference databases (below). Local runs are the right
  default for reproducibility and for RAVEN's offline/HPC use; web APIs (EBI InterProScan, TCDB web
  BLAST) are rate-limited and impractical at proteome scale.
* **DeepLoc is unchanged.** It is already an external step the user runs (web server *or* local
  install) and whose output the localisation module consumes. Transport scoring **reuses that existing
  output** — it adds no new prediction dependency.
* **Optional online / no-binary paths.** (a) **Bring-your-own annotation:** accept a pre-computed
  transporter-annotation table from any source (eggNOG-mapper, InterProScan, a web service), bypassing
  the binaries entirely — this also preserves the module's predictor-agnostic philosophy. (b) A small
  **InterProScan/TCDB REST fallback** for tiny inputs or environments without the binaries.

## New reference data (both repos)

Provisioned through the existing data channels (raven-toolbox `data.py` / raven-data release assets;
RAVEN `downloadRavenBinaries.m` data + `checkInstallation.m`):

* **Transporter HMMs** — a curated set of transporter-family Pfam accessions compiled into one `.hmm`
  (the accession list is the only thing to maintain; `hmmpress` it for `hmmsearch`).
* **TCDB sequences** — the TCDB FASTA (small, ~tens of MB) built into a DIAMOND db, plus the
  TC-number → substrate-class table.
* **Substrate-class ontology** — a small, bundled mapping from metabolite identity (name/ChEBI/KEGG)
  and from TC/family substrate descriptors to a coarse shared class (sugars / amino acids / organic
  acids / ions / nucleotides / lipids). This is the join key for `substrate_match`.

## raven-toolbox (Python) plan

New module **`src/raven_toolbox/localization/transport_evidence.py`**:

```python
@dataclass
class TransporterAnnotation:
    gene: str
    confidence: float            # max of Pfam/TCDB hit strengths, 0..1
    families: list[str]          # e.g. ["PF00153"]
    tc_numbers: list[str]        # e.g. ["2.A.29.2.4"]
    substrate_classes: set[str]  # coarse classes
    mechanism: str | None        # "uniport"|"symport"|"antiport"|None

def annotate_transporters(
    proteome: str | Path,                 # FASTA
    *, hmm_db: Path | None = None,        # default: provisioned transporter HMMs
    tcdb_db: Path | None = None,          # default: provisioned TCDB DIAMOND db
    table: pd.DataFrame | None = None,    # bring-your-own: skip the binaries
    threads: int = 1,
) -> dict[str, TransporterAnnotation]: ...

def evidence_aware_transport_cost(
    model: cobra.Model,
    annotation: Mapping[str, TransporterAnnotation],
    gene_compartments: Mapping[str, set[str]],   # from the existing DeepLoc scores
    *, base_cost: float = 0.5, unsupported_floor: float = 0.5,
) -> dict[str, float]:                            # {metabolite_base: cost} for predict_localization
    ...
```

* Reuse `binaries.py` to resolve `hmmsearch`/`diamond`; mirror the invocation patterns in
  `reconstruction/kegg/hmm.py` and `reconstruction/homology/blast.py`. (`pyhmmer` is a viable
  pip-only alternative to shelling out, if we prefer no external HMMER for the Python side.)
* Integration: `predict_localization` already accepts the resulting mapping via `transport_cost=`.
  Optionally add a convenience `transport_evidence=` parameter that runs the layer internally.
* **Extension:** the current cost mapping is keyed by metabolite base; add optional keying by
  `(metabolite_base, frozenset(compartment_pair))` so membrane-specific evidence is honoured.
* Tests: `tests/test_transport_evidence.py` — annotation parsing, `compartment_match`/`substrate_match`
  logic, and an end-to-end `predict_localization` run where a supported transport is retained and an
  unsupported one is dropped. Gate binary-dependent tests with `shutil.which` / `pytest.importorskip`.

## RAVEN (MATLAB) plan

New functions under `localization/` (or a new `transporters/`), reusing the homology/KEGG binary
wrappers:

```matlab
annotation   = annotateTransporters(fastaFile, varargin)
%   wraps getDiamond-style DIAMOND vs TCDB + hmmsearch vs Pfam (the getKEGGModelForOrganism HMMER
%   pattern); name-value: 'hmmDB','tcdbDB','table' (bring-your-own), 'cores'

transportCost = scoreTransportEvidence(model, annotation, geneComps, varargin)
%   returns a transportCost vector (numel(model.mets)) to pass straight to predictLocalization;
%   geneComps from the existing DeepLoc-based localisation
```

* Reuse `downloadRavenBinaries.m` (DIAMOND + HMMER already bundled) and add the two databases to the
  data provisioning + `checkInstallation.m`.
* **COBRA name-collision check:** RAVEN `.m` filenames must not clash with COBRA Toolbox function
  names — verify `annotateTransporters`/`scoreTransportEvidence` (and any helper) are unique before
  committing.
* Integration: pass the vector to the existing `predictLocalization(... ,'transportCost', v)`. If
  membrane-specific costs are wanted, extend `predictLocalization` to accept a per-(met,compartment)
  cost (currently per-met) — optional, additive.
* Tests in `testing/` following the existing transporter/localisation test pattern; gate on
  `tBinaries`-style binary availability.

## Phased rollout

1. **Family scan + DeepLoc placement** — `hmmsearch` over the transporter HMM set + compartment from
   DeepLoc. Carrier-general and all-membrane from day one (the dropped transports span c↔mito 27,
   c↔extracellular 20, c↔ER 7, c↔peroxisome 7 — no membrane is privileged).
2. **TCDB substrate specificity + mechanism** — `diamond` vs TCDB → substrate-matched scoring.
3. **Consensus / refinement** — combine family + TCDB + orthology (EggNOG/KEGG, already computed in
   reconstruction) + DeepLoc; add directionality; resolve conflicts.

**Status.** The **coarse-first pipeline is implemented** in
`raven_toolbox.localization.transport_evidence`:

* `evidence_aware_transport_cost` — the scoring core (per-metabolite `transport_cost` mapping both
  assignment functions already accept).
* `annotate_proteome` — the **`hmmsearch` (Pfam families) + `diamond` (TCDB) back-end**: scans a
  proteome FASTA against the transporter Pfam HMM db and the TCDB DIAMOND db (both auto-downloaded from
  the raven-data `transporters-*` release), mapping families/TC-numbers to coarse substrate classes via
  the curated `transporter_tables`. `annotate_transporters` still takes a pre-computed table.
* `default_substrate_of` — the **model-side** coarse classifier (metabolite name → substrate class),
  so a metabolite and a transporter meet in the shared vocabulary.
* `SubstrateOntology` + `substrate_chebi` — the **specific-substrate layer**: TCDB's curated
  `TC-ID → substrate ChEBI` table + the ChEBI `is_a`/protonation graph give a graded
  metabolite→substrate roll-up (exact 1.0, decaying by hop; an optional `sibling_weight` also credits
  chemical *relatives* of the cargo) that `evidence_aware_transport_cost` layers on top of the coarse
  class; `annotate_proteome` fills `TransporterAnnotation.substrate_chebi`.

Databases are built by `scripts/build_transporter_data.py` (Pfam HMMs + TCDB DB + the TCDB-substrate
and ChEBI-ontology tables) and the pipeline is **validated** on yeast-GEM (see *Validation → Result*).

## Validation

Reuse this study's own benchmark
([`analyse_carvefungi_transports.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/scripts/analyse_carvefungi_transports.py)):
curated-transport precision/recall + the functional (essentiality) test, before vs after. Success:

* the **kept**-transport curated-match rate rises *above* the dropped rate (the cut becomes
  *selective*, no longer ~equal at 41% vs 42%);
* the 5 individually-essential transports are retained;
* the gains **reproduce on a non-fungal model** (e.g. AraCore) — the organism-agnosticism check.

**Result (yeast-GEM,
[`validate_transport_evidence.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/scripts/validate_transport_evidence.py)).**
Annotating the yeast proteome finds 364 transporter genes; scoring all 1337 metabolites against them
(`base_cost` 0.5, carrier compartments from DeepLoc) makes the cut **selective**: of the metabolites
the evidence keeps cheap (387), **70 %** are curated yeast-GEM transports, versus only **39 %** among
the unsupported set the full prior would drop — where a blanket penalty treats both alike at the 48 %
base rate. Every essential cytosol↔mito shuttle checked — (S)-malate, citrate, 2-oxoglutarate,
oxaloacetate, NADP(+)/NADPH, L-serine, 2-dehydropantoate, pyruvate, PEP — is fully evidenced (cost 0,
retained). Coarse-stage recall is 42 % — but most of that miss is **gene-detection / localisation**
bound rather than substrate-bound (see the ChEBI-layer note). The non-fungal reproduction (AraCore)
remains to run.

**Selectivity on yeast-GEM.** On yeast-GEM (a larger ground truth than the 138-transport carve), the
same script measures how selective the cut is — the curated rate among *kept* (would-retain) vs
*dropped* (would-drop) metabolites (base rate 48 %):

| approach | kept curated | dropped curated | recall | selective? |
|---|--:|--:|--:|:--|
| no evidence (base rate) | 48 % | 48 % | — | no (indiscriminate) |
| evidence: coarse | 70 % | 39 % | 42 % | yes |
| evidence: + ChEBI | 70 % | 37 % | 47 % | yes |
| evidence: + ChEBI + `sibling_weight` 0.5 | 52 % | 38 % | 75 % | recall-tilted |

Without transporter evidence the kept and dropped sets are both ~48 % curated (indiscriminate).
Evidence weighting makes the kept set ~1.8× more curated than the dropped set and retains **all** the
individually essential carriers (cost 0). `sibling_weight` trades specificity for recall (precision
70 → 52 %, recall 47 → 75 %) and is off by default.

**On the actual carve — feasibility-respecting.** `analyse_carvefungi_transports.py` benchmarks against
a *genuine*, functional, gene-annotated reconstruction (`build_reference_carve_model.py` drives
CarveFungi's own EggNOG scoring + its `carve_model` MILP + its gene-annotation step end to end — not a
bare reaction-id cache, which drops CarveFungi's uptake reactions and each reaction's solved direction
and so cannot grow standalone or support a real feasibility check). The regenerated model: 991
reactions, 280 genes, 591 GPR-annotated, growth 0.758 on a defined minimal medium (9.7 % MIP gap).

On its 170 inter-compartment transports (58 match curated yeast-GEM, 11 individually essential), "ours"
is a **feasibility-respecting** reduction: rank unsupported transports (cost ≥ 0.35) worst-evidence
first, tentatively knock each one out (bounds → 0,0) and re-run FBA, keep it knocked out if growth
survives (`--min-growth-fraction`, default **0.9** of native), otherwise restore it — feasibility
overriding missing evidence, one reaction at a time (a greedy upper bound on what a joint solve could
drop; see the caveats below). **Why 0.9, not "still alive":** an earlier pass used a 1 %-of-native
floor (only "not dead") and produced a symptom worth naming explicitly — achieved growth *fell* as
`sibling_weight` rose (0.168 → 0.019), even though evidence can only ever *add* credit, never remove it
(`evidence_aware_transport_cost` combines coarse/ChEBI/sibling with `max`, and the set of
evidence-protected transports was verified strictly non-shrinking as sibling weight increases). That
was real, but it was the *greedy search's* artefact, not the evidence's: a 1 % floor lets a long chain
of individually-small hits compound into severe, undetected growth erosion, and *which* chain gets
tried first depends on the exact cost values, which shift across variants — different variants sacrifice
different (but each individually replaceable) transports on the way to the same weak bar. Raising the
floor to 0.9 makes the reduction growth-*preserving*, not just growth-nonzero, and removes that
path-dependent noise from the reported growth (swept in the reference commit: at floor ≥ 0.7 the
achieved growth is effectively equal — within numerical noise — across every variant tested):

| approach | transports kept | curated replicated | essential kept | spurious kept | growth |
|---|--:|--:|--:|--:|--:|
| CarveFungi (native) | 170 | 58/58 | 11/11 | 112 | 0.758 |
| ours: coarse | 110 | 56/58 | **11/11** | 54 | 0.752 |
| ours: + ChEBI | 111 | 56/58 | 10/11 | 55 | 0.752 |
| ours: + ChEBI + sibling 0.3–0.7 | 110–113 | 56–57/58 | 10/11 | 54–56 | 0.752 |
| ours: + ChEBI + sibling 1.0 | 114 | **58/58** | **11/11** | 56 | 0.752 |

CarveFungi's native carve is **bloated** — 112 of its 170 transports are spurious (non-curated), because
it never minimises transport. Every evidence-aware variant cuts the transport network by a third
(170→110–114) and spurious transports by **~50 %** (112→54–56) while keeping **97–100 %** of curated
and **91–100 %** of individually-essential transports, at **99 % of native growth** (0.752 / 0.758) in
every variant — full curated+essential retention at both the `coarse` baseline and at `sibling` 1.0.

**On growth units — a caveat, and the right standard to hold this to.** These "growth" values are the
regenerated carve's *own* FBA objective on its artificially-constructed biomass reaction and its
hand-set minimal medium — not a calibrated growth **rate** in h⁻¹. yeast-GEM's own growth-rate
validation
([`growth.py`](https://github.com/SysBioChalmers/yeast-GEM/blob/0b717e7dd5ca8a3b1b074f8055a736c2e9ec33ee/code/python/yeastgem/model_tests/growth.py))
is the right standard to aspire to: it fixes *experimentally measured* chemostat uptake rates (glucose/
O₂/NH₃) reaction-by-reaction and reports R² against 32 real growth-rate observations. It does not
transfer directly to this carve (different reaction-id namespace, a `carve_model`-derived biomass
equation with hand-picked stoichiometric weights, no calibrated exchange bounds) — so the only claim
made here is the scale-invariant one, fraction of *this same model's own* native optimum, which is
exactly what both the 1 %-floor symptom and its 0.9-floor fix are about.

**Remaining caveat — greedy path-dependence.** This is a one-at-a-time, worst-evidence-first reduction,
not a joint MILP, so it is an upper bound on what feasibility alone would allow to drop, and *which*
essential transport survives can still depend on removal order even at a growth-preserving floor — one
essential transport (an ergosterol-precursor cytosol↔ER shuttle, `r_1754`) is retained at `coarse` and
at `sibling` 1.0 but dropped in between, because enabling ChEBI shifts *other* transports' relative
costs enough to change the removal order, and carve-local redundancy lets the greedy pass sacrifice it
under some orderings but not others. A joint solve would not have this ambiguity; at the 0.9 floor its
effect is confined to this single essential/curated count, not to growth. Reproduce:
`build_reference_carve_model.py --carvefungi-dir <clone> --out <path>` (once; ~20 min CPLEX carve) then
`analyse_carvefungi_transports.py --model <path> --yeast-gem <yeast-GEM>` (`--min-growth-fraction` to
sweep the floor).

**ChEBI layer (yeast-GEM).** Adding the graded ChEBI roll-up on top of the coarse class lifts the
selective cut's **recall from 42 % to 47 %** (kept 387→426) at steady 70 % precision. Two details make
it work: the roll-up is graded by hop distance (a metabolite that *is* the curated substrate scores
1.0, a subtype or protonation/tautomer variant less), and — essential — TCDB annotates substrates with
*secondary/deprecated* ChEBI ids that carry no `is_a` edges of their own, so they are **normalised onto
their connected primary id** (~19 k `alt_id` mappings shipped in `chebi_relations.tsv.gz`) before any
walk; without that, even glucose↔glucose scored 0 and only 17 % of assigned substrates reached a model
metabolite (51 % after). The remaining miss is **gene-detection / membrane-localisation** bound, not
substrate-bound: ~58 % of curated transports have no detected, correctly-localised carrier of matching
cargo, so substrate precision cannot conjure one — forcing the ChEBI verdict (`strict`, no coarse
rescue) *cuts* recall to 31 % without raising precision, confirming the coarse classes already fit
yeast's promiscuous MFS/MCF/ABC families. The layer is kept because it never regresses
(strongest-evidence-wins) and it exposes each carrier's *specific* substrate ChEBIs — decisive for
models with **narrow-specificity transporters**, where a coarse class collapses cargo the roll-up keeps
distinct (a curated *hexose* carrier scores D-glucose 0.80, D-galactose 0.90, D-fructose/D-mannose 0.65
— all of which one "sugar" class treats alike).

**Substrate discrimination (intrinsic, model-free —
[`analyse_substrate_discrimination.py`](https://github.com/SysBioChalmers/raven-toolbox/blob/develop/scripts/analyse_substrate_discrimination.py)).**
Where that payoff comes from, isolated from the detection/localisation confound: within each coarse
TCDB family, sibling transport systems carry *different* specific substrates. Over every (curated
substrate, same-family system) pair, the ChEBI roll-up **rules out ~99 %** of the non-carriers —
unrelated cargo scores 0, and even the chemically-related siblings score only ~0.78 (vs 1.00 for a
true carrier) — whereas the coarse class, shared across the whole family, rules out **0 %** (it scores
every sibling identically). So the layer supplies exactly the substrate resolution a coarse class
cannot: latent on yeast's broad MFS/MCF/ABC families, decisive on a proteome of narrow-specificity
transporters (specific sugar/amino-acid permeases, ion channels).

## Open questions / risks

* **Substrate matching** (metabolite → coarse class, and TC/family descriptor → the same class) is the
  hardest piece; start coarse and expand. A wrong match wrongly cheapens a transport.
* **Transporter HMM accession list** needs curation and periodic refresh (Pfam releases).
* **Annotation completeness varies by organism** — keep the unsupported-transport prior *mild and
  tunable*; never hard-forbid a transport for lack of evidence.
* **Per-(met, compartment-pair) cost keying** is an additive extension on both sides; scalar/per-met
  works first.
* **Binary availability:** `hmmsearch` + `diamond` are bundled cross-platform in both repos; the
  bring-your-own-annotation mode covers any environment where they are not.
