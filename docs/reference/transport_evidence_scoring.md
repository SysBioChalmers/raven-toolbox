# Evidence-aware transport scoring — design & implementation plan

A plan to replace the **blanket** inter-compartment transport penalty in the localisation assignment
with a **per-transport, evidence-aware** cost, in both **RAVEN** (MATLAB) and **raven-toolbox**
(Python). Carrier-general (any transporter family, any membrane) and organism-agnostic (sequence-only).

## Motivation

The assignment step — RAVEN [`predictLocalization`](https://github.com/SysBioChalmers/RAVEN/blob/main/localization/predictLocalization.m)
and raven-toolbox `predict_localization` — adds a transport reaction whenever a gene is placed away
from a metabolite it must reach, charging a fixed `transportCost`/`transport_cost` per transport.
A blanket
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
and ChEBI-ontology tables) and the pipeline is **validated** on yeast-GEM (see *Validation* below).

## Validation

The pipeline is validated on real curated ground truth using only the shipped public API — no
benchmark-specific scaffolding — in two studies under `docs/studies/`:

* [**Replicate yeast-GEM**](../studies/replicate_yeast_gem.md) — flatten curated yeast-GEM to one
  compartment with `merge_compartments`, then run `annotate_proteome` -> `evidence_aware_transport_cost`
  -> `assign_compartments` and score the recovered compartmentalisation (reaction- and gene-level
  agreement, added-transport count, functional connectivity, growth) against the curated original.
* [**vs CarveFungi**](../studies/vs_carvefungi.md) — the same evidence-aware assignment head-to-head
  with CarveFungi's own compartment prediction, scored against yeast-GEM at gene level and against
  CarveFungi's own placement.

Both drive the real functions on the yeast proteome + DeepLoc inputs in `data/deeploc/`; see each
study for the current numbers.

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
