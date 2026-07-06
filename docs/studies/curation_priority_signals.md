<!-- Design catalog: signals for ranking which compartment assignments / transports need manual
curation. Generated from a multi-lens brainstorm + synthesis; verify computability against the code
before implementing. Intended to back a curation_priority() feature on AssignmentProposal. -->

# Curation-Priority Catalog for `assign_compartments`

## 1. What a curation-priority score is

A **curation-priority score** is a single 0–1 "please-check-me" number attached to every reaction placement and every added transport in an `AssignmentProposal`, so a curator can sort ~thousands of calls and review only the highest-value top few percent. Each signal below is a detector that returns, per target, a raw 0–1 suspicion value plus a human-readable reason. The composite is a **weighted, saturating combination** — not a plain sum. Two families feed it: *evidence-uncertainty* signals (was the localization call itself shaky?) and *topology/functionality* signals (does the placed network betray a mislocation — islands, transport sandwiches, blocked flux, fragile certificates?). The multiplicative "high-stakes × uncertain" gates (essential **and** evidence-weak; override **and** prior-violating) are what separate a merely-uncertain call from a *consequential* one. Compute cheap signals over all targets first; run the expensive FBA/FVA signals only on the candidates the cheap layer already flagged, then blend into a final rank and surface a ranked table with per-flag explanations and a suggested curator action.

---

## 2. Deduplicated signal catalog

Merges applied: the seed "reaction between two transports" appears four times in the brainstorm (`reaction_flanked_by_added_transports`, `transport_sandwich`, `reaction_between_two_transports`, `reaction_stranded_between_transports`) → merged into one canonical **Transport Sandwich** signal plus its transport-side twin **Removable-by-Replacement**. `functionality_override_on_weak_evidence` + `functionality_override` → **Override Against Confident Evidence**. Cross-species signals (×3) → **Ortholog Compartment Outlier**. Multi-loc duplicate signals (`inactive_multiloc_duplicate`, `mono_multi_flip`, `futile_cycle_only_flux`) → **Dead Multi-Loc Duplicate**. Certification single-point-of-failure signals (`certification_sensitive_singleton`, `sole_lifeline_transport`, `certification_floor_fragility`) → **Fragile Certificate Load-Bearer**. Impermeant-cargo signals (`implausible_cargo_chemistry`, `transport_creates_biologically_absent_pool`) → **Impermeant Cargo**.

### Lens A — GPR evidence disagreement (cheap, no FBA)

**A1. Isozyme Compartment Split** · *target: placement*
The reaction's OR-linked isozymes point to different compartments; the placement silently kept one and discarded the losers. Real dual-localization looks identical to a wrong pick.
- **Fires when:** ≥2 OR-linked genes each have a scored row and their argmax compartments partition into ≥2 distinct compartments, each clearing a confidence floor (e.g. ≥0.5).
- **Compute:** Parse GPR into isozyme clauses; per clause take argmax compartment; fire if the distinct-compartment set has size ≥2. Score = `1 − (evidence for assigned comp / total isozyme evidence)`. Down-weight if the reaction was multi-localized (both kept).
- **Cost:** cheap · **FP risk:** moderate (real duals — e.g. cytosolic+mito aminoacyl-tRNA synthetases). · **Action:** open the GPR; confirm whether the discarded isozyme's compartment is the correct one, or whether the placement should be multi-localized.

**A2. Complex Subunit Disagreement** · *target: placement*
AND-linked subunits of one physical complex must assemble in one compartment; if their DeepLoc calls disagree, one call is wrong (stronger error signal than A1, since subunits can't function apart).
- **Fires when:** within one AND-clause, ≥2 subunit genes have scored rows whose argmax compartments differ.
- **Compute:** per AND-clause, gather subunit scores; fire on argmax disagreement. Priority = `1 − mean(subunit score in assigned comp)`. Weight down a dissenter with low `raw_confidence` (predictor noise).
- **Cost:** cheap · **FP risk:** low-moderate (loose "complexes" that aren't obligate; moonlighting subunits). · **Action:** verify co-assembly; check the dissenting subunit against UniProt.

**A3. Solver Tie-Break** · *target: gene*
When a gene's top-two compartment scores are near-equal, the MILP's pick is arbitrary (can flip run-to-run under warm start). Essentially no evidence-based preference.
- **Fires when:** `score[0] − score[1] < 0.05` **and** the two compartments differ in consequence (not two cytosol-like comps) **and** the gated reactions carry flux.
- **Compute:** sort each gene's score row; flag tight margins. Rank by inverse margin × flux-importance of gated reactions.
- **Cost:** cheap · **FP risk:** moderate (harmless if both comps are equivalent). · **Action:** resolve with a second predictor / UniProt — any weak extra evidence fixes it deterministically.

### Lens B — Evidence strength & predictor capability (cheap)

**B1. Override Against Confident Evidence** · *target: placement* — **highest-value evidence class**
The method's headline move: placing a reaction *against* its gene's top DeepLoc compartment because flux demanded it. Trustworthy when it overrode weak evidence; suspicious when it overrode a **sharp, confident** call — strong localization evidence was sacrificed to make the network grow.
- **Fires when:** assigned compartment ≠ gene-consensus argmax `c_top`, **and** `c_top` was high-confidence (high `raw_confidence`, wide top-1−top-2 margin).
- **Compute:** per reaction, `c_top`, its `raw_confidence`, and margin. Priority = `raw_confidence(c_top) × margin(c_top) × (1 − score(assigned_comp))`. Confirm with the method's internal functionality-override flag if exposed. Rank by `top_score − placed_score`.
- **Cost:** cheap · **FP risk:** moderate (legit retro-translocation / dual targeting DeepLoc can't see — e.g. FUM1, ACO1). Mitigate by combining with a topology flag (transport sandwich). · **Action:** ask *why* flux forced this off a confident call — is a neighbor/gap-fill mislocated, or is this a real relocation?

**B2. Untrusted-Compartment Placement** · *target: placement*
DeepLoc systematically can't resolve lumen-vs-membrane for most organelles (only the mito m/mm split is trusted per `DEEPLOC_COMPARTMENT_TRUST`). A placement into a low-trust compartment rests on unreliable evidence *however confident the raw probability looks* — a predictor-capability limit distinct from per-call confidence.
- **Fires when:** assigned compartment has low trust (g, ce, lp, v, non-mito membrane splits, and especially **n**).
- **Compute:** `priority = 1 − trust[assigned_comp]`; aggregate to reaction as min trust over its comps; combine multiplicatively with `raw_confidence`.
- **Cost:** cheap · **FP risk:** low-moderate (trust table is yeast/DeepLoc-specific — expose an override). · **Action:** seek orthogonal evidence (DeepTMHMM for membrane, curated DB) for low-trust comps.

**B3. Diffuse / Borderline Placement** · *target: placement*
Rolls the gene-level borderline-margin and diffuse-distribution uncertainty up to the specific reaction placed on it, flagging when the near-tied runner-up is a *different* compartment (one coin-flip from flipping).
- **Fires when:** consensus gene has top-1−top-2 margin < 0.15 **or** normalized score entropy > 0.7.
- **Compute:** per gene, margin and entropy over positive scores; reaction priority = `max` over genes of `0.5·(1−margin) + 0.5·entropy`; flag when runner-up ≠ assigned. Report top-2 comps + margin.
- **Cost:** cheap · **FP risk:** low as an uncertainty measure; many borderline calls are still right — sharpen with neighbor context (Lens D). · **Action:** check whether the runner-up compartment matches the pathway neighbors.

**B4. No-Evidence, Placed-by-Function** · *target: placement*
Catalyzing genes are entirely absent from DeepLoc — compartment chosen purely by the connectivity MILP (`proposal.unplaced_reactions`). The least evidentially-supported placements in the whole proposal.
- **Fires when:** reaction in `unplaced_reactions`, or has no genes but was still moved.
- **Compute:** read `unplaced_reactions` directly. Base priority 1.0; refine up if it created added transports, is essential by single-reaction-deletion, or its EC/subsystem neighbors sit elsewhere. Down-weight spontaneous/transport/exchange reactions (expected to lack GPR).
- **Cost:** cheap · **FP risk:** low as a flag (absence is factual); refinement separates risky from benign. · **Action:** find external localization evidence (orthologs, UniProt, COMPARTMENTS) for the unscored gene.

### Lens C — The seed motif: transport bridging (cheap→medium)

**C1. Transport Sandwich** · *target: both* — **the user's seed signal**
A reaction whose substrate is imported by one added transport and whose product is exported by another is doing round-trip shuttling that vanishes if the reaction itself moved to the neighboring compartment. Two added transports flanking one reaction is the classic fingerprint of a single mislocated reaction (or pathway fragment).
- **Fires when:** a placed reaction R at compartment X has ≥1 **reactant** base with an `added_transport (m_in, X)` **and** ≥1 **product** base with an `added_transport (m_out, X)` — after excluding currency metabolites.
- **Compute:** build `S = {(base, comp)}` from `added_transports`; for each R at X map metabolites to `(base, X)`; fire if a reactant base and a product base both ∈ S for X. Score = fraction of R's non-cofactor metabolites transported across X's boundary, higher when both sides transported and when R's own placement is weak (× B1/B3). Optionally confirm the transports carry flux (structural-only → lower). Report the incoming/outgoing metabolites and the candidate alternative compartment.
- **Cost:** cheap (medium if flux-confirmed) · **FP risk:** moderate — genuine compartmentalized reactions import a substrate and export a product; **currency metabolites (ATP, NAD(H), H⁺, CoA) must be excluded** or it fires on everything. · **Action:** try re-placing R (and immediate neighbors) into the transports' source compartment; if the network still certifies with both transports gone, the placement was wrong.

**C2. Removable-by-Replacement** · *target: transport*
The transport-side view of C1: count how many added transports would *disappear* if one adjacent reaction were re-placed. Targets the highest-leverage transport-error class.
- **Fires when:** a movable reaction is bracketed by an added import and an added export; ranks each transport by transports-removed × placement-weakness of the bracketed reaction.
- **Compute:** bipartite incidence of `added_transports` ↔ consuming/producing reactions; for each bracketed reaction count transports removed by moving it to `default_compartment`. Rank by `(transports removed) × (1 − gene confidence / override flag)`.
- **Cost:** cheap · **FP risk:** low-moderate. · **Action:** re-place the bracketed reaction; confirm growth holds with fewer transports.

**C3. Same-Pathway Split by Flux-Carrying Transport** · *target: both*
A flux-carrying added transport whose two endpoints host reactions of the **same subsystem/EC-pathway** in **different** compartments signals a pathway split the method bridged rather than co-located. Strongest when a single reaction is flanked by two such transports.
- **Fires when:** added transport `(base, c)` carries pFBA/FVA flux **and** the reactions producing/consuming `base` on each side share a subsystem/EC tag but sit in different compartments.
- **Compute:** for each flux-carrying added transport, collect movable reactions touching `base` per side; group by subsystem/EC; fire when both sides share a pathway tag. Rank by shares-pathway ∧ carries-flux ∧ flank-count ≥ 2. Down-weight currency cargo.
- **Cost:** medium · **FP risk:** moderate (genuine inter-compartment shuttles; noisy/missing pathway tags). · **Action:** re-place the minority arm to eliminate the transport.

### Lens D — Topological outliers (cheap→medium, no FBA)

**D1. Neighbor / Pathway Compartment Outlier** · *target: placement*
Merges pairwise-neighbor and subsystem-singleton logic: a reaction should sit where its metabolic neighbors and pathway co-members sit. A lone stray is an outlier.
- **Fires when:** (a) >70% of R's non-currency metabolite-neighbor reactions sit in a single compartment Y ≠ X, **or** (b) R is the sole occupant of X among ≥3 reactions sharing its subsystem tag, the rest concentrating in Y.
- **Compute:** strip currency metabolites; tally neighbor compartments and subsystem-group compartment distribution. Score = `neighbor_fraction_in_Y × log(1 + neighbor_count)` (ignore R with <2 neighbors). Report X, Y, fraction.
- **Cost:** cheap-medium · **FP risk:** currency metabolites create false neighbors (must strip); coarse/missing subsystem tags; genuine interface reactions legitimately disagree with one side. · **Action:** if R is the only pathway member not in Y, re-place it into Y.

**D2. Sole Compartment Occupant** · *target: placement*
A placement making R the only inhabitant of a compartment is suspicious — real organelles host many reactions; a near-empty compartment usually means a score-only placement no metabolic context supports, and it necessarily incurs transports.
- **Fires when:** R's compartment X holds ≤2 placed reactions, especially when X also required added transports to connect R.
- **Compute:** count placed reactions per compartment; score `1/occupant_count`, boosted if R crosses X's boundary via added transports (from C1). Respect known-small compartments (peroxisome, lipid particle) via a per-compartment threshold.
- **Cost:** cheap · **FP risk:** legitimately sparse organelles false-fire; small draft GEMs are sparse everywhere. · **Action:** confirm the compartment is expected to be near-empty; else move R to where its metabolites are.

**D3. Metabolite Pool Fan-Out** · *target: transport*
A base metabolite forced into many compartments (high transport fan-out) is a symptom that reactions touching it were scattered inconsistently.
- **Fires when:** `fanout[m] = |{c : (m,c) ∈ added_transports}|` is in the top ~5% for non-currency metabolites.
- **Compute:** count distinct compartments per base metabolite in `added_transports`; cross-ref the compartment spread of reactions touching m. Score = fanout, normalized to down-weight currency.
- **Cost:** cheap · **FP risk:** currency/ubiquitous metabolites dominate — **must exclude**. · **Action:** find the stray reaction pinning m to an odd compartment; re-placing it often removes one transport.

**D4. Connected-Component Bridge** · *target: transport*
An added transport whose removal disconnects a small reaction island from the network is a load-bearing bridge holding up a possibly-mislocated island driven by scores, not metabolic coherence.
- **Fires when:** an added transport is a cut edge isolating a component of ≤K reactions.
- **Compute:** build the compartmentalized bipartite metabolite-reaction graph; test each added transport for bridge/articulation status. Score inversely with island size; boost if the transport carries no FBA flux. Report the isolated reactions.
- **Cost:** medium (articulation analysis at genome scale) · **FP risk:** genuine peripheral pathways attach via one real transporter. · **Action:** dissolve the island by re-placing it into the default compartment; drop the bridge.

**D5. Nuclear Metabolic Reaction** · *target: placement*
Small-molecule/central/energy/lipid metabolism almost never runs in the nucleus; DeepLoc over-calls nuclear (documented: certified over-connects nucleus 159 vs curated 58; trust for **n** ≈ 0.18). A high-prior error class.
- **Fires when:** reaction placed in nucleus **and** its metabolites are not nucleotide/NA precursors **and** it is not transport/exchange/spontaneous.
- **Compute:** filter placements to nucleus id; exclude nucleotide/DNA/RNA subsystems and transport/exchange/spontaneous. Rank by `trust[n] × #genes`; escalate if the gene has a plausible cytosolic score.
- **Cost:** cheap · **FP risk:** low-moderate (a few real nuclear steps — NAD synthesis, nucleotide salvage — handled by the exclusion). · **Action:** almost always move to cytosol; verify the gene's second-best compartment.

### Lens E — Biochemical plausibility of transports (cheap→medium)

**E1. Impermeant Cargo** · *target: transport*
Merges the cofactor-chemistry and biologically-absent-pool signals. Certain species essentially never cross membranes as free molecules: CoA / acyl-CoA thioesters, NAD(P)(H), FAD/FMN, bulky charged phosphorylated intermediates, tRNA-/enzyme-bound species. An added transport of such cargo means an adjacent reaction should regenerate the cofactor **in-compartment** — the transport is a mathematical patch where reality uses a shuttle or a separate pool.
- **Fires when:** transported base metabolite matches a curated impermeant class (CoA moiety, dinucleotide cofactor, FAD/FMN, |charge| ≥ 2 with MW above a threshold, "tRNA"/"-bound" in name).
- **Compute:** curated blocklist + chemistry rules from formula/charge/name; match each `added_transport`. Score = class-implausibility weight; escalate flux-carrying (load-bearing implausibility) over structural (dead weight). Cross-ref E2 before hard-flagging the few genuine translocases (ADP/ATP, pyruvate, carnitine shuttle).
- **Cost:** cheap · **FP risk:** low-moderate (real shuttles exist — say "is there a shuttle?", not "wrong"). · **Action:** confirm the cofactor is regenerated in the destination; if so, one adjacent reaction is mislocated — trace it.

**E2. No Carrier Evidence** · *target: transport*
The method adds transports on purely topological grounds with zero check that a carrier exists. A crossing for a substrate with no annotated transporter family is the commonest false membrane crossing in draft GEMs.
- **Fires when:** the transported base metabolite has no model gene annotated with a transporter signature (TCDB, transport Pfam, EC 7.x, subsystem "Transport") matching its substrate class.
- **Compute:** resolve metabolite chemical class; scan gene annotations for a matching transporter; flag if none. Rank by `1 − fraction of organism transporters that plausibly match`.
- **Cost:** medium · **FP risk:** moderate (incomplete transporter annotation — use as evidence-of-absence weighting, not a hard reject; pair with existing native transports of the same substrate). · **Action:** check TCDB/Pfam; if no carrier, ask whether the pathway should be re-placed to avoid the crossing.

**E3. Wrong-Destination / Orphan-Import Transport** · *target: transport*
Merges wrong-destination and no-source-production. A transport importing a metabolite into a compartment where nothing produces it — and which is a biologically odd home for that cargo — is a dead-end patch: either a missing biosynthetic reaction or, more often, a mislocated consumer.
- **Fires when:** in the destination compartment the transport is the **only** producer of the base metabolite (no local biosynthesis) **and** the destination is an outlier vs the metabolite's compartment-occurrence consensus.
- **Compute:** per `added_transport (base, c)`, list local producers of `base` in c other than the transport; if none, flag. Build a per-metabolite compartment prior (reference DB or consensus of reactions using it); flag rare destinations. Rank by `#consumers depending on it × absence of local production`.
- **Cost:** medium · **FP risk:** moderate (import is the legitimate way a compartment gets what it can't synthesize — strongest combined with E1/E2). · **Action:** decide missing-local-synthesis vs misplaced-consumer; the transport is a symptom.

**E4. Redundant Parallel Transport** · *target: transport*
Two added transports (or an added one duplicating a native one) moving the same base metabolite across the same compartment pair are mutually redundant — a free deletion candidate that also flags pool over-fragmentation.
- **Fires when:** >1 transport (added or native) moves the same base metabolite between the same unordered compartment pair.
- **Compute:** group **all** transport reactions by `(base, {comp_a, comp_b})`; any group >1 flags each added member. Confirm via FVA that each carries the pool's flux alone. Score by group size.
- **Cost:** cheap · **FP risk:** very low for exact duplicates (FVA removes opposite-direction cycle pairs). · **Action:** delete the redundant added transport; confirm growth unchanged.

### Lens F — Functional / FBA-based signals (medium→expensive)

**F1. Blocked-Everywhere Placement** · *target: placement*
A placed reaction carrying zero flux in **every** certification medium is likely misplaced — its substrates/products don't co-localize with it. Had it been placed where its metabolites are, it would carry flux. (Subtract reactions already blocked pre-assignment to isolate placement-induced blocks.)
- **Fires when:** reaction ∈ `find_blocked_reactions` on the intersection of all media, and was **not** blocked in the pre-assignment draft.
- **Compute:** materialize with `apply_assignment`; per medium `_apply_medium` then `find_blocked_reactions`; intersect. A blocked reaction whose base metabolites are *not* served by any added transport (truly stranded) is more suspect.
- **Cost:** medium · **FP risk:** moderate (pre-existing draft gaps; conditionally-used reactions off-panel — the pre-assignment subtraction handles the former). · **Action:** check whether the reaction's substrates exist in its assigned compartment; if only elsewhere, re-place following substrate localization.

**F2. Dead Multi-Loc Duplicate** · *target: placement*
Merges inactive-duplicate, mono/multi-flip, and futile-cycle-only-flux. A multi-localized copy that carries no productive flux (loopless-FVA span < eps in every medium, or active only through a thermodynamically-infeasible loop) is a dead duplicate harvesting a compartment's score for nothing — the exact unsoundness the enrichment step tries to prevent.
- **Fires when:** a reaction has ≥2 placements and at least one duplicate's loopless-FVA span < `multi_localize_eps` across all media, **or** it carries flux in standard but not loopless FVA.
- **Compute:** for `len(placements[rid])>1`, run loopless FVA (`cycleFreeFlux`) over the duplicate ids `{rid}_{c}` at a growth-holding fraction, per medium — re-running `_enrich_multilocalization`'s own acceptance gate as a final audit. Prioritize duplicates whose added-compartment score barely clears `multi_localize_threshold (0.7)` and whose justifying FVA flux is near `eps`.
- **Cost:** medium-expensive · **FP risk:** low-moderate by design (some duals are genuinely needed only on an off-panel medium; widening the panel reduces this). · **Action:** decide real dual-localization vs scoring artifact; drop the copy or find the medium that justifies it.

**F3. Medium-Conditional Functionality** · *target: placement*
A placement functional on only a subset of the carbon-source panel points to a pathway correctly placed for one metabolic mode but stranded for another (works on glucose, dead on ethanol/oleate because its compartment breaks respiration/β-oxidation).
- **Fires when:** the fraction of media in which the reaction is unblocked/carries FVA flux is low, **and** this media-restriction is new after assignment (placement-induced), not intrinsic to the draft.
- **Compute:** per-medium unblocked/FVA-active vector; score `1 − media_active/media_total`; baseline against the draft's own per-medium activity. Weight reactions whose subsystem matches a condition-specific pathway (respiration, β-oxidation, glyoxylate).
- **Cost:** medium · **FP risk:** low-moderate (many reactions are legitimately condition-specific — the placement-induced delta is the discriminator). · **Action:** check whether the compartment strands the pathway on the failing medium.

**F4. Essential × Uncertain (High-Stakes)** · *target: placement* — **the key multiplicative gate**
A placement whose deletion kills growth is load-bearing for the whole certified result; if it is *also* on shaky localization evidence, one curation error there invalidates the model's functionality claim. High stakes × uncertain = highest curation value.
- **Fires when:** `single_reaction_deletion` drops biomass below the floor (essential) **and** the placement is evidence-weak (low `raw_confidence`, borderline margin, low-trust compartment, or a functionality override).
- **Compute:** `single_reaction_deletion` over movable placements — to save time, restrict to reactions already flagged by the cheap score-only signals (A–E). Intersect essentiality with any Lens B / B1 uncertainty flag. Priority = essential ∧ any-uncertainty.
- **Cost:** expensive · **FP risk:** low for the ranking purpose (essentiality is exact; the only "false positive" is an essential reaction whose evidence is actually strong — cheap insurance to glance at). · **Action:** confirm this enzyme's compartment against literature/UniProt — a wrong call here is the model's single biggest functional risk.

**F5. Fragile Certificate Load-Bearer** · *target: both*
Merges certification-sensitive singleton, sole-lifeline transport, tiny-essential-flux, and floor-fragility. A placement or added transport whose individual removal flips the proposal from certified to uncertified is the sole thing propping up connectivity there. Compounded with the documented growth-floor caveat, these are where a single wrong compartment silently converts "certified" into non-growing. **Tiny-essential-flux** (essential yet `|flux|/biomass < ~1e-3`) is a specific leaked-flux fingerprint.
- **Fires when:** setting this one reaction/transport's bounds to 0 drops biomass below the floor on some medium **and** no parallel route exists; escalate when biomass slack over the floor is small, or when an essential transport carries only trace flux.
- **Compute:** for each added transport and each essential-candidate placement, zero its bounds and re-run `_certify` across media; flag those making `ok=False`. Compute `slack = growths[medium] − min_growth`; low slack = fragile certificate. For essential transports, `|flux|/biomass_flux` from pFBA — trace flux flags a leaked-flux artifact. Cross with B1/gapfill flags: load-bearing **and** override/gap-filled = top of queue.
- **Cost:** expensive · **FP risk:** low (tests are exact; the nuance is degeneracy — some singletons are essential because the draft has a single route, so pair with an evidence-weakness flag). · **Action:** manually verify these single points of failure first; a placement error here breaks the whole certified model.

**F6. Impossible Cofactor Compartment** · *target: placement*
A reaction's EC class / cofactor implies a hard requirement (O₂ for oxidases/monooxygenases EC 1.13/1.14; heme/ubiquinone for the respiratory chain). Placed in a compartment that cannot supply it, the placement contradicts prior knowledge even at high DeepLoc confidence.
- **Fires when:** the reaction's EC or metabolite set marks a compartment-restricted cofactor requirement, and that co-substrate is neither produced in nor transported into the assigned compartment.
- **Compute:** small curated EC-prefix → cofactor table; check the required base-key exists in the assigned compartment (produced or via added transport); if absent, flag. Cross-check by blocking the requirement in that compartment and confirming the reaction goes blocked in FVA.
- **Cost:** medium · **FP risk:** moderate (coarse EC→cofactor mapping; ubiquitous cofactors — fire only when the cofactor is genuinely compartment-restricted). · **Action:** if the compartment can't supply the cofactor, move the reaction (and likely its neighbors) to the canonical compartment for its EC class.

### Lens G — External priors (cheap when data present; expensive for orthology)

**G1. Contradicts Curated Localization DB** · *target: placement*
When a curated source (UniProt/COMPARTMENTS/HPA) was one of the score sources, a placement landing where the curated source explicitly does **not** list the gene — while it *does* list another compartment — is a direct contradiction of authoritative prior knowledge, stronger than any DeepLoc-only signal.
- **Fires when:** the assigned compartment is absent from the gene's non-empty curated authoritative set.
- **Compute:** require pre-fusion per-source `LocalizationScores`; take the curated source's non-zero compartments per gene; flag assigned ∉ set when set non-empty. Score by curated-source confidence; reduce severity if it was a functionality override (the method knew it fought evidence). Guard with a gene-id-overlap sanity check.
- **Cost:** cheap · **FP risk:** low when the DB is trustworthy (moderate if id-mapping is imperfect). · **Action:** treat the DB as a ground-truth candidate — either the placement is wrong or it's a real dual-localization the DB captured and the model missed.

**G2. Ortholog Compartment Outlier** · *target: both*
Merges the three cross-species signals. With orthology across organisms, a reaction/gene whose ortholog is consistently in compartment Y in sibling species but landed in X here is a conservation-breaking outlier — a strong prior single-organism DeepLoc can't provide.
- **Fires when:** an orthogroup is available and this organism's placement is the minority against a strong consensus (≥70–75% of orthologs agree on Y, mean `raw_confidence` high).
- **Compute:** group placements by orthogroup; consensus = confidence-weighted mode; fire on `assigned ≠ consensus` with strong support. Priority = `consensus_support × consensus_confidence × (1 − this organism's own-placement evidence)`. Silently skip when single-organism.
- **Cost:** expensive (requires multi-organism run) · **FP risk:** moderate (real lineage-specific relocalization; orthogroup-mapping errors — treat as prior, not verdict; down-weight orthogroups that themselves disagree). · **Action:** check whether this species truly relocated the enzyme vs an assignment error; else align to consensus.

**G3. Override × Prior Violation (Compounded)** · *target: placement* — **most conservative, highest-precision**
The method's own "I fought the evidence" flag (functionality override) **combined with** an external-prior violation (impossible cofactor compartment **or** curated-DB contradiction) isolates placements most likely to be a topology artifact, not real biology — usually the pathway should have been relocated wholesale, or a gap-fill/transport is missing.
- **Fires when:** a placement is a functionality override **and** it also fires F6 or G1.
- **Compute:** intersect the override set with F6/G1 hits. Priority = override severity × external-prior-violation strength.
- **Cost:** cheap (given the component signals) · **FP risk:** low by construction (requires two independent flags to co-fire). · **Action:** find the upstream missing reaction/transport that forced the reaction into a biologically wrong compartment; fixing it lets the reaction return home.

### Lens H — Repair-forced & gap-fill provenance (cheap→medium)

**H1. Growth-Gap Pin** · *target: placement*
`_diagnose_growth_gap` pins the sole movable producer of a blocked biomass precursor to the biomass compartment (protected from confinement). Placed by a last-resort single-producer heuristic, not evidence — if that gene is really elsewhere, the network is missing a transporter/isozyme and the pin masks a real gap.
- **Fires when:** a reaction is the unique movable producer of a biomass precursor and was forced to the biomass compartment (reconstruct `gap_pinned`), especially placed against its gene's top evidence.
- **Compute:** find biomass precursors blocked in the pre-repair placement; the movable reaction touching that precursor base forced to `bio_comp` is a pin. Rank by proximity of `growths[__primary__]` to `min_growth` (load-bearing + fragile).
- **Cost:** medium · **FP risk:** low-moderate. · **Action:** verify the pinned reaction is the correct producer, or whether an isozyme/transporter should supply the precursor.

**H2. Confinement-Forced Group** · *target: placement*
`_diagnose_confinement` drags movable reactions into a pinned reaction's compartment when they share a **non-transportable** metabolite. If that metabolite *should* have been transportable (a curation error in the transportable set), the whole forced group is mislocated as a downstream consequence — review the group as a unit.
- **Fires when:** a reaction's compartment was set by the confinement repair (in the `forced` map or a co-location group), especially where its own gene's top compartment differs.
- **Compute:** re-run `_diagnose_confinement` on final placements; group by shared non-transportable base; flag members whose compartment = the forced target. Also flag the shared confined metabolite. Rank by `#reactions dragged × evidence disagreement of the dragged members`.
- **Cost:** medium · **FP risk:** moderate (co-locating a genuinely confined pathway is usually correct — the signal is *should this metabolite be transportable?*). · **Action:** check the shared metabolite's transportability; if it should have a transporter, free the whole group.

**H3. Arbitrary Gap-Fill** · *target: placement*
`proposal.added_reactions` are pulled from a universal DB purely to restore biomass — **no** localization evidence, **no** organism gene support, chosen by a parsimonious objective with many equal-cost optima (the specific pick is arbitrary). Inherently the most uncertain additions in the proposal.
- **Fires when:** the reaction id ∈ `added_reactions`.
- **Compute:** enumerate directly. Escalate by (a) essentiality (`single_reaction_deletion` below floor = most suspect), (b) non-uniqueness (re-run gapfill with `iterations>1`), (c) whether the reaction's chemistry/compartment is even consistent with an organism lacking the gene. Rank essential, non-unique gap-fills first.
- **Cost:** medium · **FP risk:** low (gap-fills are uncertain by definition). · **Action:** confirm the organism plausibly has this reaction; find the real gene or replace the arbitrary universal pick.

---

## 3. Recommended first implementation

Ship a composite from **seven signals** chosen for the best cost/value ratio: all cheap to compute over the full proposal except two targeted expensive checks run only on already-flagged candidates. This covers the four highest-value error modes — evidence overridden, the seed transport motif, external-prior contradiction, and high-stakes essentiality — without a genome-scale FBA sweep on turn one.

| # | Signal | Target | Cost | Why first |
|---|--------|--------|------|-----------|
| 1 | **B1 Override Against Confident Evidence** | placement | cheap | The method's headline risk; sharp confident evidence sacrificed for flux. |
| 2 | **C1 Transport Sandwich** | both | cheap | The user's seed hypothesis; direct mislocation fingerprint. |
| 3 | **B4 No-Evidence Placed-by-Function** | placement | cheap | Zero-evidence placements read straight from `unplaced_reactions`. |
| 4 | **E1 Impermeant Cargo** | transport | cheap | Highest-precision transport-error class; pure chemistry lookup. |
| 5 | **D1 Neighbor/Pathway Outlier** | placement | cheap | Catches mislocations B1/C1 miss; independent topological view. |
| 6 | **H3 Arbitrary Gap-Fill** | placement | cheap | Enumerated directly; inherently uncertain additions. |
| 7 | **F4 Essential × Uncertain** | placement | expensive (gated) | Run `single_reaction_deletion` only on candidates 1–6, then multiply. |

**How to combine.** Give each target a raw score per signal in [0,1], then:

```
base   = 1 − Π_i (1 − w_i · s_i)          # noisy-OR: any strong signal lifts the score,
                                           # multiple weak ones accumulate, saturates at 1
stakes = 1 + β · essential_flag            # F4 multiplier, β ≈ 0.6
priority = clip(base · stakes, 0, 1)
```

Suggested weights `w_i` (reflecting precision, higher = more trusted): B1 0.9, C1 0.85, E1 0.85, H3 0.8, B4 0.7, D1 0.6. Noisy-OR is deliberate over a weighted sum: it rewards **corroboration** (a reaction that is both an override *and* transport-sandwiched should rank far above either alone) and won't let one noisy currency-metabolite hit dominate. The `stakes` multiplier implements the "high-stakes × uncertain" principle — essentiality alone never flags a reaction, but it *amplifies* an already-uncertain one.

**Two hard preprocessing rules** (or the whole ranking is dominated by noise):
1. **Currency-metabolite exclusion list** (ATP/ADP, NAD(P)(H), H⁺, H₂O, CoA, Pi, PPi, FAD) applied to C1, D1, and all transport-fan-out logic. Without it, the seed signal fires on nearly every reaction.
2. **Subtract the pre-assignment draft baseline** for any FBA-derived signal (F1, F3, F4) so pre-existing gaps aren't blamed on placement.

**How to present the output.** A single **ranked table**, most-suspect first, one row per target:

| priority | target_id | type | placement (X→?) | top flags | suggested action |
|---|---|---|---|---|---|
| 0.94 | `R_ACOAH_m` | placement | m (gene top: c, 0.97) | Override-vs-confident · Transport-sandwich (m_in acACP, m_out 3hbcoa) · Essential | Try re-placing to cytosol; check if both flanking transports vanish and growth holds |

Columns: composite priority, target id, placement with the gene's overridden top compartment in parentheses, the top 2–3 contributing flags as chips, and the single concrete curator action from the highest-weighted firing signal. Below the table, a **per-flag expander** for each row listing every firing signal with its raw score and one-line reason (which isozyme dissented, which metabolites are transported, the essentiality delta), plus the candidate alternative compartment so the curator can act without re-deriving anything. Default the view to the **top 3–5%** by priority, with the full list one click away — that slice is where the override, seed-motif, and essential-uncertain catches concentrate, and it is the few-percent the user asked to surface first.

## 4. v1 status — shipped, with one design correction from genome-scale testing

`curation_priority()` ships in [`localization/curation.py`](../../src/raven_toolbox/localization/curation.py) (tests in `tests/test_localization_curation.py`). It implements the section-3 composite — override, sandwich, no-evidence, neighbour-outlier, gap-fill, impermeant, and the gated essentiality multiplier — combined by noisy-OR with the currency exclusion. An adversarial multi-lens review + genome-scale run on yeast-GEM (2296 placements) forced several corrections over the naive first cut; the substantive ones:

- **The topology signals had to be gated by evidence weakness — the key finding.** The seed hypothesis was that "substrate imported, product exported" (sandwich) and "neighbours in another compartment" (neighbour-outlier) are *rare, suspicious* motifs. At genome scale they are the opposite: this method relocates reactions into organelles **by bridging their metabolites with transports**, so the sandwich motif is the *normal* output — it fired on 47 % of reactions, neighbour-outlier on 43 %, flagging ~70 % of the model and burying the real errors. The evidence distribution is bimodal: 86 % of placements have a confident gene supporting the assigned compartment, 14 % have none. The fix multiplies both topology signals by **weakness = 1 − DeepLoc-support(assigned compartment)**, so they only speak when the placement is *also* evidentially unsupported — the true signature of a questionable call. This took the flagged set from ~70 % to ~20 % of reactions with a clean, separable top few percent (priority ≥ 0.9: 0.3 % of reactions; the top-8 are the fully-corroborated override+sandwich+neighbour+essential rows), and made the ranking informative rather than saturated. `override` is inherently evidence-based already; `no_evidence`/`gap-fill`/`impermeant` are independent of the gate.
- **Essentiality is measured against the materialised model's own achievable growth, not the assignment's `min_growth` floor.** The floor can exceed what the applied model actually reaches (the documented growth-floor gap); comparing knockout growth to it would flag every candidate and collapse the stakes multiplier uniformly. The gate now takes the applied baseline and skips itself when it cannot discriminate.
- **override uses `support(top) − support(assigned)`, not a top1−top2 margin**, which collapses when two organelle sub-compartments (mito matrix vs membrane) both score high — exactly the confident-override cases the signal exists to surface; and it fires only when the top compartment is *not* one of the reaction's placements, so genuine dual-localisation (FUM1/ACO1/tRNA-synthetases) does not false-fire.
- **Topology signals score per (reaction, compartment) over the full placement list**, so dual-localised second copies are surfaced; **currency/impermeant matching normalises the base key** (`M_` prefix, `[compartment]` tag) so the exclusion holds on id-keyed models; **impermeant now includes free NAD(P)H/FAD/FMN** (shuttle-only cargo), catching e.g. `NADPH→v`, `acetyl-CoA→p` transports on yeast.

Rule 2 of the section-3 preprocessing (subtract the draft FBA baseline) applies to flux-*delta* signals (F1/F3); the shipped essentiality gate is a lethality test against the applied baseline, which is the same intent realised more directly. Open follow-ups: rank the flat impermeant band (all free-CoA/cofactor transports tie at 0.85); optionally expose the per-signal raw scores for the UI expander; validate the ranking against curated yeast-GEM relocation errors once independent ground truth (P2) is in.

**Coupling (`n_affected` / `affected`).** Curation is not independent: re-placing or removing one reaction changes the metabolic context of everything sharing its metabolites at that compartment, and removing a transport starves whatever depends on it. Each row therefore carries its **one-hop blast radius** — the other placed reactions that curating it would touch (for a placement: reactions sharing a non-currency metabolite with it at the same compartment; for a transport: reactions using the transported species at either end of the bridge). `n_affected` is the full count (sort by it to separate high-impact curations that deserve care or a batched review from isolated ones safe to fix alone — on genome-scale yeast this ranges 0 to ~112, median 3, with ~160 rows coupling to ≥10 reactions), and `affected` lists up to a dozen coupled ids, the ones also flagged in the table first. The coupling is deliberately one-hop, not the transitive connected component, so the number stays interpretable; each coupled reaction has its own radius the curator can follow.
