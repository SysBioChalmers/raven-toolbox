# Homology cut-off calibration — protocol

:::{admonition} Protocol, not results
:class: warning

This is the design. Nothing has been measured yet; no default should be changed
on the strength of this page. It exists so the method is agreed — and the
success criterion fixed — *before* the numbers arrive, which is the only way the
chosen cut-off is a finding rather than a curve fitted after the fact.
:::

Homology-based reconstruction decides which template reactions transfer to a new
organism on the strength of sequence hits. The thresholds that make that decision
have never been calibrated in either toolbox: RAVEN and raven-toolbox share the
same three values, and neither records where they came from.

## Which parameter is actually under study

Not the aligner's cut-off. `run_blast` / `getBlast` filter at collection time —
1e-4 in both, after
[raven-toolbox#91](https://github.com/SysBioChalmers/raven-toolbox/pull/91) —
but `get_model_from_homology` / `getModelFromHomology` then filter the hit table
again:

```python
keep = (hits.evalue <= max_evalue) & (hits.align_len >= min_align_len) & (hits.identity >= min_identity)
```

with defaults `1e-30`, `200`, `40` on both sides. Since 1e-30 is twenty-six
orders stricter than the collection threshold, **every hit that survives model
building would have survived collection anyway**: the aligner cut-off is inert
for the default pipeline, and tuning it would change nothing.

The levers are therefore `max_evalue`, `min_align_len` and `min_identity` — plus
`strictness` / `bidirectional`, which decide whether a hit must be reciprocal.

## Ground truth

Two sources, chosen because their weaknesses do not overlap.

### Primary: KEGG, across organisms

Build the template organism's model from KEGG, use it to reconstruct a target
organism by homology, and compare against that target's *own* KEGG-derived
model. Both models carry KEGG reaction ids, so the comparison needs no
identifier mapping — the reason this design is worth preferring over comparing
against a curated GEM in a foreign namespace. It also scales: any organism in
KEGG can be a target, which is what makes a distance series affordable.

:::{admonition} This measures agreement with KEGG, not correctness
:class: caution

KEGG's KO assignments come from KOFAM HMMs and SSDB, and SSDB is itself built
from all-against-all BLAST. The "ground truth" is therefore *not* independent of
sequence similarity, and a cut-off tuned to maximise agreement with it is tuned
to reproduce KEGG's own homology calls. That is informative — KEGG's calls are
curated and widely used — but it is not an accuracy estimate, and the results
document must say so in the same breath as any number it reports.
:::

A second limitation: a KEGG organism model is annotation-driven, so it is a
*lower bound* on real metabolic content. Reactions the homology draft adds beyond
it are not automatically wrong. Metrics keep missing and extra separate rather
than folding both into one score.

### Secondary: curated GEM pairs

The non-circular check. Use **yeast-GEM** as the template and reconstruct
organisms for which a curated GEM already exists in the same namespace —
**hanpo-GEM** (*Hansenula polymorpha*) and **rhto-GEM** (*Rhodotorula
toruloides*). Curation is independent of BLAST, and the shared namespace again
avoids id mapping. Few data points, but they answer the question the KEGG sweep
cannot: does the tuned optimum survive contact with a model somebody checked by
hand?

If the two ground truths disagree about the optimum, that disagreement is the
most interesting result available and belongs in the write-up rather than being
averaged away.

## Measure at the hit level first

Going straight to models buffers the signal. A reaction transfers if *any* of its
template genes has a surviving hit, so reaction-set agreement is insensitive to
the thresholds across a wide range, and the sweep would look flat for reasons
that have nothing to do with the parameter being right.

KEGG gives gene→KO for both organisms, so each hit can be judged directly:

> for a hit between template gene *g*ₛ and target gene *g*ₜ, do they share a KO?

That yields precision, recall and F1 of **ortholog calls** as a function of the
thresholds — the quantity the parameters actually control — and it is cheap to
compute. Model-level metrics then show how far that propagates.

## Organisms

A distance series, because the optimum is a function of divergence and a single
target would produce a number that generalises to nothing. Template throughout:
*S. cerevisiae* (`sce`).

| Band | Proposed targets | Why |
|---|---|---|
| Close | `kla`, `zro` | same family; nearly everything should transfer |
| Medium | `yli`, `cal` | the working range for real reconstructions |
| Distant | `ncr`, `ani` | filamentous fungi; where thresholds start to bite |
| Very distant | `ath`, `eco` | the regime where homology transfer should mostly *fail*, and a cut-off that still transfers freely is wrong |

Confirm the KEGG organism codes against the organism table before running
(`organisms_in_domain`, `stream_organism_gene_ko`); the codes above are proposed
from memory and one wrong code invalidates a row.

The very-distant band matters more than it looks. Every threshold study risks
optimising recall until it admits everything; including pairs where the right
answer is "few reactions transfer" keeps the objective honest.

## The sweep is cheap

Run the aligner **once per organism pair** and cache the hit table. Every
threshold combination is then pure post-processing of that table — no realignment
— so the cost of the study is one BLAST run per pair (minutes), not one per
grid point.

Proposed grid, one parameter at a time first, then a local grid around whatever
that suggests:

| Parameter | Values | Default |
|---|---|---|
| `max_evalue` | 1e-100, 1e-50, 1e-30, 1e-20, 1e-10, 1e-5 | 1e-30 |
| `min_align_len` | 50, 100, 150, 200, 300 | 200 |
| `min_identity` | 20, 30, 40, 50, 60 | 40 |
| `bidirectional` | true, false | true |

Repeat the whole sweep on a **DIAMOND** hit table for the same pairs. If the two
aligners imply different optima, the thresholds are compensating for aligner
behaviour rather than describing biology, and the defaults should differ per
aligner — which neither toolbox currently allows.

## Metrics

Per organism pair and threshold combination:

**Hit level** (against KO sharing) — precision, recall, F1 of ortholog calls.

**Reaction level** (against the target's KEGG model) — recall (fraction of the
KEGG model's reactions recovered), extra count (reactions transferred that KEGG
does not list, reported separately, *not* as false positives), Jaccard for a
single summary figure.

**GPR level** — of the reactions present in both, the fraction whose gene sets
agree. This is where a threshold change shows up first: a reaction survives on
one gene while its complex quietly loses members.

Report per band, not just pooled. A parameter that helps `kla` and hurts `ncr`
is a different finding from one that helps everything.

## Success criterion, fixed in advance

The objective is stated here so it cannot be chosen after seeing the curves:

> **Maximise hit-level F1 against KO sharing, averaged over the medium and
> distant bands, subject to the very-distant pairs transferring no more
> reactions than the current defaults do.**

Close pairs are excluded from the objective because almost any threshold works
there; very-distant pairs act as a constraint rather than a target. A change to
the defaults is proposed only if it improves that objective by a margin larger
than the spread across organisms within a band — otherwise the honest conclusion
is "the current values are as good as anything measured", which is a perfectly
good result and one this study should be willing to reach.

## Confounds to control

- **E-values depend on search space.** BLAST's e-value scales with query length
  times database size, so a cut-off tuned on one proteome does not transfer to a
  larger or smaller one. Report proteome sizes alongside every result, and check
  whether `min_identity` and bitscore — which are size-independent — give a
  flatter optimum across the series. A plausible outcome of this study is that
  `max_evalue` is the wrong knob and identity is the right one.
- **Template completeness.** A reaction missing from the *template* can never
  transfer, and would be scored as a miss. Report the template model's coverage
  of each target's KEGG model as a ceiling on achievable recall.
- **Bidirectionality interacts with the thresholds** — a strict cut-off applied
  to both directions is stricter than it looks. Sweep it explicitly rather than
  holding it fixed.
- **KEGG version.** Pin one release for the whole study and record it; KO
  assignments move between releases.

## Deliverables

1. A results document in this directory, following
   [kegg_hmm_cutoff_calibration.md](kegg_hmm_cutoff_calibration.md): method,
   measurements, chosen defaults with rationale, and a cross-validation section.
2. A reproducible driver under `scripts/`, taking the KEGG release and organism
   list as arguments.
3. If the defaults move, they move **away from RAVEN's** — the same situation as
   the KO-assignment cut-offs, where measurement justified diverging (1e-30
   against RAVEN's 1e-50, gene ratio 0.9 against 0.8). That divergence needs
   recording in the raven-gecko-parity ledger, which currently has no status for
   "deliberately different, with evidence", and a back-port proposal so MATLAB
   can follow if the evidence convinces.
4. A parity scenario pinning whichever values are chosen, so the two toolboxes
   cannot drift apart on them again.
