# Localization (Phase 7) — critical review + redesign

Design note for porting RAVEN's `predictLocalization.m` to raven-toolbox. **Status: proposal,
not yet implemented**; the user-facing API and the algorithmic choices are settled here
before code lands, because RAVEN's all-or-nothing redistribution doesn't match how the
function actually gets used.

## 1. What RAVEN's `predictLocalization` does, and where it falls short

The MATLAB function takes a single-compartment (or to-be-merged) model + a per-gene per-
compartment score table (`GSS`) and assigns every reaction to a compartment by simulated
annealing, maximising `Σ gene_score(g,c) − transport_cost`. Decisions baked into that
shape that we should not transcribe blindly:

1. **Destroys existing compartmentalization.** "If the model contains several
   compartments they will be merged." Existing curated assignments — typically the most
   reliable signal in a draft GEM — are thrown out and recreated from the predictor
   alone. This is fine for a single-compartment Phase-3 draft but obviously wrong as a
   default for a refinement step on a multi-compartment model.
2. **Closed-model requirement.** "The input model should also not include any exchange,
   demand or sink reactions, otherwise this function would not provide any results."
   The caller must pre-process; the function has no graceful path for open models.
3. **Silent connectivity surgery.** "Reactions that are unconnected are removed and
   saved in removedRxns" — and the iterative "metabolite must be produced somewhere"
   loop quietly deletes any reaction whose substrate has no producer. **This is the
   "model must be complete" assumption**: a draft with dangling substrates loses
   reactions without the caller necessarily noticing.
4. **Reactions without GPRs get fake genes** (`&&FAKE&&N`) scored 0.5 everywhere — they
   carry no signal but still have to be placed.
5. **Mono-localization for genes** ("a gene can only be assigned to one compartment …
   This is a simplification to keep the problem size down"). Biologically wrong for the
   many dual-localized enzymes (mitochondrial/cytosolic isoforms etc.); a constraint of
   the solver, not the biology.
6. **Simulated annealing** with a `maxTime` cap — stochastic, no optimality certificate,
   slow on large models.
7. **`expandModel` first**, splitting isozyme reactions before optimisation. Changes the
   reaction set rather than handling isozymes in scoring.
8. **No partial-update mode.** Cannot say *"I added these 50 reactions, only place
   those, keep the existing assignments"*.
9. **No suggest-only mode.** Mutates the model in place; no preview.
10. **Plotting and the predictor are baked into the function** (`plotResults`,
    GSS-shape input).

Items 1–3 and 8–9 are exactly what the user flagged. Items 5–7 are also worth
reconsidering on a modern stack.

> **Final API (post-implementation):** the proposal below was refined during
> implementation. The shipped API:
>
> * `reactions_to_relocate` is a **required** caller-passed set of IDs (no
>   `notes['localization']='uncertain'` auto-detect — one mechanism is clearer).
> * **Multi-compartment is the default scoring model.** No `multi_compartment_genes`
>   boolean. The highest-scoring compartment a gene lands in is "free"; every
>   additional compartment costs `multi_compartment_penalty` *plus* its (typically
>   lower) predictor score is its own implicit penalty. Pick a large penalty for
>   effectively mono-localised genes.
> * `mergeCompartments` and `copyToComps` are ported separately as
>   `raven_toolbox.manipulation.merge_compartments` / `copy_to_compartment` (they're useful
>   independently of `predict_localization` — for flattening for analysis or building
>   dual-localised pathways).
> * `mapCompartments` is **not** ported — its main use case overlaps with
>   `compare_models`.

## 2. Proposed `predict_localization` for raven-toolbox

Decompose the function into independent concerns:

* **A pluggable predictor input** — a `pandas.DataFrame` of *gene × compartment* scores
  (higher = stronger evidence). Loaders for WoLF PSORT, DeepLoc, … convert their output
  to this table; the algorithm is predictor-agnostic. Ships with `load_wolfpsort()` and
  documents the `gene_id × compartment → score` contract so users can plug in DeepLoc /
  TargetP / their own.
* **A deterministic MILP** instead of simulated annealing.
* **An explicit scope of reactions to (re-)place** — by default *all reactions whose
  compartment is unset or marked uncertain*; user can pass `reactions_to_relocate=[…]`
  to lock everything else.
* **Existing compartments respected by default** — they are the high-confidence prior,
  not noise to be merged away.
* **Open models tolerated.** Boundary reactions (exchange/sink/demand) are pinned to
  their current compartment automatically and excluded from relocation.
* **Incomplete models tolerated.** Dangling metabolites are allowed; no silent reaction
  removal. The MILP places what it can and reports the rest in `unplaced_reactions`,
  not by deleting them.
* **`apply=False`** returns a `LocalizationProposal` (a diff: which reactions move,
  which transports are added) without mutating the input.

### 2.1 API sketch

```python
@dataclass
class LocalizationScores:
    """Per-gene compartment scores. Index = gene_id; columns = compartment ids.
    Missing genes / NaN scores fall back to a uniform prior."""
    df: pd.DataFrame

@dataclass
class LocalizationProposal:
    """What the predictor proposes, before applying it."""
    moved: pd.DataFrame              # rxn_id, from_compartment, to_compartment, score_delta
    added_transports: pd.DataFrame   # met_id, from, to, cost
    unplaced_reactions: list[str]    # couldn't be placed (e.g. no scored genes)
    objective: float

@dataclass
class LocalizationResult:
    """Result of applying a proposal (or running with apply=True directly)."""
    model: cobra.Model
    gene_compartments: dict[str, str | set[str]]    # set when multi_compartment_genes=True
    added_transports: list[cobra.Reaction]
    proposal: LocalizationProposal

def predict_localization(
    model: cobra.Model,
    scores: LocalizationScores,
    *,
    # === The two requested features ===
    reactions_to_relocate: Iterable[str] | None = None,   # None ⇒ pick automatically (see below)
    keep_existing: bool = True,                            # don't merge compartments first
    apply: bool = True,                                    # False ⇒ return proposal only
    # === Other knobs ===
    default_compartment: str = "c",
    transport_cost: float | Mapping[str, float] = 0.5,
    multi_compartment_genes: bool = False,                 # relax mono-localization
    require_producibility: bool = False,                   # off ⇒ don't drop unproducible mets
    mip_gap: float = 0.001,
    time_limit: float | None = None,
) -> LocalizationResult | LocalizationProposal:
    """Assign reactions to subcellular compartments via a deterministic MILP.

    By default (``reactions_to_relocate=None``, ``keep_existing=True``), only reactions
    that lack a compartment assignment (or are marked uncertain via a notes flag) are
    placed; reactions already in a compartment stay put. Pass an explicit
    ``reactions_to_relocate`` list to override.

    With ``apply=False`` the function returns a :class:`LocalizationProposal` describing
    what would change — useful for review or for diffing against a curator's choices.
    """
```

### 2.2 What gets placed by default

| reaction state | default behaviour |
|---|---|
| has a (single) compartment assignment in the existing model | **pinned** (kept as-is) |
| boundary reaction (exchange / sink / demand) | **pinned** to its current compartment |
| no compartment set OR marked uncertain (e.g. `notes['localization'] == 'uncertain'`) | **relocated** |
| caller passed `reactions_to_relocate=[…]` | **only those** relocated; rest pinned |

This handles both user requests in one design:

* *Incomplete model*: dangling mets / missing producers are tolerated; reactions are
  placed if any scored gene supports them, and unplaceable ones land in
  `unplaced_reactions`.
* *Selective re-localization*: pass the subset; defaults respect what's there.

### 2.3 MILP formulation (mono-localization variant, `multi_compartment_genes=False`)

Let `R*` = reactions to relocate, `G*` = unique genes appearing in those reactions, `C`
= compartments, `M` = metabolites occurring in `R*` reactions.

Variables (all binary):
* `x[r, c]` for `r ∈ R*`, `c ∈ C` — reaction `r` placed in compartment `c`.
* `y[g, c]` for `g ∈ G*`, `c ∈ C` — gene `g` assigned to compartment `c`.
* `t[m, c]` for `m ∈ M`, `c ∈ C ∖ {default}` — transport of `m` between `default` and
  `c`.

Constraints:
* Each relocated reaction goes to exactly one compartment: `Σ_c x[r,c] = 1`.
* Each gene goes to exactly one compartment: `Σ_c y[g,c] = 1`.
* Gene-reaction coupling: if reaction `r` requires gene `g` (per its GPR), then placing
  `r` in `c` requires `g` in `c`: `x[r,c] ≤ y[g,c]` for the AND-clause; OR-clauses use
  the standard linearisation (any isozyme satisfying the placement is enough).
* Metabolite presence in each compartment: if any reaction touching `m` is placed in
  `c`, then `m` must exist in `c` — modelled by adding a per-compartment "demand for
  transport" when needed. Concretely: if a metabolite participates in reactions placed
  in `c ≠ default`, a transport `t[m,c]` is required unless another reaction in `c`
  balances it.

Objective:
```
maximise  Σ_{g,c} gene_score(g,c) · y[g,c]   −   Σ_{m,c} transport_cost(m) · t[m,c]
```

With `multi_compartment_genes=True`, the `Σ_c y[g,c] = 1` constraint is relaxed (gene
can be in several compartments), and gene-reaction coupling becomes
`x[r,c] ≤ y[g,c]` per compartment as before — but the gene now contributes its
compartment-specific score once per assigned compartment. (The cost of allowing this is
that the same gene's score is double-counted across compartments; an alternative is to
include a small penalty per extra compartment, controllable via a `multi_compartment_penalty`
argument — leave this to a later iteration.)

### 2.4 Why MILP rather than simulated annealing

RAVEN's SA was chosen in an era when MATLAB MILP solvers were less accessible. Today
we have Gurobi already wired in. The MILP gives:
* deterministic, reproducible answers (important for science),
* an explicit optimality certificate / gap,
* faster solves on the small/medium problems this is (≪ 30k binaries even for a
  whole-genome model with 10 compartments),
* the same well-understood `mip_gap` / `time_limit` controls the rest of raven-toolbox uses,
* graceful degradation: at the time limit, return the best incumbent rather than the
  "current SA state" with no quality guarantee.

The SA's main practical benefit was *good-enough answers without a license-encumbered
MILP solver*. We've already documented that genome-scale (f)tINIT needs Gurobi
([docs/init_solver_benchmark.md](init_solver_benchmark.md)); reusing the same backend
is consistent. (For users on open-source solvers, GLPK on this MILP should be tractable
because it's smaller than ftINIT — the per-reaction big-M is at most `|C|`, not 1000.)

## 3. Pluggable predictors

The algorithm consumes a `gene × compartment` score table. raven-toolbox ships:
* `load_wolfpsort(path)` — RAVEN-compatible WoLF PSORT output → score table.
* (later) a thin adapter for **DeepLoc** (TSV with per-class probabilities).
* The format is open: a user can build the DataFrame from any source.

A small `omics.scoring` helper could compute a uniform "no-evidence" row for genes
absent from the predictor's output (RAVEN's 0.5-everywhere fallback).

## 4. What this delivers for the user's two requests

* **Incomplete model OK** — `require_producibility=False` (default) means unbalanced
  metabolites are tolerated; reactions are placed when their genes have signal, and
  unplaceable ones are *reported*, not silently deleted. No "model must be complete"
  precondition.
* **Selective re-localization** — `reactions_to_relocate=[…]` pins everything else.
  With `apply=False` the user gets a diff to review; with `apply=True` only the
  selected reactions move. The default (no `reactions_to_relocate` argument) places
  only the uncompartmentalised / uncertain reactions, which is the natural "update the
  draft" workflow.

## 5. Open questions to resolve before implementing

1. **How does raven-toolbox mark "uncertain" compartmentalization?** Proposal: a
   `notes['localization'] = 'uncertain'` flag on reactions; or a passed-in set of ids.
   The "auto-pick-what-to-place" mode reads this.
2. **Multi-compartment gene scoring**: simple multi-counting (the same score in every
   compartment the gene lands in) vs a penalty per extra compartment. Start with simple
   multi-counting; add the penalty later if users ask.
3. **DeepLoc adapter** — ship now or later? DeepLoc has a stable TSV format; a 30-line
   loader is trivial when there's a real use case. Defer until requested.

## 6. Implementation plan

1. `localization/scores.py` — `LocalizationScores`, `load_wolfpsort`.
2. `localization/milp.py` — `predict_localization` (MILP).
3. `localization/apply.py` — convert a proposal into a re-compartmentalised cobra model
   (added transports, updated metabolite compartments).
4. Tests on a small handcrafted model + a regression test against the RAVEN paper's
   yeast example.
5. PROGRESS / PLAN updates.
