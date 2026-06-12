# ftINIT (Phase 4d) — critical review and implementation plan

Critical read of RAVEN's ftINIT before porting (memory: *review very critically, don't
transcribe*). Covers `ftINIT.m`, `ftINITInternalAlg.m`, `getINITSteps.m`,
`INITStepDesc.m`, `prepINITModel.m`, `mergeLinear.m`, `groupRxnScores.m`,
`rescaleModelForINIT.m`, `ftINITFillGaps*.m` (~3000 lines).

## 1. What ftINIT is, vs the tINIT we already have

tINIT (Phase 4c, `init/`) solves **one** MILP: pick the reaction subset that best
matches expression scores while staying flux-consistent. We ported its core as
`run_init` (+ `get_init_model`). ftINIT ("fast tINIT") keeps the *objective* but
attacks the **runtime** of the MILP on genome-scale models (Human-GEM: ~12k rxns)
with four ideas:

1. **One-time preprocessing (`prepINITModel` → `prepData`).** All omics-independent
   work — simplification, essential-reaction discovery via tasks, reaction
   classification, linear merging, scaling — is done once per *template* model and
   reused across every sample. (RAVEN comment: prep is ~1 h; each model is then fast.)
2. **Linear merging (`mergeLinear`).** Reactions connected through a degree-2
   metabolite (one producer, one consumer) are contracted into one, summing scores.
   Human-GEM 11888 → 7922 rxns — a ~33 % smaller MILP, losslessly.
3. **Staged MILP (`getINITSteps`, default `'1+1'`).** Instead of one big MILP, run
   two smaller ones: step 1 decides only GPR-associated reactions (transport/
   exchange/spontaneous/no-GPR *ignored* = left in, not in the problem); step 2 fixes
   step-1 reactions as **essential** and decides the remaining no-GPR reactions. The
   `'full'` series is a single MILP ≈ classic tINIT.
4. **A cheaper MILP formulation (`ftINITInternalAlg`).** Positive-score reactions use
   a **continuous** on-indicator instead of a binary (see §3) — roughly halving the
   integer count, the dominant cost.

After the MILP: map merged groups back to original reactions, delete the rest,
**gap-fill so all metabolic tasks are feasible** (`ftINITFillGapsForAllTasks`),
re-add exchanges, optionally `removeLowScoreGenes`.

## 2. Architecture (data flow)

```
template model ──prepINITModel──> prepData {refModel, refModelWithBM, minModel(merged+scaled),
  (once)                                     groupIds, origRxnIds, essentialRxns,
                                             essentialMetsForTasks, taskStruct, toIgnore* masks}
                                                │
omics (tissue) ──scoreComplexModel──> rxnScores │  (per sample)
                                                ▼
                              ftINIT step loop (getINITSteps)
                              ├ step i: groupRxnScores(mask_i) ─> ftINITInternalAlg (MILP) ─> on/off
                              │         prev results folded in as 'essential'
                              └ ...
                                                ▼
                         delete off-reactions (mapped through groupIds)
                                                ▼
                         ftINITFillGapsForAllTasks (per-task gap-fill MILP)
                                                ▼
                         re-add exchanges, removeLowScoreGenes ─> model
```

## 3. The core MILP (`ftINITInternalAlg`) — formulation

Reactions are partitioned into **six categories**, each with its own variables.
`forceOnLim = 0.1` (min flux to count a reaction "on"); big-M `= 100` (rev split) and
`1000`/own-ub elsewhere; `intTol = 1e-7`.

| category | indicator | key constraint | binary? |
|---|---|---|---|
| **PosIrrev** (score>0, irrev) | `Y∈[0,1]` continuous | `v ≥ forceOnLim·Y` | **no** |
| **PosRev** (score>0, rev) | `Y∈[0,1]` + dir-binary | split `v=v⁺−v⁻`; `v⁺+v⁻ ≥ 0.1·Y`; binary picks direction | 1 (dir only) |
| **NegIrrev** (score<0, irrev) | `Y∈{0,1}` | `v ≤ 100·Y` | **yes** |
| **NegRev** (score<0, rev) | `Y∈{0,1}` | split; `v⁺+v⁻ ≤ 0.1·Y` | yes |
| **EssIrrev** (forced on) | — | `lb = min(0.99·|prev flux|, 0.1)` | no |
| **EssRev** (forced on, rare) | dir-binary | split + force ≥0.1 | 1 |

**The key insight:** a *positive*-score reaction needs no binary — a continuous `Y`
suffices, because the objective *maximises* `Σ score·Y` so it pushes `Y→1` whenever
flux can flow; the `v ≥ 0.1·Y` gate makes `Y>0` impossible without flux. Only
*negative*-score reactions need a true binary (the objective would otherwise leave
their `Y` at 0 for free). Reactions with **score 0 are left out of the problem
entirely** (always present, can carry flux) — this is how the `toIgnore*` categories
are excluded (their score is zeroed by `groupRxnScores`).

Objective: minimise `Σ(−score)·Y − prodWeight·Σ mon` (metabolite-production bonus).
`prodWeight = 5` (passed by `ftINIT.m`, overriding the docstring default 0.5).

We already have a *simpler* version of this in `run_init`: binaries for all,
`eps·x ≤ v ≤ ub·x`, reversible split, `prod_weight` sinks. ftINIT's formulation is
the optimised superset.

## 4. Critical findings (what NOT to transcribe)

**Magic numbers** — collected, to be named constants and (where they bite) calibrated
like we did for the HMM cut-offs (K15):

| value | where | role | concern |
|---|---|---|---|
| `forceOnLim = 0.1` | InternalAlg | min "on" flux | scale-dependent; RAVEN itself muses "test 0.01" |
| big-M `100` / `1000` | InternalAlg / FillGapsMILP | on/off gates | must exceed real fluxes; ties to ±1000 clamp |
| `intTol 1e-7` / `1e-9` | InternalAlg / FillGapsMILP | integer tol | solver-specific (Gurobi v10 numerics) |
| `prodWeight = 5` | ftINIT.m | met-production bonus | "has not been evaluated" (their words) |
| score clamps `±0.1`, `0.01` | ftINIT.m / groupRxnScores | avoid 0 scores | MILP can't have exactly-0 score (binary flips freely) |
| `maxStoichDiff = 25` | rescale | coeff-ratio cap | docstring says 250 — **doc/code disagree** |
| MIPGap schedule `0.0004→0.003`, abs-gap `10/20`, TimeLimit `120/5000` | getINITSteps | per-step solver tuning | Gurobi-tuned; HiGHS will differ |
| int-extraction `>1e-3`, `>0.5` | FillGapsMILP / InternalAlg | read binaries | tolerance-driven |

**RAVEN/Gurobi coupling & quirks**
- `ftINITFillGapsMILP` **hard-errors on glpk-via-COBRA**; "only tested with Gurobi".
  raven-toolbox uses **optlang**, so we are solver-agnostic — but the magic numbers and
  `Seed=26` reproducibility were Gurobi-tuned. Validation must expect *equally optimal
  but not identical* reaction sets (alternative optima); compare on objective value
  and task feasibility, not exact rxn identity.
- Several quirks worth fixing upstream too: `rescaleModelForINIT` doc/code mismatch;
  `ftINITFillGaps` has a dead default-score branch referencing an undefined `models`;
  `prepData` does **not** store `origRxnIds` (downstream relies on implicit ordering —
  we will store it explicitly); the 2-column `b` (ranged mass-balance RHS) is not a
  cobra concept and must be built as optlang ranged constraints.

**Genuine simplification/improvement opportunities** (log in IMPROVEMENTS, back-port):
- `simplifyModel`'s ~1 h reversibility-removal pass is an FVA in disguise; cobra's
  optimised (optionally parallel) FVA does it far faster — the single biggest prep
  speedup available to us.
- Carry `rxnScores` as an aligned array, not a smuggled model field threaded through
  `simplifyModel`/`removeReactions`/`mergeModels`.
- Expose the magic numbers as parameters with documented, calibrated defaults rather
  than burying them; `forceOnLim`/`prodWeight`/big-M are scale-dependent and deserve
  the same empirical treatment as K15.
- `removeLowScoreGenes` (gene pruning by isozyme/complex role) is a clean standalone
  utility, reusable beyond ftINIT.

**Prerequisite gap in our codebase:** `check_tasks` currently returns only
pass/feasible per task. `prepINITModel` needs **essential-reaction discovery** (which
reactions, and in which direction, are essential for each task) and the task
metabolite set. That work was explicitly deferred from 4a/4c to here.

## 5. What cobra/raven-toolbox already gives us

- MILP via **optlang** (HiGHS bundled, Gurobi/CPLEX if licensed) — `run_init`
  already builds INIT-style MILPs this way; ftINIT extends that builder.
- FVA / pFBA / `find_blocked_reactions` (`simplifyModel`, gap-fill feasibility).
- `check_tasks` (relaxed mass-balance via constraint bounds = RAVEN's `b`), to be
  extended with essential-reaction output.
- `connect_blocked_reactions` (connectivity MILP) — related to but not the same as
  the task gap-fill MILP (which adds min-cost template reactions to make a task
  feasible); the task gap-filler is new.
- `score_reactions_from_genes` / `gene_scores_from_expression` (scoring).

## 6. Implementation plan — phased sub-steps

Correctness first on small models (the RAVEN `tinitTests.m` cases define reaction
scores directly via single-gene expression — ideal oracles), then the speed layers.

| sub-phase | deliverable | builds on | risk |
|---|---|---|---|
| **4d.0** | **Test oracles**: port the `tinitTests.m` toy models + `getINITSteps`-style score→expression helper, so every later sub-phase is checked against RAVEN's expected on/off sets. | tasks, score | low |
| **4d.1** | **Essential-reaction discovery**: extend `check_tasks` to return the per-task essential-reaction matrix + direction + task-metabolite set (the missing `prepINITModel` step-3 input). | tasks/check.py | med |
| **4d.2** | **prepData preprocessing** (`prep_init_model`): simplify (cobra FVA for the reversibility pass), essential rxns (4d.1), reaction-classification masks (exchange/import/transport/spont/no-GPR), **linear merge** (`merge_linear` + `group_ids`, storing `orig_rxn_ids`), rescale. | manipulation/simplify, 4d.1 | **high** (merge is fiddly) |
| **4d.3** | **Core staged MILP** (`ftinit_internal` + step loop): the 6-category formulation extending `run_init`; `INITStep` dataclass + `get_init_steps('1+1'/'2+1'/'full'/…)`; prev-results-as-essential folding; per-step MIPGap retry schedule. Start **without** metabolomics. | init/init.py | **high** |
| **4d.4** | **Task gap-filling** (`ftinit_fill_gaps_for_tasks`): per-task, LP-feasibility-gated, min-cost template-reaction MILP (ranged `b` constraints in optlang). | gapfilling, 4d.1 | med |
| **4d.5** | **Assembly + `removeLowScoreGenes`**: map merged groups back, delete, re-add exchanges, prune low-score genes; top-level `ftinit(...)` entry. | score, 4d.2–4 | low |
| **4d.6** | **Metabolomics block** — **deferred** (decided 2026-05-26). The linear merge eliminates degree-2 detected metabolites, so a clean flux-based bonus is impossible; it needs RAVEN's producer-group-mapping + `mon`/`vnrbm`/`vnrvm`/`vnim` negative-producer force-flux blocks — the most intricate MILP in ftINIT, for its least-used input (not in the RNA-seq/single-cell/HPA ranking; "not designed for metabolomics only"), with a randomness-laden oracle (T0008). `ftinit(metabolomics=…)` raises `NotImplementedError`. | 4d.3 | med |
| **4d.7** | **Validation + calibration**: run `'full'` vs `'1+1'` on a real template (Human-GEM or a smaller curated GEM), compare to RAVEN by objective/task-feasibility/size (not exact rxns); calibrate `forceOnLim`/`prodWeight`/big-M; document like K15. | all | med |

Suggested module layout under `init/`: `prep.py` (4d.1–2), `ftinit.py` (4d.3 + top
level), `taskfill.py` (4d.4), `genes.py` (`remove_low_score_genes`), reusing
`init.py`/`score.py`/`build.py`.

## 7. Decisions (locked 2026-05-25)

1. **Solver: agnostic.** Write the MILP to optlang's generic interface; let cobra's
   configured solver decide. Calibrate the default constants on **HiGHS** (open,
   bundled) but do not hardcode anything Gurobi-specific. Validation compares on
   objective/task-feasibility, never exact reaction identity.
2. **Correctness-first.** Implement the single-step `'full'` MILP (extend `run_init`,
   **no** linear merge) and match RAVEN on the toy oracles *first*; then add linear
   merge and `'1+1'`/`'2+1'` staging as separately-verified speed layers. → reorders
   the sub-phases below: **4d.3 (full MILP) precedes 4d.2 (merge)**.
3. **Metabolomics deferred to 4d.6.**
4. **Validation includes Human-GEM** (in 4d.7): toy oracles + small GEM for fast
   iteration, then Human-GEM + a tissue dataset for the real ftINIT regime.

### Resulting execution order

`4d.0 oracles → 4d.1 essential-rxn discovery → 4d.3 full MILP (correctness gate) →
4d.2 linear merge → 4d.3b '1+1'/'2+1' staging → 4d.4 task gap-fill → 4d.5 assembly +
removeLowScoreGenes → 4d.6 metabolomics → 4d.7 calibration + Human-GEM validation.`

The §6 table lists deliverables by id; this is the order they'll actually be built.
