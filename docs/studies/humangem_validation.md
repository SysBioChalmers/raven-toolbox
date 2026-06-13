# Human-GEM cell-type model validation: raven-toolbox vs RAVEN

Validation of raven-toolbox's tINIT/ftINIT against MATLAB RAVEN on a real genome-scale
reconstruction (Human-GEM) using the Hart2015 RNA-seq dataset (5 cell lines: DLD1,
GBM, HCT116, HELA, RPE1). The goal is functional equivalence — do raven-toolbox and RAVEN
extract the *same* context-specific reaction sets from the same inputs?

## Method

* **Template & inputs.** RAVEN built the ftINIT reference model from Human-GEM
  (`prepHumanModelForftINIT`: remove drug/exchange/artificial reactions, set
  spontaneous/custom lists) and exported it as `raven_refModel.xml` (10198 reactions).
  raven-toolbox builds on that *same* exported model, so the candidate reaction universe is
  identical and set comparison is exact.
* **Scoring.** Gene scores from `log2(TPM+1)`-style expression via
  `gene_scores_from_expression`, mapped to reactions through the GPR
  (`score_reactions_from_genes`), matching RAVEN's `getExprForRxnScore`.
* **ftINIT.** Series `1+1` (2 staged MILP steps). RAVEN run via `ftINIT.m` with Gurobi;
  raven-toolbox via `raven_toolbox.init.ftinit` with Gurobi (`mip_gap=0.001`, `time_limit=600`).
* **tINIT.** raven-toolbox `get_init_model` (classic single-MILP INIT) on HCT116, compared to
  the ftINIT result for the same cell line.
* **Tasks.** Two raven-toolbox ftINIT variants: *no-task* (expression only) and
  *task-constrained* (essential metabolic tasks, `metabolicTasks_Essential.txt`, force
  task-essential reactions to be kept). RAVEN's reference is task-constrained.
* **Solver.** Gurobi 13.0.1 for both tools.

## Engineering findings (raven-toolbox tractability)

Getting ftINIT to run at genome scale surfaced three issues, all now fixed and matching
RAVEN's design:

1. **O(n²) constraint construction.** Building the steady-state balances with Python
   `sum()` re-canonicalises a growing sympy expression at each term; hub metabolites
   (ATP/H⁺/H₂O in ~10³ reactions) made one constraint take ~minutes (≈154 s total build,
   benchmark: 1500-term `sum` = 59 s vs `optlang.symbolics.add` = 0.01 s). Fixed by
   building flat term lists once per reaction and summing with `optlang.symbolics.add`
   (in both ftINIT and tINIT).
2. **Big-M too loose.** The on/off indicator constraints used each reaction's own bound
   (±1000) as big-M; with `force_on=0.1` that is a ~10⁴ ratio → very weak LP relaxation
   → Gurobi never closes the gap. RAVEN uses a fixed big-M = 100. Adopted.
3. **Stoichiometric rescaling.** A fixed big-M=100 is only valid if no reaction needs
   flux ≫100; ported RAVEN's `rescaleModelForINIT` (cap each reaction's coefficient
   dynamic range at 25×, normalise mean |coeff| to 1) into `prep_init_model`. Without it
   the staged MILP is infeasible (step-1 caps transports that step-0 used freely).

Net effect: a full ftINIT cell-line solve went from *not finishing* to ~200 s,
comparable to RAVEN.

## Results

### Reaction counts

| cell line | RAVEN ftINIT | raven-toolbox ftINIT (no-task) | raven-toolbox ftINIT (task) |
|-----------|-------------:|--------------------------:|-----------------------:|
| DLD1      | 7782 | 7744 | 7774 |
| GBM       | 7668 | 7667 | 7680 |
| HCT116    | 7780 | 7752 | 7776 |
| HELA      | 7832 | 7789 | 7816 |
| RPE1      | 7569 | 7564 | 7570 |

Counts agree within ~0.5 % everywhere; the task-constrained run is closest (e.g. RPE1
7570 vs 7569, HCT116 7776 vs 7780). raven-toolbox tINIT (HCT116) gives 6024 reactions — a
smaller model, as expected from the different (classic INIT) objective.

### Agreement — raven-toolbox (no-task) ftINIT vs RAVEN ftINIT

| cell line | shared | only raven-toolbox | only RAVEN | Jaccard |
|-----------|-------:|--------------:|-----------:|--------:|
| DLD1   | 7667 |  77 | 115 | 0.976 |
| GBM    | 7562 | 105 | 106 | 0.973 |
| HCT116 | 7675 |  77 | 105 | 0.977 |
| HELA   | 7707 |  82 | 125 | 0.974 |
| RPE1   | 7470 |  94 |  99 | 0.975 |

**~97.5 % of reactions are identical** between the two independent implementations, even
though this run is *expression-only* while RAVEN's reference is task-constrained. The
"only RAVEN" surplus (≈99–125) is expected to include task-essential reactions that the
task-constrained run (below) recovers.

### Agreement — raven-toolbox (task-constrained) ftINIT vs RAVEN ftINIT

| cell line | shared | only raven-toolbox | only RAVEN | Jaccard |
|-----------|-------:|--------------:|-----------:|--------:|
| DLD1   | 7699 | 75 | 83 | 0.980 |
| GBM    | 7588 | 92 | 80 | 0.978 |
| HCT116 | 7696 | 80 | 84 | 0.979 |
| HELA   | 7735 | 81 | 97 | 0.978 |
| RPE1   | 7493 | 77 | 76 | 0.980 |

Adding the essential metabolic tasks (same task list RAVEN uses) raises agreement to
**Jaccard 0.978–0.980** and makes the disagreement symmetric (only-raven-toolbox ≈ only-RAVEN
≈ 80), confirming the prediction: the task constraints recover RAVEN's task-essential
reactions. The residual ≈80 reactions each way out of ~7700 is at the level expected from
MIP-gap tolerance (both accept near-optimal incumbents) and alternate optima.

### raven-toolbox tINIT vs ftINIT (HCT116)

tINIT 6024 rxns vs ftINIT 7752; shared 5957, Jaccard 0.762. tINIT is nearly a subset
(only 67 reactions unique to it) — the two methods agree on a common core, with ftINIT
keeping more (its staged formulation and task handling are less aggressive at removal).
This matches the expected tINIT/ftINIT relationship rather than indicating a defect.

## Conclusions

From identical inputs on a genome-scale human reconstruction, raven-toolbox reproduces RAVEN's
ftINIT reaction selection to **97.5 % (no-task) and 98 % (task-constrained) set identity**
across five cell lines — strong evidence of functional equivalence between the two
independent implementations. Agreement is symmetric and the residual (~80 reactions each
way) is consistent with MIP near-optimality and alternate optima rather than any
systematic divergence.

Reaching genome-scale tractability required matching RAVEN's numerical-conditioning
choices and fixing optlang-specific construction costs (see *Engineering findings*):
fixed big-M = 100, `rescaleModelForINIT`, `optlang.symbolics.add` instead of Python
`sum()` in every MILP build (ftINIT, tINIT, and the gap-fill). With these, a
task-constrained cell-line model builds in ~15–25 min (dominated by the
essential-forced staged MILP) and a no-task one in ~3 min, comparable to RAVEN.
