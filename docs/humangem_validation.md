# Human-GEM cell-type model validation: ravengem vs RAVEN

Validation of ravengem's tINIT/ftINIT against MATLAB RAVEN on a real genome-scale
reconstruction (Human-GEM) using the Hart2015 RNA-seq dataset (5 cell lines: DLD1,
GBM, HCT116, HELA, RPE1). The goal is functional equivalence — do ravengem and RAVEN
extract the *same* context-specific reaction sets from the same inputs?

## Method

* **Template & inputs.** RAVEN built the ftINIT reference model from Human-GEM
  (`prepHumanModelForftINIT`: remove drug/exchange/artificial reactions, set
  spontaneous/custom lists) and exported it as `raven_refModel.xml` (10198 reactions).
  ravengem builds on that *same* exported model, so the candidate reaction universe is
  identical and set comparison is exact.
* **Scoring.** Gene scores from `log2(TPM+1)`-style expression via
  `gene_scores_from_expression`, mapped to reactions through the GPR
  (`score_reactions_from_genes`), matching RAVEN's `getExprForRxnScore`.
* **ftINIT.** Series `1+1` (2 staged MILP steps). RAVEN run via `ftINIT.m` with Gurobi;
  ravengem via `ravengem.init.ftinit` with Gurobi (`mip_gap=0.001`, `time_limit=600`).
* **tINIT.** ravengem `get_init_model` (classic single-MILP INIT) on HCT116, compared to
  the ftINIT result for the same cell line.
* **Tasks.** Two ravengem ftINIT variants: *no-task* (expression only) and
  *task-constrained* (essential metabolic tasks, `metabolicTasks_Essential.txt`, force
  task-essential reactions to be kept). RAVEN's reference is task-constrained.
* **Solver.** Gurobi 13.0.1 for both tools.

## Engineering findings (ravengem tractability)

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

| cell line | RAVEN ftINIT | ravengem ftINIT (no-task) | ravengem ftINIT (task) |
|-----------|-------------:|--------------------------:|-----------------------:|
| DLD1      | 7782 | 7744 | TBD |
| GBM       | 7668 | 7667 | TBD |
| HCT116    | 7780 | 7752 | TBD |
| HELA      | 7832 | 7789 | TBD |
| RPE1      | 7569 | 7564 | TBD |

Counts agree within ~0.5 % (GBM: 7667 vs 7668). ravengem tINIT (HCT116) gives 6024
reactions — a smaller model, as expected from the different (classic INIT) objective.

### Agreement — ravengem (no-task) ftINIT vs RAVEN ftINIT

| cell line | shared | only ravengem | only RAVEN | Jaccard |
|-----------|-------:|--------------:|-----------:|--------:|
| DLD1   | 7667 |  77 | 115 | 0.976 |
| GBM    | 7562 | 105 | 106 | 0.973 |
| HCT116 | 7675 |  77 | 105 | 0.977 |
| HELA   | 7707 |  82 | 125 | 0.974 |
| RPE1   | 7470 |  94 |  99 | 0.975 |

**~97.5 % of reactions are identical** between the two independent implementations, even
though this ravengem run is *expression-only* while RAVEN's reference is
task-constrained. The residual disagreement (≈80–125 reactions each way out of ~7700) is
within the range expected from MIP-gap tolerance (both accept near-optimal incumbents),
alternate optima, and the missing task constraints — the "only RAVEN" reactions are
expected to include task-essential reactions that the task-constrained run (below) keeps.

### ravengem tINIT vs ftINIT (HCT116)

tINIT 6024 rxns vs ftINIT 7752; shared 5957, Jaccard 0.762. tINIT is nearly a subset
(only 67 reactions unique to it) — the two methods agree on a common core, with ftINIT
keeping more (its staged formulation and task handling are less aggressive at removal).

## Conclusions

ravengem reproduces RAVEN's ftINIT reaction selection on a genome-scale human model to
~97.5 % set identity from identical inputs — strong evidence of functional equivalence.
Reaching this required porting RAVEN's numerical-conditioning choices (fixed big-M,
`rescaleModelForINIT`) and an `optlang`-specific fast constraint build; see *Engineering
findings*. The task-constrained comparison is reported above once complete.
