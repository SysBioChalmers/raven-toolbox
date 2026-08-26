# Model manipulation parameter benchmarks

Functions: `raven_toolbox.manipulation`

Date: 2026-06-20.

---

## `remove_genes` — `blocked_reactions` policy

**Parameters tested:** `'keep'` (MATLAB default), `'remove'` (Python default)

Test model: e_coli_core (95 reactions, 72 metabolites, 137 genes).
Target gene: `b1779` (encodes GAPD, glyceraldehyde-3-phosphate dehydrogenase).
GAPD is known to be essential for aerobic growth on glucose in E. coli.

| `blocked_reactions` | Reactions after removal | Predicted growth | Correct? |
|---|---|---|---|
| `'keep'` (MATLAB) | 95 (unchanged) | 0.874 mmol gDW⁻¹ h⁻¹ | No — GAPD kept with empty gene rule |
| `'remove'` (Python) | 94 | 0.000 | Yes — GAPD absent, glycolysis broken |

With `'keep'`: the reaction's stoichiometry is preserved after the gene rule is cleared.
The FBA objective (growth) treats it as a freely available reaction with no gene
constraint, giving full growth. This is the correct behaviour when using `remove_genes`
for annotation/curation tasks (preserving reaction capacity while editing gene rules),
but it is wrong for gene-essentiality prediction.

With `'remove'`: single-gene reactions (GAPD, catalysed by b1779 only) are deleted.
The broken glycolytic pathway correctly produces zero growth.

**Decision: keep `blocked_reactions='remove'`.** Python's default is correct for the
primary use case (gene essentiality, metabolic engineering). Users who need MATLAB's
`'keep'` behaviour (annotation curation) should pass `blocked_reactions='keep'`
explicitly — add a migration note to the docstring.

---

## `allow_excretion` in `get_init_model`

**Context:** `get_init_model` used to default to `allow_excretion=True` while the
higher-level wrappers `run_init` and `run_ftinit` default to `allow_excretion=False`.
This was an inconsistency in default values, not a difference in algorithm.

**Parameters tested:** `True` (former `get_init_model` default), `False` (now the default)

Test model: synthetic dead-end model with a blocked metabolite `r_BD`.
Test condition: `prod_weight=0.5` (the default) and `prod_weight=0` (edge case).

| `allow_excretion` | `prod_weight` | Objective | `r_BD` kept? | Correct? |
|---|---|---|---|---|
| `True` | 0.5 | 28 | Yes | same |
| `False` | 0.5 | 28 | Yes | same |
| `True` | 0 | 28 | Yes | — |
| `False` | 0 | 16 | No | — |

With `prod_weight=0.5` (the default), `allow_excretion` has **no effect**. The INIT
formulation creates sink variables `s_{met}` for each metabolite with a weight
`prod_weight × score`; when `prod_weight > 0` these sinks absorb net metabolite
production regardless of the `allow_excretion` flag. The flag only activates an
additional exchange reaction pathway that is redundant when sinks are present.

With `prod_weight=0` the sinks are removed, so `allow_excretion=True` is the only
way to allow net metabolite production — and gives a different (larger) objective.

**Decision: `get_init_model` default changed to `allow_excretion=False`** (done
2026-06-20, `6f3b57c`) to match `run_init` and `run_ftinit`. This has no effect on
results in the default workflow (`prod_weight=0.5`) but removes the confusing
inconsistency. The docstring note that `allow_excretion` only has an effect when
`prod_weight=0` was added in the same commit.

---

## `constrain_reversible_reactions` — `eps`

**Parameter:** `eps=1e-9` (Python default, no MATLAB equivalent)

`eps` is the absolute tolerance used to detect whether a reversible reaction's
forward or backward optimum is functionally zero (below `eps`) before applying
the bound constraint. An `eps` that is too small would leave near-zero optima
unconstrained; too large would wrongly constrain reactions with real small fluxes.

**Status: untested** — no sensitivity benchmark run. The value `1e-9` matches
Gurobi's default primal feasibility tolerance and is consistent with cobrapy
conventions. No known failures on yeast-GEM or iJO1366 in current tests.

Proposed future test: run `constrain_reversible_reactions` on iJO1366 at eps=1e-12,
1e-9, 1e-6 and compare the number of reactions constrained and the resulting FBA
objective.
