# Cross-language parity checks

These tests answer one question: **does raven-toolbox still agree with MATLAB
RAVEN?** Until now that agreement was reported in study documents and nothing
would fail if it stopped being true.

Run them with:

```bash
pytest -m parity
```

They are marked so the default suite stays fast and self-contained; nothing here
is required for raven-toolbox to be correct on its own terms.

## What "agree" means

"Identical output" is achievable for some functions and meaningless for others,
so the checks are written in three tiers. Putting a function in the wrong tier
produces either a flaky test or a vacuous one.

### Tier 1 — exact

Deterministic transformations of a model or a file: I/O round-trips, task-list
parsing, gene-association normalisation, elemental balance, identifier sorting,
model merging, reversibility splitting, KEGG table parsing. Values must match.
A disagreement means one implementation is wrong.

Note that *exact* means semantically exact, not byte-identical. RAVEN and
raven-toolbox both write valid YAML but differ in key order and quoting style,
which carries no meaning; comparing bytes would test the serialiser's habits
rather than the model.

### Tier 2 — set-level

Mixed-integer results: INIT/ftINIT extraction, gap-filling, compartment
assignment. These problems have many optima of equal objective value, so a
different answer is not a wrong answer. Compare overlap (Jaccard, containment)
against a recorded baseline band, and hold the solver fixed — Gurobi and GLPK
can legitimately land on different optima of the same problem.

### Tier 3 — statistical

Flux and random sampling. Two runs of the *same* implementation differ. Compare
distributions at a fixed seed, never individual samples.

## Where the reference values come from

Two sources, with different trade-offs:

1. **RAVEN's own artefacts** (`test_yaml_interop.py`). RAVEN's repository
   contains models written by MATLAB RAVEN itself. Reading them with
   raven-toolbox is a genuine cross-language check that needs no MATLAB
   installation — but it does need a RAVEN checkout, so these tests skip unless
   `RAVEN_ROOT` points at one:

   ```bash
   git clone --depth 1 -b develop3 https://github.com/SysBioChalmers/RAVEN
   RAVEN_ROOT=$PWD/RAVEN pytest -m parity
   ```

   RAVEN is GPL and raven-toolbox is MIT, so those files are **read where they
   are** and never copied into this repository.

2. **Recorded oracles** (`test_oracles.py`). For behaviour that cannot be
   inferred from a file — what `checkTasks` reports, how `getElementalBalance`
   grades a reaction — MATLAB has to be run once and its answers recorded.
   `scripts/parity/generate_oracles.m` does that, reading the fixtures in
   `tests/data/parity/` (authored here, so no GPL material is involved) and
   writing JSON into `tests/data/parity/oracles/`. Those tests skip when the
   oracle file is absent, so a checkout without them still passes.

   To regenerate, in MATLAB with RAVEN on the path:

   ```matlab
   cd scripts/parity
   generate_oracles
   ```

   Commit the resulting JSON, and note the RAVEN commit it came from in
   `oracles/PROVENANCE.md`.

## Determinism

`test_determinism.py` is not cross-language: it checks that raven-toolbox gives
the *same* answer twice. It lives here because it protects the same property the
parity tiers do — several recent fixes made compartment placement and gap-filling
deterministic, and nothing would have caught a regression.
