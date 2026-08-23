# Roadmap

Where raven-toolbox and its documentation are going, in the order the work has to happen.
This page is the **plan**: phases, what each one delivers, what it depends on, and how you
know it is finished. [todo.md](todo.md) is the **backlog**: the item-level detail behind
each phase.

Two repositories are in play:

* **raven-toolbox** (this repo) — code, tests, CI, maintainer reference documents.
* **[raven-docs](https://github.com/edkerk/raven-docs)** — the dual-language documentation
  site (MkDocs + Material on Read the Docs) that serves both MATLAB RAVEN and
  raven-toolbox, with both toolboxes as git submodules. All user-facing prose belongs there.

Phases are **dependency-ordered, not calendar-dated**. Phase 0 blocks nothing and should
happen immediately; phases 1 and 2 run in parallel; phase 3 depends on decisions D2 and D4
below. There is no external 1.0 anchor (decision D5) — 1.0 ships when phase 4's exit
criteria are met, so **phases 3 and 4 may overlap freely** rather than running in sequence.

---

## The problem this roadmap solves

The port itself is essentially done — the functional scope of MATLAB RAVEN is covered, with
three deliberate omissions, and the results are validated on Human-GEM, yeast and a
multi-organism set. What is *not* done is everything that makes that work usable and
trustworthy by someone who is not the author:

1. **A user cannot find out how the two toolboxes differ.** The one page that tries to say
   so is factually wrong (it names Python functions that do not exist).
2. **Nothing enforces the parity claims.** They live in prose in study documents; no test
   would fail if the two implementations diverged tomorrow.
3. **There is no Python tutorial anywhere.** The site's guides are 100 % MATLAB while
   claiming to be dual-language.
4. **No one can tell which functions are mature.** A five-times-validated pipeline and a
   thin wrapper look identical from the outside.

Each phase below closes one of those, in the order that makes the next one cheaper.

---

## Phase 0 — Stop publishing things that are false

*No dependencies. Small. Do it first.*

The documentation site currently misleads readers, and every day it stays up is a day
someone copies a function name that does not exist.

* Remove or correct the invented Python names in `differences.md` (`get_blast`,
  `fill_gaps`, `import_model`, `export_to_excel_format`, `solve_lp`, …).
* Drop the claim that RAVEN Excel round-trips between the two toolboxes — Python has no
  Excel import by design.
* Correct `protocol/index.md`, which promises protocols "in both MATLAB and Python with
  MATLAB/Python tabs" when the site has no Python content at all.
* Unbreak the `matlab-vs-python.md` reference — it is both an in-page link and a `mkdocs.yml`
  nav entry pointing at a file that does not exist.

**Exit criterion:** nothing on the published site is factually wrong. Interim honesty
("Python examples coming") is acceptable; incorrect specifics are not.

## Phase 1 — Make the differences documentation true *and* self-maintaining

*Depends on: phase 0. Runs parallel to phase 2.*

Fixing the names by hand solves today and rots by next release. The deliverable is a
**generator**, and a set of small pages it feeds.

* A build-time generator in raven-docs that emits the MATLAB ↔ Python mapping from the two
  submodules' API data and **fails the build on any hand-written function name that does
  not resolve**. This is the durable fix, and the reason phase 1 is worth more than a
  correction pass.
* The `differences/` section, one page per question: `index`, `mapping`, `python-only`,
  `matlab-only`, `behaviour`, `parity`.
* Behind it, in this repo: extend [migration.md](migration.md) with the reverse direction
  (121 of 185 public names have no row today) so the generator has something to read.

**Exit criterion:** a wrong or stale function name cannot reach the published site, and
every public function of either toolbox appears in the mapping with a counterpart or an
explicit "no equivalent, because …".

**Note:** `differences/behaviour.md` and `differences/parity.md` can be *drafted* here from
what is already documented, but they are only *finished* by phase 2 — that is what supplies
the evidence.

## Phase 2 — Make parity a test, not a claim

*Depends on: phase 0 only. Runs parallel to phase 1; supplies phase 1's last two pages.*

Today parity rests on hand-transcribed oracles, one YAML gate, and prose in study
documents. The deliverable is a harness that fails a build when the two implementations
disagree.

* **Tier the contract first** — exact / set-level / statistical — because "identical output"
  is not achievable uniformly, and pretending otherwise produces either flaky tests or
  vacuous ones. This tiering is also the content of `differences/parity.md`.
* Golden-fixture harness: committed small inputs, a MATLAB driver that writes oracles, and
  `pytest -m parity` that runs green with no MATLAB installed.
* Determinism regressions for the recently-fixed placement and gap-fill paths (#76, #83,
  `c239d2e`), which currently have no guard.
* A nightly job on a licensed runner for the Gurobi-dependent tier-2 checks, so the study
  numbers stop being refreshed by hand.

**Exit criterion:** every parity claim on the site is backed by a test that would fail if it
stopped being true, and the behaviour table is populated from measured divergences rather
than from memory.

## Phase 3 — A new tutorial, dual-language from the first line

*Depends on: decisions D1 and D2 below; benefits from phase 2's fixtures.*

The site has one protocol — the published hanpo-GEM reconstruction — and it is MATLAB-only.
**This roadmap does not port it.** It stays as it is: a faithful rendering of a published
protocol, MATLAB, untouched, alongside the untouched legacy tutorials 1–5.

Instead we **define a new tutorial**, written dual-language from the start. It may
legitimately resemble the hanpo-GEM protocol in structure and in several steps — homology
draft, biomass, gap-filling, simulation, curation are the natural arc of any reconstruction
— but it is its own exercise with its own organism and its own runnable scripts in both
languages.

Deliberately **out of scope for now:** the SLIME lipid-curation treatment. Parallel work is
under way on those methods, so the new tutorial does not depend on them and does not attempt
a Python equivalent. This removes what would otherwise have been the one hard blocker.

**Two reconstruction tutorials, one per route** (decision D1). RAVEN has two ways to build a
draft, they suit different situations, and no current guide covers the second at all:

1. **Reconstruction by homology — a non-conventional yeast.** A curated template plus
   BLAST/DIAMOND, structurally parallel to the existing protocol but its own organism and
   inputs. Written first: the arc is well understood, the inputs are small, and no new
   infrastructure is needed. **It stops short of lipid curation** — biomass is defined
   without the SLIME treatment, which keeps the tutorial clear of the parallel method work.
2. **Reconstruction *de novo* from KEGG — a small prokaryote.** The route with no coverage
   anywhere on the site, and a raven-toolbox strength (artefact builders, HMM cutoff
   calibration). Second, because it carries an infrastructure dependency the homology route
   does not: the KEGG artefacts have to be fetched, so the tutorial has to say clearly what
   is downloaded, how large it is, and what can be skipped. The existing artefact/manifest
   tooling and [maintaining_kegg_data.md](../maintenance/maintaining_kegg_data.md) are the
   starting point.

Then, in order:

3. **GEM extraction (ftINIT).** Adapt the existing material, which is already in good shape,
   rather than writing it fresh.
4. **GEM comparison.** Already listed as planned on the site.
5. Then, as capacity allows: compartmentalisation, confidence tracking, ecModel handoff —
   see [todo.md](todo.md) §3.

Both reconstruction tutorials are dual-language, both lanes runnable, both scripts
version-controlled and transcluded into the pages so neither can rot.

**Exit criterion:** a reader can complete a full reconstruction in Python, following the
site, without consulting the source; and every code block on the page comes from a script
that CI executes.

## Phase 4 — Maturity, and what 1.0 means

*Draws evidence from phases 1–3, but may run alongside phase 3 (decision D5). Gates 1.0.*

* Coverage reporting with a floor, and a first pass over the modules no test obviously
  touches.
* A per-function maturity table — **stable** / **provisional** / **experimental** — against
  a written rubric, published on the site.
* The known functional gaps stated in one place (ftINIT metabolomics, no Excel import, no
  MetaCyc, no dynamic FBA, Gurobi in practice for genome-scale work).
* API pages and docstrings for the subsystems that have neither (`confidence`, `biomass`,
  `conditions`, `curation`, `annotation`, `manifest`, `data`, `binaries`) — thin docstrings
  here render as thin pages on the site.
* Resolve the open 💡 proposals in [improvements.md](improvements.md): implement or formally
  drop.
* A stability policy: what counts as public API, the deprecation window, and how `notes` /
  `annotation` key names are versioned.

**Exit criterion:** every public function carries a maturity label backed by evidence, and
the API surface is one we are willing to freeze.

## Continuous

Not phased — these run alongside:

* Land or close the open PRs ([#75](https://github.com/SysBioChalmers/raven-toolbox/pull/75),
  [#85](https://github.com/SysBioChalmers/raven-toolbox/pull/85)).
* Triage the 14 untracked scripts in `scripts/` — several are the cross-language drivers
  phase 2 needs.
* Track the upstream optlang/HiGHS and GLPK blockers.
* Back-port queue to MATLAB RAVEN: `assignCompartments.m`, then FS4 / B2 / G1 / G5.

---

## Decisions this roadmap is waiting on

| # | Decision | Blocks | Status |
|---|---|---|---|
| **D1** | **Vehicle for the new tutorials** — **both routes**: reconstruction by homology for a non-conventional yeast (first), and *de novo* from KEGG for a small prokaryote (second). Neither depends on SLIME lipid methods. Still to pick within this: the specific organisms, templates, and the published phenotype each is validated against — the way *H. polymorpha*'s methylotrophy anchors the existing protocol. | Phase 3 | **Decided** (2026-08-20); organisms TBD |
| **D2** | **How Python-only material presents** on a site whose convention is linked MATLAB/Python tabs defaulting to MATLAB. Compartmentalisation and confidence tracking have no MATLAB counterpart. Badge-and-explain, or a separate grouping? | Phase 3 | **Open** |
| **D3** | **SLIME lipid curation** — out of scope for the new tutorial while the parallel method work proceeds. | — | **Decided** (2026-08-20) |
| **D4** | **Where the tutorial's runnable scripts live** — a dedicated tutorial repository, the model repository, or this repo's `scripts/`. Determines what CI can execute. | Phase 3 | **Open** |
| **D5** | **What anchors 1.0** — nothing external. 1.0 ships when phase 4's exit criteria are met, so phases 3 and 4 may overlap. | Phase 4 | **Decided** (2026-08-20) |

Record D1–D2 and D4 in the raven-docs decisions log, since they shape the site.
