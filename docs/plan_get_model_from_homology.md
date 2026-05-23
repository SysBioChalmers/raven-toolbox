# Design plan: `get_model_from_homology` (port of RAVEN `getModelFromHomology`)

Detailed design for the core of homology reconstruction (Phase 3a). Companion to
PLAN.md §2.3a. Goal: faithful *intent*, but a clearer, more robust implementation —
RAVEN's own comments flag several soft spots, and the user asked for logic
improvements.

---

## 1. What it is trying to do (the intent)

Build a draft GEM for a **new organism** by transferring reactions from one or more
**template GEMs**, using **orthology** (bidirectional protein BLAST/DIAMOND hits) to
map template genes → new-organism genes and rewrite each reaction's GPR accordingly.

A reaction is transferred if at least one of its genes has an ortholog in the new
organism; its GPR is rewritten so template genes become their new-organism orthologs.
When several templates are given, an optional **preferred order** makes each gene's
reactions come from a single (highest-priority) template; otherwise reactions from all
templates are merged. Metabolites are unified across templates by **name[compartment]**.

## 2. RAVEN's algorithm, in steps

1. Strip cross-species gene metadata from templates (`geneShortNames`, `geneMiriams`,
   `proteins`, `geneFrom`); standardize template grRules.
2. Validate: gene ids without stray `()`; blast `fromId`s match template ids; ≥5 % of a
   template's genes appear in the hits (else the FASTA/model id styles disagree).
3. **Filter hits** by `evalue ≤ maxE`, `aligLen ≥ minLen`, `identity ≥ minIde`.
4. Drop template reactions with no genes (and genes with no reactions).
5. `onlyGenesInModels`: drop hits whose template gene isn't in a model.
6. `strictness==3`: reduce to best-by-**E-value** hits per `fromGene`.
7. Build sparse matrices `allTo` (new→template hits) and `allFrom` (template→new) per
   template, indexed by gene position.
8. **Ortholog map** = per strictness: `1`/`3` → `allTo & allFrom` (reciprocal);
   `2` → just `allFrom` or `allTo` (one-directional, per `mapNewGenesToOld`).
9. Drop mappings to genes not in the model; simplify.
10. Per template, keep reactions associated with a mapped gene (keeping AND-complexes
    where *any* subunit mapped).
11. `preferredOrder`: remove genes already claimed by an earlier template so each gene's
    reactions come from one template.
12. **Rewrite GPRs** by `regexprep` string substitution: replace each template gene with
    its new ortholog(s) — `(new1 or new2)` if several; template genes with no ortholog
    are renamed `OLD_<modelid>_<gene>` (kept so AND-complexes survive).
13. `mergeModels(..., 'metNames')`; then `regexprep` away `OLD_…` genes that ended up in
    `or` relations; set id/name/notes/confidence; standardize; delete unused genes.

## 3. Proposed Python design

```python
def get_model_from_homology(
    models,                     # list[cobra.Model], or one model
    hits,                       # bidirectional hits DataFrame (run_blast / make_ortholog_hits)
    model_for,                  # target organism id
    *,
    preferred_order=None,       # list[model_id]; if set, each gene's reactions from one model
    bidirectional=True,         # require reciprocal hits (RBH-style)        [improvement A]
    best_hits_only=False,       # keep only best-scoring hit per gene first  [improvement A]
    map_direction="new_to_old", # used only when bidirectional=False
    score="bitscore",           # best-hit criterion: "bitscore" | "evalue" [improvement D]
    complex_policy="flag",      # AND-subunits lacking orthologs: flag|keep|drop [improvement C]
                                # default "flag" = RAVEN-compatible (OLD_<model>_<gene>)
    only_genes_in_models=False,
    max_evalue=1e-30, min_align_len=200, min_identity=40,
) -> HomologyResult            # .model + .gene_map + .reaction_sources  [improvement F]
```

**Data flow (all on the hits DataFrame + cobra GPR AST, no sparse-matrix juggling):**

1. Validate + filter hits (`max_evalue`/`min_align_len`/`min_identity`) — a DataFrame mask.
2. Optionally reduce to best hit per gene by `score` (groupby-idxmax).
3. Build the **ortholog map** `new_gene -> {model_id: {template_gene, ...}}`:
   - `bidirectional`: inner-join the two directions on (new_gene, template_gene) — a pandas
     merge, i.e. reciprocal hits. (`bidirectional` + `best_hits_only` = reciprocal best hits.)
   - else: take the one direction (`map_direction`).
4. For each template, for each reaction: rewrite the GPR on the **cobra GPR AST**
   (improvement B) — substitute each leaf gene by the OR of its new orthologs; handle
   unmapped leaves per `complex_policy`; keep the reaction iff ≥1 leaf mapped.
5. Resolve `preferred_order`: a new gene's reactions come from the first template (in order)
   that maps it (improvement E — a dict lookup, not matrix index math).
6. Transfer the surviving reactions into the draft, unifying metabolites by `name[comp]`
   (reuse `merge_models` / `add_reactions_from_model`); record provenance (improvement F).
7. Finalize: id=`model_for`, name, per-reaction note/confidence, prune unused genes,
   lint GPRs (`find_non_dnf_grrules`).

`make_ortholog_hits` (PLAN §2.3a) feeds the same DataFrame, so the whole thing is
testable without BLAST.

## 4. Logic improvements (the point of this exercise)

| # | Improvement | Why it's better than RAVEN |
|---|---|---|
| **H1** | **Two orthogonal params** (`bidirectional`, `best_hits_only`) replace the overloaded `strictness` 1/2/3. | RAVEN's `strictness` conflates two independent axes — *directionality* (one-way vs reciprocal) and *best-hit filtering*. Splitting them is clearer and exposes all 4 combinations (incl. one-way best-hit). RAVEN mapping documented: `1`→bidir; `2`→one-way; `3`→bidir+best-hits = **reciprocal best hits (RBH)**. |
| **H2** | **AST-based GPR rewriting** (cobra `GPR`) instead of `regexprep` on rule strings. | RAVEN's author flags the string substitution as uncertain ("I hope that it's ok…") and needs a follow-up regex to clean `OLD_… or` leftovers. Substituting leaves on the parsed boolean tree is robust (no partial-match hazards, no cleanup pass) and reuses the AST machinery already in `expand_model`/`find_non_dnf_grrules`. |
| **H3** | **Explicit complex policy** for AND-subunits lacking an ortholog: `keep` (drop the unmapped subunit, transfer the reaction), `drop` (require a fully-mapped AND-clause, else drop the reaction), `flag` (RAVEN's `OLD_<model>_<gene>` placeholder). | RAVEN's behaviour is the implicit, fragile `OLD_`+regex path its comments distrust. Making it an explicit, AST-correct policy (OR = keep isozyme branches that mapped; AND = configurable) is principled and lets the user choose optimism vs strictness. `flag` preserves RAVEN compatibility. |
| **H4** | **Best-hit selection by `bitscore`** (default), E-value optional. | E-value depends on database size and ties at ~0 for strong hits; **bitscore is database-size-independent** and the standard criterion for reciprocal-best-hit orthology. RAVEN uses min E-value only. |
| **H5** | **DataFrame ortholog map** replaces the `allGenes`/`allTo`/`allFrom` sparse-matrix + `sub2ind` index juggling. | The reciprocal mapping and preferred-order resolution become a pandas merge + dict lookup — far less error-prone than the position-index gymnastics RAVEN itself annotates with uncertainty. |
| **H6** | **Provenance** in the result: per reaction, which template + which ortholog pairs supported it (`reaction.notes`), and a returned `gene_map`. | RAVEN returns only `hitGenes` (flat old/new lists) and a fixed note. Structured provenance aids curation and debugging. |

**RAVEN-compatibility:** accept a `strictness=` alias that sets `bidirectional`/`best_hits_only`
(1/2/3) and default `complex_policy="flag"` *iff* `strictness` is passed, so legacy calls reproduce
RAVEN behaviour; otherwise use the clearer defaults above.

## 5. Resolved: `complex_policy` default

**Decided: default `complex_policy="flag"`** (RAVEN-compatible) — for an AND-complex reaction with a
subunit lacking an ortholog, keep the reaction and mark the missing subunit (`OLD_<model>_<gene>`),
matching current RAVEN output. `keep` (drop the unmapped subunit) and `drop` (require a fully-mapped
AND-clause) remain available for users wanting cleaner or higher-confidence drafts. Implementation
note: even under `flag`, the GPR is rebuilt on the **AST** (improvement H2) — the `OLD_` marking and
the removal of `OLD_` genes left in `or` branches are done as AST operations, not the regex passes
RAVEN uses.

## 6. Test strategy

- **Core, no BLAST:** drive with `make_ortholog_hits` + small template models; assert reactions
  transfer, GPRs rewrite correctly (one-to-one, one-to-many isozymes, AND-complex under each
  `complex_policy`), `bidirectional`/`best_hits_only` select the right hits, `preferred_order`
  routes a gene's reactions to one template, metabolites unify by name[comp].
- **Edge cases:** unmapped subunit policies; gene mapping to multiple orthologs (`(a or b)`);
  reaction supported by only one isozyme; the ≥5 % overlap validation; reciprocal vs one-way.
