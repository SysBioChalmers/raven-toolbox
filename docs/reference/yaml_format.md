# RAVEN / cobrapy YAML model format

This document describes the YAML format produced and consumed by

- **cobrapy** ([`cobra.io.{load,save}_yaml_model`](https://github.com/opencobra/cobrapy))
- **raven-toolbox** (`raven_toolbox.io.yaml.{read,write}_yaml_model`, see [API](api/index.md))
- **RAVEN MATLAB** (`readYAMLmodel.m` / `writeYAMLmodel.m` in the [RAVEN repo](https://github.com/SysBioChalmers/RAVEN/tree/feat/yeast-gem-shared/io))

The same file can be round-tripped through any of the three. cobrapy is the canonical core; raven-toolbox and RAVEN MATLAB add namespaced extensions (RAVEN curation fields, MIRIAM cross-refs already covered by cobrapy's `annotation`, and the GECKO `ec-*` sections) without disturbing the cobra-known shape.

---

## At a glance

```yaml
!!omap
- metaData: !!omap
    - id: yeastGEM_develop
    - name: The Consensus Genome-Scale Metabolic Model of Yeast
    - version: 9.0.0
    - date: 2026-05-27
    - taxonomy: taxonomy/559292
- metabolites:
    - !!omap
      - id: s_0001
      - name: ATP
      - compartment: c
      - charge: -4
      - formula: C10H16N5O13P3
      - annotation: !!omap
          - kegg.compound: C00002
          - smiles: "[O-]P(=O)([O-])OP(=O)([O-])O..."
      - inchis: InChI=1S/C10H16N5O13P3/...
      - deltaG: -2768.1
- reactions:
    - !!omap
      - id: r_0001
      - name: hexokinase
      - metabolites: !!omap
          - s_0001: -1.0
          - s_0568: -1.0
          - s_0394: 1.0
          - s_0423: 1.0
      - lower_bound: 0.0
      - upper_bound: 1000.0
      - gene_reaction_rule: YGL253W or YCL040W or YFR053C
      - subsystem: Glycolysis / Gluconeogenesis
      - notes: "MetaNetX ID curated (PR #220)"
      - annotation: !!omap
          - kegg.reaction: R00299
          - sbo: SBO:0000176
      - eccodes: 2.7.1.1
      - deltaG: -17.39
      - confidence_score: 2.0
- genes:
    - !!omap
      - id: YGL253W
      - name: HXK2
      - annotation: !!omap
          - uniprot: P04807
- compartments: !!omap
    - c: cytoplasm
    - e: extracellular
```

Three structural rules are non-obvious and worth pointing out before the field-by-field detail:

1. The whole document is one **ordered mapping** — `!!omap` — at the root. Every nested map that should preserve key order is also `!!omap` (metaData, each metabolite / reaction / gene entry, `annotation`, `metabolites`, `compartments`, and the ec sections).
2. Each metabolite, reaction, and gene is **one `- !!omap` element** of a list. Inside that mapping, every field is written as `- key: value`. This is cobrapy's native shape and is what RAVEN MATLAB's reader keys off.
3. Strings are **unquoted by default**; quotes appear only when YAML would otherwise misparse the value (leading `-`, `[`, `?` or `:`; embedded `: ` or ` #`; values that look like `true` / `false` / `null`).

---

## Top-level layout

```
!!omap
- metaData: !!omap    # optional; RAVEN extension
- metabolites:        # required
- reactions:          # required
- genes:              # required (may be `genes: []`)
- compartments: !!omap # required
- gecko_light: <bool> # optional; GECKO extension
- ec-rxns:            # optional; GECKO extension
- ec-enzymes:         # optional; GECKO extension
```

| Key | Required | Source | Notes |
|-----|----------|--------|-------|
| `metaData` | optional | RAVEN | Provenance block. Holds `id`, `name`, `version`, `date`, `taxonomy`, optionally `givenName` / `familyName` / `email` / `organization` / `note` / `sourceUrl`, plus `defaultLB` / `defaultUB`. Cobrapy ignores this block (no semantic loss for the core model). |
| `metabolites` | yes | cobra core | Ordered list of `- !!omap` entries. |
| `reactions` | yes | cobra core | Ordered list of `- !!omap` entries. |
| `genes` | yes | cobra core | Ordered list; may be `genes: []` for a model with no genes. |
| `compartments` | yes | cobra core | `!!omap` of `<code>: <full name>`. |
| `gecko_light` | optional | GECKO | Scalar boolean. Cobrapy / raven-toolbox emit this at the top level; the older spelling `geckoLight` inside `metaData` is still accepted on read. |
| `ec-rxns` | optional | GECKO | Per-reaction kcat / source / enzymes coupling table. |
| `ec-enzymes` | optional | GECKO | Per-enzyme MW / sequence / concentration table. |

Cobrapy writes `id` / `name` / `version` at the root level instead of inside `metaData`. The RAVEN readers accept both placements; the RAVEN writers normalize to the `metaData` form.

---

## Metabolite entry

Field order (cobra-core first, then RAVEN extensions):

```yaml
- !!omap
  - id: s_0001               # required
  - name: ATP                # cobra
  - compartment: c           # cobra
  - charge: -4               # cobra (number)
  - formula: C10H16N5O13P3   # cobra
  - notes: "free-text"       # cobra
  - annotation: !!omap       # cobra (MIRIAM + smiles)
      - kegg.compound: C00002
      - chebi:
          - CHEBI:15422
          - CHEBI:30616
      - sbo: SBO:0000247
      - smiles: "OC1=NC..."  # quoted when it contains [ ] : etc.
  - inchis: "InChI=1S/..."   # RAVEN extension
  - deltaG: -2768.1          # RAVEN extension
  - metFrom: KEGG            # RAVEN extension
```

Cobrapy emits exactly the first seven keys (the cobra-core block). raven-toolbox and RAVEN MATLAB additionally emit `inchis`, `deltaG`, and `metFrom` when those fields are populated. On read, cobrapy puts the RAVEN extensions on the metabolite as attribute fall-through; raven-toolbox captures them into `metabolite.notes` (keyed by their YAML name); RAVEN MATLAB stores them on `model.inchis` / `model.metDeltaG` / `model.metFrom`.

Annotation entries with multiple values are emitted as a YAML list (`chebi:` then several `-` items). Single-value entries are emitted inline (`kegg.compound: C00002`). SMILES strings live inside the annotation block under the `smiles` key — not as a top-level metabolite field, which is the historical RAVEN MATLAB shape and is still accepted on read for backward compatibility.

---

## Reaction entry

```yaml
- !!omap
  - id: r_0001                                    # required
  - name: hexokinase                              # cobra
  - metabolites: !!omap                           # cobra (sorted by met id)
      - s_0001: -1.0
      - s_0394: 1.0
  - lower_bound: 0.0                              # cobra (number)
  - upper_bound: 1000.0                           # cobra (number)
  - gene_reaction_rule: YGL253W or YCL040W        # cobra
  - objective_coefficient: 1                      # cobra; omitted when 0
  - subsystem: Glycolysis / Gluconeogenesis       # cobra
  - notes: "MetaNetX ID curated (PR #220)"        # cobra
  - annotation: !!omap                            # cobra
      - kegg.reaction: R00299
      - sbo: SBO:0000176
  - eccodes:                                      # RAVEN extension
      - 2.7.1.1
      - 2.7.1.2
  - references: "PMID:12345"                      # RAVEN extension
  - rxnFrom: KEGG                                 # RAVEN extension
  - deltaG: -17.39                                # RAVEN extension
  - confidence_score: 2.0                         # RAVEN extension
```

Some fields are conditional:

- `objective_coefficient` is only written when non-zero (cobrapy convention).
- The `metabolites` block uses `!!omap []` (flow-style empty omap) when the reaction has no metabolites — this keeps the file a valid YAML 1.2 document.
- `eccodes` is written inline (`eccodes: 2.7.1.1`) when there is exactly one code, and as a list when there are several. Same for `references`.

**Notes key naming.** Cobrapy and the current raven-toolbox / RAVEN MATLAB writers use **`notes`**. Pre-`feat/yeast-gem-shared` yeast-GEM files used `rxnNotes`; both readers accept that as a legacy alias.

**Bounds typing.** Bounds are emitted as floats with an explicit decimal point (`1000.0`, `-1000.0`), matching Python's float repr and cobrapy's output.

---

## Gene entry

```yaml
- !!omap
  - id: YGL253W              # required
  - name: HXK2               # cobra; omitted when empty
  - annotation: !!omap       # cobra
      - uniprot: P04807
      - ncbigene: 856421
  - protein: P04807          # RAVEN extension
```

Empty names (`name: ''`) are not emitted (matches RAVEN MATLAB's historical behavior).

---

## Compartments

```yaml
- compartments: !!omap
    - c: cytoplasm
    - e: extracellular
    - m: mitochondrion
```

Just an `!!omap` of `<short code>: <human-readable name>` pairs. Compartments don't carry their own MIRIAMs in the current format.

---

## metaData (RAVEN extension)

```yaml
- metaData: !!omap
    - id: yeastGEM_develop
    - name: The Consensus Genome-Scale Metabolic Model of Yeast
    - version: 9.0.0
    - date: 2026-05-27
    - defaultLB: -1000.0
    - defaultUB: 1000.0
    - givenName: Eduard
    - familyName: Kerkhoven
    - email: eduardk@chalmers.se
    - organization: Chalmers University of Technology
    - taxonomy: taxonomy/559292
    - note: "Saccharomyces cerevisiae - strain S288C"
    - sourceUrl: https://github.com/SysBioChalmers/yeast-GEM
```

Pure provenance. Cobrapy ignores the block; raven-toolbox keeps the verbatim dictionary on `model.notes['metaData']` and additionally lifts `id` / `name` / `version` to `model.id` / `model.name` / `model.notes['version']` so cobra-shape accessors find them. RAVEN MATLAB populates `model.id` / `model.name` / `model.version` / `model.annotation.*` from the same fields.

`date` is preserved across round-trips when present on the model; otherwise the writer fills in `YYYY-MM-DD` of the current date.

---

## GECKO sections

For enzyme-constrained models, three additional top-level keys carry the EC layer:

```yaml
- gecko_light: false        # true for the "light" formulation
- ec-rxns:
    - !!omap
      - id: r_0001
      - kcat: 25.3
      - source: brenda
      - notes: ""
      - eccodes: 2.7.1.1
      - enzymes: !!omap
          - P04807: 1.0
- ec-enzymes:
    - !!omap
      - genes: YGL253W
      - enzymes: P04807
      - mw: 53942
      - sequence: "MVHLGPK..."
      - concs: .nan
```

These map onto `model.ec` in RAVEN MATLAB and `raven_toolbox.io.ec_data.EcData` (attached as `model.ec`) in raven-toolbox. Cobrapy ignores the sections.

The older spelling `geckoLight` inside `metaData` is also accepted on read.

---

## Annotations

The `annotation` block uses MIRIAM-style namespace keys. Cobrapy treats the block as a free-form dictionary; raven-toolbox preserves it verbatim through `cobra.Metabolite.annotation` / `Reaction.annotation` / `Gene.annotation`; RAVEN MATLAB maps it to `model.metMiriams` / `rxnMiriams` / `geneMiriams`.

- A single value is written inline: `kegg.compound: C00002`.
- Multiple values are written as a YAML list:

  ```yaml
  - chebi:
      - CHEBI:15422
      - CHEBI:30616
  ```

- The `smiles` key inside a metabolite's `annotation` carries the SMILES string (cobrapy convention). RAVEN MATLAB historically emitted `smiles` as a metabolite top-level field; both readers still accept that, but writes are normalized to the annotation block.
- The `sbo` key carries the Systems Biology Ontology term assigned by `assignSBOterms` / `add_sbo_terms`.

---

## Numbers, strings, quoting

**Numbers.** Whole-number floats are written with an explicit `.0` (`1000.0`, `-1000.0`, `0.0`). Other floats use up to 15 significant digits (`-17.39`, `-2768.1`). `NaN` is encoded as `.nan`; `+Inf` / `-Inf` as `.inf` / `-.inf` (YAML 1.2 conventions).

**Strings.** Default style is bare (no quotes). The writer falls back to double-quoted style when the value:

- starts with `-`, `?`, `:`, or any flow indicator (`[`, `]`, `{`, `}`, `,`, `&`, `*`, `!`, `|`, `>`, `%`, `@`, `` ` ``, `#`);
- contains `: ` (would otherwise be parsed as a key/value), ` #` (comment), or one of `[`, `]`, `{`, `}`;
- has leading or trailing whitespace;
- spells a YAML reserved word case-insensitively (`true`, `false`, `null`, `yes`, `no`, `on`, `off`, `~`).

In a double-quoted string, only `\` and `"` are escaped. Other characters (including Unicode and newlines if the underlying model permitted them) are passed through.

---

## Tooling interoperability matrix

| File written by ↓ \ Reader → | cobrapy | raven-toolbox | RAVEN MATLAB |
|---|---|---|---|
| cobrapy (`save_yaml_model`) | full | full + extras land in `notes` via attribute fall-through | works for root-level `id` / `name` / `version` (added in this release) |
| raven-toolbox (`write_yaml_model`) | core (no `metaData`-derived `id`); RAVEN extras live as unknown top-level keys but don't break parsing | full | full |
| RAVEN MATLAB (`writeYAMLmodel`) | core (no `metaData`-derived `id`); RAVEN extras land via attribute fall-through | full | full |

"Full" = every field read back into its canonical position on the model object; "core" = cobrapy-known fields, RAVEN extensions ignored or kept on the object as attribute fall-through (`reaction.eccodes` etc., not re-emitted on save). A round-trip through cobrapy is therefore **lossy for RAVEN extensions** — only the core fields survive `cobrapy.load → cobrapy.save`. Round-trips through raven-toolbox or RAVEN MATLAB are lossless.

---

## What round-tripping looks like

Loading `yeast-GEM.yml` (2748 metabolites, 4102 reactions, 1143 genes) and re-writing it through any of the three tools preserves every documented piece of content:

| Count | After round-trip |
|---|---|
| metabolites | 2748 / 2748 |
| reactions | 4102 / 4102 |
| genes | 1143 / 1143 |
| reactions with eccodes | 2411 |
| reactions with deltaG | 3984 |
| metabolites with deltaG | 2696 |
| metabolites with SMILES | 1788 |
| reactions with notes (rxnNotes) | 1443 |

(Cobrapy round-trips give 2748 / 4102 / 1143 for the core but drop the RAVEN extensions in the rightmost column — that's the documented loss.)

---

## What changed in `feat/yeast-gem-shared`

- raven-toolbox writer no longer drops `!!omap` tags (was producing files RAVEN MATLAB's reader couldn't load).
- raven-toolbox now preserves `eccodes` and accepts the legacy `rxnNotes` reaction key on read.
- RAVEN MATLAB writer reorders metabolite / reaction fields to match cobrapy.
- RAVEN MATLAB writer renames the reaction `rxnNotes` key to `notes` and emits SMILES inside the annotation block (still accepts both shapes on read).
- RAVEN MATLAB writer's `preserveQuotes` default is now `false`; values that need quoting (SMILES with `[O-]`, leading flow indicators, booleans, `: `-containing strings) are quoted defensively per value.
- RAVEN MATLAB writer emits whole-number bounds as `1000.0` (matches cobrapy / Python float repr) instead of `1000`.
- RAVEN MATLAB reader accepts cobrapy's root-level `id` / `name` / `version` / `gecko_light`, the `!!omap`-tagged `metaData` header, and `notes` (canonical) in addition to `rxnNotes` (legacy).
- Empty `reaction.metabolites` blocks are emitted as `!!omap []` (valid YAML 1.2) rather than an empty `!!omap` with no value.
- Document-start marker `---` dropped to match cobrapy's bare `!!omap` root.

These changes are byte-stable for cobrapy and raven-toolbox users; existing yeast-GEM YAML files continue to load. The first time a yeast-GEM curation pass rewrites the file with the new MATLAB writer, the diff will look large (because of the reordering and quote-style changes) but the model content is unchanged.
