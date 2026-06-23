# yeast-GEM protein sequences for DeepLoc 2.1

Ready-to-upload protein FASTA for **every gene in yeast-GEM**, for predicting subcellular
localisation with [DeepLoc 2.1](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/).

| file | sequences |
|---|---|
| `yeast-GEM_proteins_001.fasta` | 500 |
| `yeast-GEM_proteins_002.fasta` | 500 |
| `yeast-GEM_proteins_003.fasta` | 143 |

Each FASTA header is the **gene id** (ORF / ordered-locus name, e.g. `YNR001C`), so DeepLoc's
`Protein_ID` output column lines up with the model — and with
`raven_toolbox.localization.load_deeploc` — directly, no remapping.

The **DeepLoc 2.1 results** for these sequences are committed alongside (`yeast-GEM_deeploc_001.csv`,
`…_002.csv`, `…_003.csv` — one per FASTA chunk) and benchmarked against yeast-GEM's curated
compartments in [`docs/studies/deeploc_yeast_benchmark.md`](../../docs/studies/deeploc_yeast_benchmark.md)
(regenerate with `scripts/benchmark_deeploc_yeast.py`).

## How to use

1. Upload each file to the [DeepLoc 2.1 web server](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/)
   — it accepts **max 500 sequences per submission**, which is why the set is split into three. (For
   the downloadable standalone tool, which has no limit, `cat`-ing the three files gives one input.)
2. Download each result as CSV.
3. Load it back into a `gene × compartment` score table:

   ```python
   from raven_toolbox.localization import load_deeploc, DEFAULT_COMPARTMENT_MAP
   scores = load_deeploc("deeploc_result.csv", compartment_map=DEFAULT_COMPARTMENT_MAP)
   # ... then predict_localization(model, scores, reactions_to_relocate=...)
   ```

## Provenance

* **Genes:** yeast-GEM `v8.7.1-57-g9376ed7` (commit `9376ed7`), 1143 genes — 1143/1143 had a reviewed
  sequence (0 missing).
* **Sequences:** UniProtKB reviewed (Swiss-Prot), *Saccharomyces cerevisiae* S288C (taxon `559292`),
  fetched 2026-06-23.
* **Generated with:**

  ```
  python scripts/prepare_deeploc_yeast.py \
      --yeast-gem <path>/yeast-GEM/model/yeast-GEM.xml \
      --out data/deeploc/yeast-GEM_proteins.fasta
  ```

Regenerate against a newer model or organism with that script (see
`raven_toolbox.localization.prepare_deeploc_input`).

## Cross-species inputs (generalisation tests)

To check the yeast-GEM result is not an artefact of yeast-GEM's curation, the same pipeline is
prepared for two independent non-yeast eukaryotes (run DeepLoc on these too, then benchmark):

* [`aracore/`](aracore/) — *Arabidopsis* AraCore: a fully independent plant model that exercises the
  **chloroplast/plastid** yeast lacked (stringent: plant is far from DeepLoc's training).
* [`humangem/`](humangem/) — Human-GEM: a human positive control. **Note the circularity caveat** —
  15% of its gene compartments were assigned by DeepLoc2 and must be excluded.
