# Human-GEM protein sequences for DeepLoc 2.1

Ready-to-upload protein FASTA for the genes of **Human-GEM** (Human1), as a human cross-check of the
DeepLoc-vs-yeast-GEM benchmark. Human is the core of DeepLoc 2.x training, so this is a *lenient*
positive control (high agreement is expected and does not by itself prove generalisation — pair it
with the stringent plant model, `../aracore/`).

| file | sequences |
|---|---|
| `Human-GEM_proteins_001.fasta` … `_006.fasta` | 500 × 5 + 339 = **2839** |

Headers are the **Ensembl gene id** (`ENSG…`), the model's gene id.

## ⚠️ Circularity caveat — exclude DeepLoc2-sourced genes when benchmarking

Human-GEM's gene localisations carry a provenance column (`model/genes.tsv`, `compDataSource`):

| source | genes | use for a DeepLoc benchmark? |
|---|--:|---|
| SwissProt (± CellAtlas) | 1926 | ✅ independent |
| CellAtlas | 442 | ✅ independent |
| **DeepLoc2** | **439** | ❌ **circular — DeepLoc grading itself** |
| (none) | 41 | — |

So **439/2848 (15%) of Human-GEM gene compartments were assigned by DeepLoc2 itself.** A fair
benchmark must **drop rows where `compDataSource == "DeepLoc2"`** and score only the ~84%
SwissProt/CellAtlas-sourced genes.

## Provenance

* **Genes:** Human-GEM v2.0.0 — Robinson JL *et al.*, *Sci. Signal.* 13(624):eaaz1482 (2020),
  [doi:10.1126/scisignal.aaz1482](https://doi.org/10.1126/scisignal.aaz1482);
  [github.com/SysBioChalmers/Human-GEM](https://github.com/SysBioChalmers/Human-GEM) (tag `v2.0.0`).
  2848 genes — **2839/2848** mapped to a sequence.
* **Sequences:** ENSG → UniProt accession via the repo's `model/genes.tsv` (`geneUniProtID`), then
  UniProtKB reviewed sequences for *H. sapiens* (taxon `9606`), fetched 2026-06-23.

## How to use

Upload each file to [DeepLoc 2.1](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/) (≤500/run),
download the CSVs, then benchmark against `genes.tsv` `compartments` **after dropping the DeepLoc2-
sourced rows** (see caveat).
