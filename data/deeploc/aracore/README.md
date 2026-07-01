# AraCore (Arabidopsis) protein sequences for DeepLoc 2.1

Ready-to-upload protein FASTA for the genes of **AraCore**, a curated *Arabidopsis thaliana* core
metabolism model — a cross-kingdom, fully independent test of whether the DeepLoc-vs-yeast-GEM
benchmark generalises (different kingdom, different lab, and AraCore predates DeepLoc 2, so its
compartments are not contaminated by DeepLoc).

| file | sequences |
|---|---|
| `AraCore_proteins_001.fasta` | 500 |
| `AraCore_proteins_002.fasta` | 160 |

Headers are the **AGI locus id** (e.g. `AT1G06680`), the model's gene id.

**Why this model.** It distinguishes **Chloroplast/plastid** (`h`) — the organelle yeast-GEM could
not exercise at all — plus Cytosol (`c`), Mitochondrion (`m`), Peroxisome (`p`), and the lumen/inter-
membrane sub-compartments. Plant proteins are under-represented in DeepLoc's human/animal-heavy
training, so this is a *stringent* generalisation test (unlike a human model, where high agreement
would just confirm DeepLoc on its own training distribution).

## Provenance

* **Genes:** AraCore v2.0 — Arnold A, Nikoloski Z, *Plant Physiol.* 165(3):1380–1391 (2014),
  [doi:10.1104/pp.114.235358](https://doi.org/10.1104/pp.114.235358); maintained at
  [github.com/pwendering/ArabidopsisCoreModel](https://github.com/pwendering/ArabidopsisCoreModel)
  (`AraCore_v2_0/AraCore_v2_0.xml`). 706 genes — **660/706** had a reviewed UniProt sequence (46
  lacked one, mostly TrEMBL-only loci).
* **Sequences:** UniProtKB reviewed, *A. thaliana* (taxon `3702`), matched by AGI locus
  (UniProt ordered-locus name, case-insensitive), fetched 2026-06-23.

## How to use

Upload each file to [DeepLoc 2.1](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/) (≤500/run),
download the result CSVs, then benchmark against AraCore's curated compartments (the model's reaction
compartmentalisation; AGI gene → its reactions' compartments).
