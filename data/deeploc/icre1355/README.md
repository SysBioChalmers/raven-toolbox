# Chlamydomonas (iCre1355) protein sequences for DeepLoc 2.1

Ready-to-upload protein FASTA for the genes of **iCre1355**, a genome-scale model of the green alga
*Chlamydomonas reinhardtii* — a second, fully independent non-yeast eukaryote (different lineage, and
it predates DeepLoc 2). Like AraCore it exercises the **chloroplast/plastid** yeast lacked, and adds
the richest algal compartmentation of the candidates (chloroplast, thylakoid lumen, mitochondrion,
nucleus, Golgi, glyoxysome, flagellum, eyespot, extra-organism).

| file | sequences |
|---|---|
| `iCre1355_proteins_001.fasta` … `_003.fasta` | 500 × 2 + 368 = **1368** |

Headers are the model's gene id (Phytozome/JGI v5.5 transcript id, e.g. `Cre01.g000350.t1.1`).

## Provenance

* **Genes:** iCre1355 — Imam S *et al.*, *Plant J.* 84(6):1239–1256 (2015),
  [doi:10.1111/tpj.13059](https://doi.org/10.1111/tpj.13059);
  [github.com/baliga-lab/Chlamy_model_iCre1355](https://github.com/baliga-lab/Chlamy_model_iCre1355)
  (`iCre1355_SBML_Matlab_files/iCre1355_auto.xml`). 1460 gene ids — **1368/1460** matched a sequence.
* **Sequences:** EnsemblPlants release-58, *C. reinhardtii* assembly **v5.5**
  (`Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.pep.all.fa.gz`), matched by the
  `gene_symbol` field (the v5.5 `Cre…` id, gene level), longest transcript per gene. Fetched 2026-06-23.

## Notes on preparation

* The model's SBML stores GPRs in the legacy `notes` element (no fbc), so cobra mangles them into the
  gene ids; the real gene ids were recovered by tokenising the GPR strings.
* **UniProt does not index the `Cre…` v5.5 ids**, so sequences come from EnsemblPlants (which keeps the
  v5.5 id in `gene_symbol`) rather than UniProt — a different source from the AraCore/Human-GEM inputs.
* The 92 unmatched genes are 41 organelle-genome-encoded proteins (`ChreCp…` chloroplast, `ChrepMp…`
  mitochondrion — absent from the nuclear proteome) and ~51 renamed/retired v5.5 gene models.

## How to use

Upload each file to [DeepLoc 2.1](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/) (≤500/run),
download the CSVs, then benchmark against iCre1355's curated compartments (DeepLoc maps to
chloroplast→plastid, cytosol, mitochondrion, nucleus, Golgi, glyoxysome→peroxisome, extra-organism;
flagellum/eyespot/thylakoid-lumen are out of DeepLoc's scope).
