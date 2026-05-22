# ravenpy

**Reconstruction, Analysis and Visualization of Metabolic Networks — in Python.**

`ravenpy` is a Python port of the [RAVEN Toolbox 2](https://github.com/SysBioChalmers/RAVEN)
(a MATLAB suite for genome-scale metabolic model reconstruction). Rather than re-implementing
the whole toolbox, `ravenpy` builds on [**cobrapy**](https://github.com/opencobra/cobrapy) for
everything cobrapy already does well (simulation, standard analyses, SBML I/O, model
manipulation) and focuses its own code on the functionality that is **unique to RAVEN**:

- **De novo reconstruction** from KEGG, MetaCyc, and protein homology (BLAST/DIAMOND).
- **Context-specific models** from omics data via tINIT / ftINIT.
- **Metabolic task** validation (`checkTasks` / `fitTasks`).
- **RAVEN-style gap-filling** against template models.
- **Omics integration** (Human Protein Atlas) and **subcellular localization** prediction.
- **RAVEN ⇄ COBRA model conversion** and RAVEN Excel/YAML I/O.

> **Status:** Pre-alpha. This repository currently contains the project scaffold and the
> [port plan](PLAN.md). See [PLAN.md](PLAN.md) for the full RAVEN→cobrapy functionality map and
> the phased roadmap.

## Design principle

The canonical in-memory object is a [`cobra.Model`](https://cobrapy.readthedocs.io) — ravenpy
functions consume and produce `cobra.Model` directly, with **no parallel RAVEN struct and no
`ravenCobraWrapper`-style adapter**. RAVEN-specific fields that cobrapy does not model natively
(e.g. `rxnMiriams`, `metDeltaG`, `rxnConfidenceScores`) live in cobra's `annotation` / `notes`
dictionaries. This avoids duplicating cobrapy's data model and keeps `ravenpy` interoperable with
the wider COBRA ecosystem. ravenpy's YAML I/O follows the cobrapy YAML standard and additionally
supports geckopy's enzyme-constrained extension keys so ecModels round-trip.

## Installation (development)

```bash
git clone https://github.com/SysBioChalmers/ravenpy
cd ravenpy
pip install -e ".[dev]"
```

## Relationship to the MATLAB RAVEN Toolbox

`ravenpy` is a derivative work of RAVEN and is therefore released under the same
**GPL-3.0-or-later** license. If you use it in scientific work, please cite the RAVEN 2 paper:

> Wang H, Marcišauskas S, Sánchez BJ, Domenzain I, Hermansson D, Agren R, Nielsen J,
> Kerkhoven EJ. (2018) RAVEN 2.0: A versatile toolbox for metabolic network reconstruction
> and a case study on *Streptomyces coelicolor*. PLoS Comput Biol 14(10): e1006541.

## License

[GPL-3.0-or-later](LICENSE)
