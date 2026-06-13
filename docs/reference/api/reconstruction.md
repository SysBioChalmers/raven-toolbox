# `raven_toolbox.reconstruction`

De novo reconstruction from KEGG and protein homology (BLAST / DIAMOND).

## Homology

Homology-based reconstruction from template models (`getModelFromHomology`, BLAST / DIAMOND).

```{eval-rst}
.. automodule:: raven_toolbox.reconstruction.homology
   :members:
   :imported-members:
```

## KEGG

KEGG-based draft reconstruction (`getKEGGModelForOrganism` and friends): download → dump
parsing → HMM libraries (maintainer build steps), then the runtime model for a KEGG species.

```{eval-rst}
.. automodule:: raven_toolbox.reconstruction.kegg
   :members:
   :imported-members:
```
