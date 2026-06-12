# Resource resolvers

Version-pinned, SHA256-verified registries that fetch the external resources raven-toolbox
depends on but does not vendor.

## `raven_toolbox.binaries`

External-tool resolver for BLAST+ / DIAMOND / HMMER (release-ZIP registry).

```{eval-rst}
.. automodule:: raven_toolbox.binaries
   :members:
```

## `raven_toolbox.data`

Data-bundle resolver (KEGG artefacts and template-model data).

```{eval-rst}
.. automodule:: raven_toolbox.data
   :members:
```

## `raven_toolbox.manifest`

Loads a shared [data/binary manifest](../../maintenance/data_manifest.md) into the two
registries above (and is consulted lazily via `$RAVEN_PYTHON_MANIFEST`).

```{eval-rst}
.. automodule:: raven_toolbox.manifest
   :members:
```
