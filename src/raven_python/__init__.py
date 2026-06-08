"""raven_python — Python counterpart of the RAVEN Toolbox, built on cobrapy.

raven_python reuses cobrapy for simulation, standard analyses, SBML I/O, and model
manipulation, and provides the RAVEN-specific functionality on top: de novo
reconstruction (KEGG / homology), context-specific modeling (tINIT / ftINIT),
metabolic task validation, connectivity gap-filling, omics integration (HPA),
sub-cellular localisation, N-model comparison, and the RAVEN-style I/O formats.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("raven-python")
except PackageNotFoundError:  # running from a source tree without an install
    __version__ = "0.0.0+unknown"
