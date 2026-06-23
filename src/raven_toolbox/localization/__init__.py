"""Sub-cellular localisation — predictor-agnostic, partial-update friendly.

:func:`predict_localization` is the MILP entry point; :func:`load_deeploc`,
:func:`load_mulocdeep` and :func:`load_compartments` parse modern predictor / evidence-database
outputs into the ``gene × compartment`` :class:`LocalizationScores` DataFrame the algorithm
consumes (pass :data:`DEFAULT_COMPARTMENT_MAP` to map labels to your model's compartment ids).

For sequence-based predictors (DeepLoc 2.1, MULocDeep), :func:`prepare_deeploc_input` writes a FASTA
of your model's gene sequences (fetched from UniProtKB, headers = gene ids) ready to run, closing the
loop with :func:`load_deeploc`.
"""
from raven_toolbox.localization.predict import (
    LocalizationProposal,
    LocalizationResult,
    apply_localization,
    predict_localization,
)
from raven_toolbox.localization.scores import (
    DEFAULT_COMPARTMENT_MAP,
    LocalizationScores,
    combine_scores,
    fetch_uniprot_localization,
    load_compartments,
    load_deeploc,
    load_mulocdeep,
    load_uniprot,
)
from raven_toolbox.localization.sequences import (
    PreparedFasta,
    fetch_protein_sequences,
    prepare_deeploc_input,
    write_fasta,
)

__all__ = [
    "DEFAULT_COMPARTMENT_MAP",
    "LocalizationProposal",
    "LocalizationResult",
    "LocalizationScores",
    "PreparedFasta",
    "apply_localization",
    "combine_scores",
    "fetch_protein_sequences",
    "fetch_uniprot_localization",
    "load_compartments",
    "load_deeploc",
    "load_mulocdeep",
    "load_uniprot",
    "prepare_deeploc_input",
    "predict_localization",
    "write_fasta",
]
