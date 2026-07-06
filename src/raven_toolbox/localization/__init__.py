"""Sub-cellular localisation — predictor-agnostic, partial-update friendly.

:func:`predict_localization` is the score-driven MILP entry point; :func:`assign_compartments`
is the functionality-constrained variant (a flux-free placement master plus real-FBA certification,
with optional gap-fill, transport pruning and sound reaction-level multi-localisation). :func:`load_deeploc`,
:func:`load_mulocdeep` and :func:`load_compartments` parse modern predictor / evidence-database
outputs into the ``gene × compartment`` :class:`LocalizationScores` DataFrame the algorithm
consumes (pass :data:`DEFAULT_COMPARTMENT_MAP` to map labels to your model's compartment ids).

For sequence-based predictors (DeepLoc 2.1, MULocDeep), :func:`prepare_deeploc_input` writes a FASTA
of your model's gene sequences (fetched from UniProtKB, headers = gene ids) ready to run, closing the
loop with :func:`load_deeploc`.
"""
from raven_toolbox.localization.assign import (
    AssignmentProposal,
    GrowthCondition,
    apply_assignment,
)
from raven_toolbox.localization.certify import assign_compartments
from raven_toolbox.localization.curation import curation_priority
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
from raven_toolbox.localization.transport_evidence import (
    TransporterAnnotation,
    annotate_transporters,
    evidence_aware_transport_cost,
)
from raven_toolbox.localization.triage import (
    DEEPLOC_COMPARTMENT_TRUST,
    ReviewReport,
    triage_localization,
)

__all__ = [
    "AssignmentProposal",
    "GrowthCondition",
    "DEFAULT_COMPARTMENT_MAP",
    "DEEPLOC_COMPARTMENT_TRUST",
    "LocalizationProposal",
    "LocalizationResult",
    "LocalizationScores",
    "PreparedFasta",
    "ReviewReport",
    "TransporterAnnotation",
    "annotate_transporters",
    "apply_assignment",
    "apply_localization",
    "assign_compartments",
    "combine_scores",
    "curation_priority",
    "evidence_aware_transport_cost",
    "fetch_protein_sequences",
    "fetch_uniprot_localization",
    "load_compartments",
    "load_deeploc",
    "load_mulocdeep",
    "load_uniprot",
    "prepare_deeploc_input",
    "predict_localization",
    "triage_localization",
    "write_fasta",
]
