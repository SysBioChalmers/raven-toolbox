#!/usr/bin/env python3
"""Validate evidence-aware transport scoring on a curated GEM (yeast-GEM by default).

The claim behind :mod:`raven_toolbox.localization.transport_evidence` is that weighting the transport
penalty by transporter evidence turns an *indiscriminate* cut into a *selective* one: real carriers get
cheap (retained), unsupported transports pay the full prior (dropped). This script measures that on a
curated model whose transport reactions are the ground truth.

It annotates the organism's proteome (``annotate_proteome`` -> hmmsearch + diamond), scores every
metabolite's transport cost (``evidence_aware_transport_cost`` with ``default_substrate_of`` and
DeepLoc-derived carrier compartments), then runs the head-to-head:

* the essential carriers a blanket transport penalty drops -- the study in
  docs/studies/carvefungi_milp_benchmark.md showed a uniform penalty drops these curated, individually
  essential shuttles (2-oxoglutarate, 2-dehydropantoate, NADP(+)/NADPH, serine); their evidence-aware
  cost should be ~0 (retained);
* the **selective-cut** metric vs the blanket penalty: curated-transport rate among *evidenced*
  (would-keep) vs *unsupported* (would-drop) metabolites, across the evidence layers (coarse ->
  +ChEBI -> +sibling). A blanket penalty is indiscriminate (keep-rate == drop-rate == base rate);
  evidence is selective when keep-rate > drop-rate.

The proteome + DeepLoc predictions ship under ``data/deeploc/``; pass the curated model with ``--model``.
Usage: ``python scripts/validate_transport_evidence.py --model path/to/yeast-GEM.xml``. ASCII-only.
"""
from __future__ import annotations

import argparse
import tempfile
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")  # quiet cobra/libsbml import chatter before it loads
import cobra  # noqa: E402
import pandas as pd  # noqa: E402

from raven_toolbox.localization import (  # noqa: E402
    DEFAULT_COMPARTMENT_MAP,
    SubstrateOntology,
    annotate_proteome,
    default_substrate_of,
    evidence_aware_transport_cost,
    load_deeploc,
)

# yeast-GEM names of carriers whose cytosol<->mito shuttle a blanket transport penalty drops (breaking
# NADPH balance / the TCA cycle); the evidence-aware cost should retain them (cost ~ 0).
ESSENTIAL_CARRIERS = (
    "(S)-malate", "citrate", "2-oxoglutarate", "oxaloacetate", "NADP(+)", "NADPH",
    "L-serine", "2-dehydropantoate", "pyruvate", "phosphoenolpyruvate",
)


def _gene_compartments(deeploc_dir: Path, keep: set[str]) -> dict[str, set[str]]:
    """gene -> {top DeepLoc compartment} for each annotated transporter gene."""
    frames = [load_deeploc(csv, compartment_map=DEFAULT_COMPARTMENT_MAP).df
              for csv in sorted(deeploc_dir.glob("*_deeploc_*.csv"))]
    sdf = pd.concat(frames)
    return {g: {sdf.loc[g].astype(float).idxmax()} for g in sdf.index if g in keep}


def _report(cost: dict[str, float], label: str, curated: set[str], base_cost: float) -> None:
    """Selective-cut line: kept (evidenced) vs dropped, each with its curated-transport rate + recall."""
    allm = set(cost)
    cur = curated & allm
    keep = {b for b, c in cost.items() if c < base_cost}   # evidenced -> cheap -> kept
    drop = allm - keep                                     # unsupported -> full prior -> dropped
    keep_rate = len(keep & cur) / max(1, len(keep))
    drop_rate = len(drop & cur) / max(1, len(drop))
    recall = len(keep & cur) / max(1, len(cur))
    print(f"  [{label:14s}] kept {len(keep):4d} (curated {keep_rate:.0%})  "
          f"dropped {len(drop):4d} (curated {drop_rate:.0%})  recall {recall:.0%}")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True, type=Path, help="curated GEM (SBML), e.g. yeast-GEM.xml")
    ap.add_argument("--data-dir", type=Path, default=Path("data/deeploc"),
                    help="dir with *_proteins_*.fasta + *_deeploc_*.csv (default: data/deeploc)")
    ap.add_argument("--base-cost", type=float, default=0.5)
    ap.add_argument("--sibling-weight", type=float, default=0.5,
                    help="ChEBI sibling credit (0=off; >0 trades specificity for recall)")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args(argv)

    # 1. annotate the proteome (all FASTA chunks concatenated)
    fastas = sorted(args.data_dir.glob("*_proteins_*.fasta"))
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as tf:
        tf.write("".join(p.read_text() for p in fastas))
        proteome = Path(tf.name)
    ann = annotate_proteome(proteome, threads=args.threads)
    proteome.unlink(missing_ok=True)
    print(f"transporter genes annotated: {len(ann)} (from {len(fastas)} proteome chunk(s))")

    # 2. carrier compartments from DeepLoc; score under the blanket approach vs evidence-aware variants
    gene_comps = _gene_compartments(args.data_dir, set(ann))
    model = cobra.io.read_sbml_model(str(args.model))
    kw = dict(base_cost=args.base_cost, base_metabolite=lambda m: m.name)
    cost_coarse = evidence_aware_transport_cost(model, ann, gene_comps,
                                                substrate_of=default_substrate_of, **kw)
    ontology = SubstrateOntology.load()
    cost_chebi = evidence_aware_transport_cost(model, ann, gene_comps,
                                               substrate_of=default_substrate_of, ontology=ontology, **kw)
    cost_sib = evidence_aware_transport_cost(model, ann, gene_comps, substrate_of=default_substrate_of,
                                             ontology=ontology, sibling_weight=args.sibling_weight, **kw)

    # 3. the carriers a blanket penalty drops but curated yeast-GEM marks essential -- evidence-aware
    #    should retain them (cost ~ 0); a blanket penalty, having no evidence, cannot tell them apart.
    print(f"\n=== essential carriers a blanket penalty drops (base_cost {args.base_cost}; "
          "lower = evidenced, retained) ===")
    for name in ESSENTIAL_CARRIERS:
        if name in cost_chebi:
            flag = "retained" if cost_chebi[name] < args.base_cost else "NOT evidenced"
            print(f"  {name:22s} cost={cost_chebi[name]:.3f}  [{flag}]")

    # 4. head-to-head: is the transport cut selective? (curated yeast-GEM transports = ground truth)
    curated = {m.name for r in model.reactions if not r.boundary
               and len({mm.compartment for mm in r.metabolites}) > 1 for m in r.metabolites}
    base_rate = len(curated & set(cost_chebi)) / max(1, len(set(cost_chebi)))
    print("\n=== transport-cut selectivity vs a blanket penalty (curated yeast-GEM = ground truth) ===")
    print(f"  [{'blanket (no evidence)':14s}] kept-curated == dropped-curated == {base_rate:.0%}  "
          "(indiscriminate -- the reference-approach behaviour)")
    _report(cost_coarse, "coarse", curated, args.base_cost)
    _report(cost_chebi, "+ChEBI", curated, args.base_cost)
    _report(cost_sib, f"+ChEBI+sib{args.sibling_weight:g}", curated, args.base_cost)


if __name__ == "__main__":
    main()
