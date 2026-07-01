#!/usr/bin/env python3
"""Validate evidence-aware transport scoring on a curated GEM (yeast-GEM by default).

The claim behind :mod:`raven_toolbox.localization.transport_evidence` is that weighting the transport
penalty by transporter evidence turns an *indiscriminate* cut into a *selective* one: real carriers get
cheap (retained), unsupported transports pay the full prior (dropped). This script measures that on a
curated model whose transport reactions are the ground truth.

It annotates the organism's proteome (``annotate_proteome`` -> hmmsearch + diamond), scores every
metabolite's transport cost (``evidence_aware_transport_cost`` with ``default_substrate_of`` and
DeepLoc-derived carrier compartments), then reports:

* the essential mitochondrial-shuttle carriers a blanket penalty wrongly drops -- their cost should be
  ~0 (fully evidenced, retained);
* the **selective-cut** metric: the curated-transport rate among *evidenced* (would-keep) metabolites
  vs among *unsupported* (would-drop) ones. Evidence is selective when keep-rate > drop-rate.

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


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True, type=Path, help="curated GEM (SBML), e.g. yeast-GEM.xml")
    ap.add_argument("--data-dir", type=Path, default=Path("data/deeploc"),
                    help="dir with *_proteins_*.fasta + *_deeploc_*.csv (default: data/deeploc)")
    ap.add_argument("--base-cost", type=float, default=0.5)
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

    # 2. carrier compartments from DeepLoc + evidence-aware cost (name-keyed metabolites)
    gene_comps = _gene_compartments(args.data_dir, set(ann))
    model = cobra.io.read_sbml_model(str(args.model))
    name_of = {m.name: m for m in model.metabolites}
    cost = evidence_aware_transport_cost(model, ann, gene_comps, substrate_of=default_substrate_of,
                                         base_cost=args.base_cost, base_metabolite=lambda m: m.name)

    # 3. essential carriers -- cost should be ~0 (evidenced, retained)
    print(f"\n=== essential carriers (base_cost {args.base_cost}; lower = more evidence, retained) ===")
    for name in ESSENTIAL_CARRIERS:
        if name in cost:
            cls = sorted(default_substrate_of(name_of[name]))
            flag = "retained" if cost[name] < args.base_cost else "NOT evidenced"
            print(f"  {name:22s} cost={cost[name]:.3f}  [{flag}]  classes={cls}")

    # 4. selective-cut metric: curated-transport rate among would-keep vs would-drop metabolites
    curated = {m.name for r in model.reactions if not r.boundary
               and len({mm.compartment for mm in r.metabolites}) > 1 for m in r.metabolites}
    allm = set(cost)
    curated &= allm
    keep = {b for b, c in cost.items() if c < args.base_cost}   # evidenced -> cheap -> kept
    drop = allm - keep                                          # unsupported -> full prior -> dropped
    keep_rate = len(keep & curated) / max(1, len(keep))
    drop_rate = len(drop & curated) / max(1, len(drop))
    base_rate = len(curated) / max(1, len(allm))
    print("\n=== selective cut (curated transports = ground truth) ===")
    print(f"metabolites={len(allm)}  curated-transported={len(curated)} ({base_rate:.0%} baseline)")
    print(f"  kept  (evidenced, cost<{args.base_cost}): {len(keep):4d}  curated-rate {keep_rate:.0%}")
    print(f"  dropped (unsupported, full cost):       {len(drop):4d}  curated-rate {drop_rate:.0%}")
    print(f"  --> selective: kept {keep_rate:.0%} vs dropped {drop_rate:.0%} "
          f"(blanket penalty would treat both at {base_rate:.0%})")
    print(f"  precision(kept are curated)={keep_rate:.0%}  "
          f"recall(curated are kept)={len(keep & curated)/max(1,len(curated)):.0%}")


if __name__ == "__main__":
    main()
