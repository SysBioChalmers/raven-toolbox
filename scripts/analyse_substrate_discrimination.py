#!/usr/bin/env python3
"""Does the ChEBI layer tell apart cargo that the coarse class collapses? (model-free, intrinsic.)

The yeast-GEM validation showed the ChEBI layer barely moves the *cost* metric there, because yeast's
dominant transporter families (MFS/MCF/ABC) are promiscuous, so a coarse class already fits. This script
isolates the layer's actual contribution — substrate *specificity* — without any organism/localisation
confound, using only the hosted artefacts (`tcdb_substrates.tsv` + `chebi_relations.tsv.gz`).

Within each coarse TCDB family (:data:`TC_FAMILY_CLASS`), different transport systems carry different
specific substrates. For every (curated substrate `s`, transport system `T`) pair *in the same family*
we score how strongly the substrate layer matches:

* **coarse class** — every system in the family shares one class, so it matches `s` to *all* of them
  (score 1.0); it cannot say which member actually carries `s`.
* **ChEBI roll-up** — :meth:`SubstrateOntology.match` scores `s` against each system's *own* curated
  substrates: ~1.0 for a true carrier, lower for a family sibling with different cargo.

So the discrimination the coarse class cannot do = the fraction of same-family *non*-carriers the ChEBI
layer scores below a threshold. High = the layer separates cargo a coarse class collapses — the payoff
expected on organisms with narrow-specificity transporters. ASCII-only.
Usage: ``python scripts/analyse_substrate_discrimination.py``
"""
from __future__ import annotations

import argparse
from collections import defaultdict

from raven_toolbox.localization.substrate_ontology import SubstrateOntology, load_tc_substrates
from raven_toolbox.localization.transporter_tables import TC_FAMILY_CLASS


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rule-out", type=float, default=0.5,
                    help="a same-family non-carrier is 'ruled out' when its ChEBI score is below this")
    args = ap.parse_args(argv)

    onto = SubstrateOntology.load()
    tc_sub = load_tc_substrates()

    # group TCDB systems by coarse family (TC-ID's first three levels), keeping curated families only
    by_family: dict[str, dict[str, frozenset[str]]] = defaultdict(dict)
    for tc_id, subs in tc_sub.items():
        fam = ".".join(tc_id.split(".")[:3])
        if fam in TC_FAMILY_CLASS:
            by_family[fam][tc_id] = subs

    # "related" = a same-family non-carrier whose cargo is chemically close enough that the roll-up
    # fires at all (score > 0) -- the *hard* decoys (e.g. glucose- vs xylose-porter, both "sugar").
    print(f"{'family':8s} {'coarse class':26s} {'sys':>4} {'ruled-out':>9} {'related':>8} {'rel~':>6}")
    tot_true: list[float] = []
    tot_false: list[float] = []
    for fam in sorted(by_family):
        systems = by_family[fam]
        if len(systems) < 2:
            continue  # need siblings to discriminate against
        subs_here = sorted({s for subs in systems.values() for s in subs})
        f_true: list[float] = []
        f_false: list[float] = []
        for s in subs_here:
            for _tc, subs in systems.items():
                score = onto.match([s], subs)
                (f_true if s in subs else f_false).append(score)
        if not f_false:
            continue
        related = [x for x in f_false if x > 0.0]
        ruled = sum(1 for x in f_false if x < args.rule_out) / len(f_false)
        cls = ",".join(sorted(TC_FAMILY_CLASS[fam]))[:26]
        print(f"{fam:8s} {cls:26s} {len(systems):4d} {ruled:8.0%} "
              f"{len(related)/len(f_false):7.0%} {sum(related)/len(related) if related else 0.0:6.2f}")
        tot_true += f_true
        tot_false += f_false

    related = [x for x in tot_false if x > 0.0]
    ruled = sum(1 for x in tot_false if x < args.rule_out) / max(1, len(tot_false))
    print("-" * 66)
    print(f"{'OVERALL':8s} {'ChEBI roll-up':26s} {'':4s} {ruled:8.0%} "
          f"{len(related)/max(1,len(tot_false)):7.0%} {sum(related)/max(1,len(related)):6.2f}")
    print(f"{'':8s} {'coarse class alone':26s} {'':4s} {0.0:8.0%} {'100%':>7s} {1.0:6.2f}"
          "  (same class -> no discrimination)")
    print(f"\n{len(tot_true)} true (substrate, carrier) pairs [ChEBI mean "
          f"{sum(tot_true)/max(1,len(tot_true)):.2f}], {len(tot_false)} same-family non-carrier pairs. "
          f"ChEBI rules out {ruled:.0%} of the non-carriers the coarse class scores identically; the "
          f"{len(related)} of them with chemically-related cargo still score only "
          f"{sum(related)/max(1,len(related)):.2f} on average (vs 1.00 for a true carrier), so even "
          "those stay separable.")


if __name__ == "__main__":
    main()
