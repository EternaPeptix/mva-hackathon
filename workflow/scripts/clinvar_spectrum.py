#!/usr/bin/env python3
"""
Measure the classification asymmetry in a gene: how often does each variant class
actually reach a definitive call?

Motivation. Our proband's second BUB1B allele classifies as a VUS under strict ACMG/AMP.
The interesting question is not whether we got that call right, but whether it is the
*normal* outcome for a missense allele in this gene. If it is, then the VUS is a property
of the gene's evidence landscape rather than a property of our variant — and that is a
statement about rare-disease genomics that can be measured rather than asserted.

Method. Read every ClinVar record over the gene span from a locally pinned ClinVar VCF,
assign each to a molecular class from the MC (molecular consequence) field, cross-tabulate
against CLNSIG, and report the rate at which each class reaches P/LP.

Deliberate choices, so the numbers can be audited:
  - Records are kept only if GENEINFO names the gene, so neighbouring genes do not leak in.
  - "Truncating" means the record carries ANY of nonsense / frameshift / splice donor /
    splice acceptor / initiator codon. A record with mixed consequences across transcripts
    is counted as truncating if any transcript truncates, which is the conservative
    direction for the claim being made (it inflates the truncating denominator).
  - CLNSIG values containing "Conflicting" are held out in their own bucket rather than
    folded into either side.
  - Review status is reported separately, because a pile of unreviewed submissions would
    be a much weaker claim than a pile of independently reviewed ones.

Usage:
    clinvar_spectrum.py --vcf ref/clinvar.vcf.gz --region 15:40150000-40230000 --gene BUB1B
"""
from __future__ import annotations
import argparse, collections, json, subprocess, sys

TRUNCATING = {
    "nonsense", "frameshift_variant", "splice_donor_variant",
    "splice_acceptor_variant", "initiator_codon_variant",
}
DEFINITIVE_P = {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}


def consequences(mc: str) -> set[str]:
    return {m.split("|")[1] for m in mc.split(",") if "|" in m}


def variant_class(mc: str) -> str:
    so = consequences(mc)
    if so & TRUNCATING:
        return "truncating"
    for name in ("missense_variant", "synonymous_variant", "intron_variant"):
        if name in so:
            return {"missense_variant": "missense",
                    "synonymous_variant": "synonymous",
                    "intron_variant": "intronic"}[name]
    return "other"


def significance(clnsig: str) -> str:
    s = clnsig.replace("_", " ").lower().strip()
    if "conflicting" in s:
        return "conflicting"
    if s in DEFINITIVE_P:
        return "P/LP"
    if "uncertain" in s:
        return "VUS"
    if "benign" in s:
        return "B/LB"
    return "other"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--region", required=True, help="chrom:start-end over the gene span")
    ap.add_argument("--gene", required=True, help="symbol as it appears in GENEINFO")
    ap.add_argument("--json", help="write machine-readable results here")
    a = ap.parse_args()

    proc = subprocess.run(
        ["bcftools", "query", "-f",
         "%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/MC\t%INFO/CLNREVSTAT\t%INFO/GENEINFO\n",
         "-r", a.region, a.vcf],
        capture_output=True, text=True)
    if proc.returncode != 0:
        print(proc.stderr, file=sys.stderr)
        return 1

    rows = [l.split("\t") for l in proc.stdout.strip().split("\n") if l.strip()]
    rows = [r for r in rows if a.gene in r[6]]
    if not rows:
        print(f"no {a.gene} records in {a.region}", file=sys.stderr)
        return 1

    table = collections.Counter((significance(r[3]), variant_class(r[4])) for r in rows)
    classes = ["truncating", "missense", "synonymous", "intronic", "other"]
    sigs = ["P/LP", "VUS", "B/LB", "conflicting", "other"]

    print(f"{a.gene}: {len(rows)} ClinVar records over {a.region}\n")
    head = f"{'':14s}" + "".join(f"{c:>13s}" for c in classes)
    print(head); print("-" * len(head))
    for s in sigs:
        cells = [table[(s, c)] for c in classes]
        if any(cells):
            print(f"{s:14s}" + "".join(f"{v:13d}" for v in cells))

    print("\nRate at which each class reaches a definitive P/LP call:")
    rates = {}
    for c in classes:
        # denominator = records that got any curated call at all
        denom = sum(table[(s, c)] for s in ("P/LP", "VUS", "B/LB"))
        if denom:
            rates[c] = table[("P/LP", c)] / denom
            print(f"  {c:12s} {table[('P/LP', c)]:5d}/{denom:<6d} = {100 * rates[c]:6.2f}%")

    if rates.get("missense", 0) and rates.get("truncating", 0):
        print(f"\n  truncating alleles are {rates['truncating'] / rates['missense']:.0f}x "
              f"more likely to reach P/LP than missense alleles")

    # A pile of unreviewed submissions would be a much weaker claim. Check.
    vus_mis = [r for r in rows
               if variant_class(r[4]) == "missense" and significance(r[3]) == "VUS"]
    if vus_mis:
        rev = collections.Counter(r[5] for r in vus_mis)
        print(f"\nReview status of the {len(vus_mis)} missense VUS:")
        for k, v in rev.most_common():
            print(f"  {k.replace('_', ' '):46s} {v:5d}  ({100 * v / len(vus_mis):5.1f}%)")
        multi = sum(v for k, v in rev.items() if "multiple" in k)
        print(f"\n  independently reviewed by >1 submitter and STILL uncertain: "
              f"{multi} ({100 * multi / len(vus_mis):.1f}%)")

    if a.json:
        with open(a.json, "w") as fh:
            json.dump({
                "gene": a.gene, "region": a.region, "records": len(rows),
                "table": {f"{s}|{c}": table[(s, c)] for s in sigs for c in classes
                          if table[(s, c)]},
                "plp_rate": rates,
            }, fh, indent=2)
        print(f"\nwrote {a.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
