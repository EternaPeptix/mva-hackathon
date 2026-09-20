#!/usr/bin/env bash
# phase_audit.sh — inspect vendor PGT/PID annotations and their recorded block spans.
#
# Recorded blocks can supply local linkage evidence. Absent or different blocks do not
# exhaust the phasing information in the reads and cannot establish trans configuration.
# Dedicated read/mate phasing is a separate analysis in run_read_phasing.py.
#
# Output: the phasing status of the two candidate alleles plus the bridging het, and the span
# distribution of recorded PGT/PID blocks. The maximum span is descriptive of this callset;
# it is not a universal physical limit or proof that additional evidence cannot bridge a pair.
#
# Usage: phase_audit.sh -v <vcf.gz> [-a chrom:pos,chrom:pos,...] [-d <separation_bp>]

set -euo pipefail

VCF=""; SITES="15:40209701,15:40216470,15:40220612"; SEP=10911
while getopts "v:a:d:" o; do
  case "$o" in
    v) VCF=$OPTARG ;; a) SITES=$OPTARG ;; d) SEP=$OPTARG ;;
    *) echo "usage: $0 -v <vcf.gz> [-a sites] [-d separation_bp]" >&2; exit 2 ;;
  esac
done
: "${VCF:?-v <vcf.gz> is required}"
command -v bcftools >/dev/null || { echo "bcftools not found" >&2; exit 1; }

echo "== phasing status at the sites of interest =="
echo "   (PGT/PID '.' means no phase annotation is recorded; GT '/' means unphased)"
IFS=',' read -ra ARR <<< "$SITES"
for s in "${ARR[@]}"; do
  bcftools query -f '   %CHROM:%POS %REF>%ALT  GT=[%GT]  PGT=[%PGT]  PID=[%PID]  DP=[%DP]  GQ=[%GQ]\n' \
    -r "$s-${s##*:}" "$VCF"
done

echo
echo "== recorded vendor phase-block spans (every PASS variant, genome-wide) =="
bcftools query -f '%CHROM\t%POS\t[%PID]\n' -i 'FILTER="PASS"' "$VCF" \
  | awk -F'\t' '$3!="." && $3!=""' \
  | SEP="$SEP" python3 -c '
import sys, os, collections, statistics
blocks = collections.defaultdict(list)
n = 0
for line in sys.stdin:
    c, p, pid = line.rstrip("\n").split("\t")
    blocks[(c, pid)].append(int(p)); n += 1
if not blocks:
    print("   no phased variants found — does this VCF carry PGT/PID?"); sys.exit()
spans = sorted(max(v) - min(v) for v in blocks.values())
sep = int(os.environ["SEP"])
print(f"   phased PASS variants            : {n:,}")
print(f"   distinct phase blocks           : {len(blocks):,}")
print(f"   block span, median              : {statistics.median(spans):,.0f} bp")
print(f"   block span, 90th percentile     : {spans[int(0.9 * len(spans))]:,} bp")
print(f"   block span, MAXIMUM genome-wide : {max(spans):,} bp")
print(f"   blocks reaching >= {sep:,} bp   : {sum(1 for s in spans if s >= sep):,}")
print()
if max(spans) < sep:
    print(f"   VERDICT: the separation to bridge ({sep:,} bp) is {sep / max(spans):.0f}x the longest")
    print( "   phase block achieved anywhere in this genome. These vendor phase blocks do not establish phase across")
    print( "   the candidate interval. Additional evidence could include parental genotypes, long reads, or")
    print( "   allele-specific long-range PCR.")
else:
    print("   VERDICT: blocks of the required span exist — inspect the locus directly.")
'
