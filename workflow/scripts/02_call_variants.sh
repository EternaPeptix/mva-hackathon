#!/usr/bin/env bash
# 02_call_variants.sh — independent variant calling from our own alignment.
#
# We call with THREE independent methods and report the intersection for the reported
# alleles. This matters because the two candidate alleles carry the entire diagnosis:
# a single-caller artefact would invalidate the result, and the three callers here fail
# in different ways (local reassembly vs. pileup vs. learned model), so agreement is
# meaningful rather than merely repeated.
#
#   1. GATK HaplotypeCaller  — local de-novo reassembly (matches the vendor's method,
#                              so disagreement with the vendor VCF is interpretable)
#   2. bcftools mpileup/call — direct pileup + multiallelic caller, no reassembly
#   3. DeepVariant           — CNN over pileup images, orthogonal error model
#
# DeepVariant publishes no arm64 image, so on Apple Silicon it runs under emulation. That
# is impractical genome-wide but perfectly fast over a locus, which is why --region exists.
set -euo pipefail

usage() { echo "usage: $0 -r REF.fa -b SAMPLE.cram -o OUTDIR -s SAMPLE [-L region] [-t threads]" >&2; exit 2; }
THREADS=$(sysctl -n hw.logicalcpu 2>/dev/null || nproc); REGION=""
while getopts "r:b:o:s:L:t:" o; do case $o in
  r) REF=$OPTARG;; b) CRAM=$OPTARG;; o) OUT=$OPTARG;; s) SAMPLE=$OPTARG;;
  L) REGION=$OPTARG;; t) THREADS=$OPTARG;; *) usage;;
esac; done
: "${REF:?}" "${CRAM:?}" "${OUT:?}" "${SAMPLE:?}"
mkdir -p "$OUT/logs"
GATK=${GATK:-gatk}
RARG=(); [ -n "$REGION" ] && RARG=(-L "$REGION")
BARG=(); [ -n "$REGION" ] && BARG=(-r "$REGION")

echo "[1/3] GATK HaplotypeCaller"
$GATK --java-options "-Xmx12g" HaplotypeCaller \
  -R "$REF" -I "$CRAM" -O "$OUT/${SAMPLE}.gatk.vcf.gz" \
  "${RARG[@]}" --native-pair-hmm-threads "$THREADS" \
  2> "$OUT/logs/gatk_hc.log"

# Hard filters below are the GATK-recommended thresholds for a single sample with no
# truth set available for VQSR. They are the same filters the vendor applied (their header
# records VariantFiltration with QD/MQ/FS), so our FILTER column is comparable to theirs.
$GATK VariantFiltration -R "$REF" -V "$OUT/${SAMPLE}.gatk.vcf.gz" \
  -O "$OUT/${SAMPLE}.gatk.filt.vcf.gz" \
  --filter-name QD2   --filter-expression "QD < 2.0" \
  --filter-name MQ40  --filter-expression "MQ < 40.0" \
  --filter-name FS60  --filter-expression "FS > 60.0" \
  --filter-name SOR3  --filter-expression "SOR > 3.0" \
  --filter-name MQRS  --filter-expression "MQRankSum < -12.5" \
  --filter-name RPRS  --filter-expression "ReadPosRankSum < -8.0" \
  2> "$OUT/logs/gatk_filt.log"

echo "[2/3] bcftools mpileup/call"
# -a AD,DP,SP keeps per-allele depths so we can audit allele balance at the reported sites.
# Baseline BAQ is left ON: it is the main defence against false SNPs beside indels, and both
# candidate alleles sit in ordinary sequence where BAQ costs us nothing.
bcftools mpileup -Ou -f "$REF" "${BARG[@]}" -a AD,DP,SP -q 20 -Q 20 --threads "$THREADS" "$CRAM" \
  | bcftools call -mv -Ou --threads "$THREADS" \
  | bcftools norm -f "$REF" -m -any -Oz -o "$OUT/${SAMPLE}.bcftools.vcf.gz" --threads "$THREADS"
bcftools index -t "$OUT/${SAMPLE}.bcftools.vcf.gz"

echo "[3/3] DeepVariant"
if [ -z "$REGION" ]; then
  echo "  ! no --region given; skipping DeepVariant (emulated genome-wide run is impractical)"
else
  docker run --rm -v "$(dirname "$(realpath "$REF")")":/ref -v "$(realpath "$OUT")":/out \
    -v "$(dirname "$(realpath "$CRAM")")":/in \
    google/deepvariant:1.6.1 /opt/deepvariant/bin/run_deepvariant \
      --model_type=WGS --ref=/ref/"$(basename "$REF")" --reads=/in/"$(basename "$CRAM")" \
      --regions="$REGION" --output_vcf=/out/"${SAMPLE}".dv.vcf.gz \
      --num_shards="$THREADS" > "$OUT/logs/deepvariant.log" 2>&1 \
    || echo "  ! DeepVariant failed — see $OUT/logs/deepvariant.log"
fi

echo "[concordance]"
# Report each caller's genotype at every site any caller emitted, so the report can state
# concordance rather than assert it.
for f in "$OUT/${SAMPLE}".gatk.filt.vcf.gz "$OUT/${SAMPLE}".bcftools.vcf.gz "$OUT/${SAMPLE}".dv.vcf.gz; do
  [ -e "$f" ] || continue
  echo "--- $(basename "$f")"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT:%DP]\n' "$f" "${BARG[@]}"
done | tee "$OUT/concordance.txt"
echo "[done] $OUT"
