#!/usr/bin/env bash
# 03_annotate.sh — VEP annotation with the fields ACMG classification actually needs.
#
# Run offline against a local cache rather than the REST API. Two reasons that matter for a
# clinical-grade result: (1) the REST API silently caps and rate-limits, so a genome-wide run
# is not reproducible; (2) a cache pins the annotation to a stated Ensembl release, so the
# result can be regenerated years later. Every consequence we report is versioned.
#
# Patient data never leaves this machine. No hosted API sees a single variant.
set -euo pipefail
usage() { echo "usage: $0 -i IN.vcf.gz -o OUT.vcf.gz -r REF.fa -c VEP_CACHE_DIR [-t threads]" >&2; exit 2; }
THREADS=$(sysctl -n hw.logicalcpu 2>/dev/null || nproc)
while getopts "i:o:r:c:t:" o; do case $o in
  i) IN=$OPTARG;; o) OUT=$OPTARG;; r) REF=$OPTARG;; c) CACHE=$OPTARG;; t) THREADS=$OPTARG;; *) usage;;
esac; done
: "${IN:?}" "${OUT:?}" "${REF:?}" "${CACHE:?}"
VEP_IMAGE=${VEP_IMAGE:-ensemblorg/ensembl-vep:release_115.2}

# --pick_allele_gene with a fixed --pick_order gives one deterministic consequence per
# allele/gene. Leaving the default ordering in place would let VEP's tie-breaks drift
# between releases; naming the order explicitly makes the choice auditable.
docker run --rm \
  -v "$(realpath "$CACHE")":/cache \
  -v "$(dirname "$(realpath "$IN")")":/in \
  -v "$(dirname "$(realpath "$OUT")")":/out \
  -v "$(dirname "$(realpath "$REF")")":/ref \
  "$VEP_IMAGE" vep \
    --input_file /in/"$(basename "$IN")" \
    --output_file /out/"$(basename "$OUT")" \
    --vcf --compress_output bgzip --force_overwrite \
    --cache --offline --dir_cache /cache --cache_version 115 \
    --assembly GRCh38 --fasta /ref/"$(basename "$REF")" \
    --species homo_sapiens --fork "$THREADS" \
    --hgvs --hgvsg --symbol --canonical --mane --biotype --domains --numbers \
    --protein --uniprot --ccds --transcript_version \
    --sift b --polyphen b \
    --af --af_gnomade --af_gnomadg --max_af \
    --pick_allele_gene \
    --pick_order mane_select,mane_plus_clinical,canonical,biotype,rank,length \
    --stats_file /out/"$(basename "${OUT%.vcf.gz}")"_vep_summary.html

bcftools index -t -f "$OUT"
echo "[done] $OUT"
echo "NOTE: SIFT/PolyPhen come from the cache. REVEL, AlphaMissense, SpliceAI and CADD are"
echo "      NOT in the base cache; Exomiser (stage 04) supplies REVEL/AlphaMissense/SpliceAI"
echo "      from its own 2512 bundle, and those are the values used for classification."
