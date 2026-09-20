#!/usr/bin/env bash
# 00_prepare_reference.sh — reconstruct the EXACT reference frame the vendor VCF was called against.
#
# The vendor VCF header declares:
#   ##reference=file://refGenome/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr.fasta
#   ##contig=<ID=15,length=101991189>
#   ##contig=<ID=15_KI270727v1_random,length=448248>
#
# Decoded, that is: the NCBI GRCh38 **no-alt analysis set with UCSC-style names, plus the
# hs38d1 decoy**, with the leading "chr" stripped from every contig name.
#
# This is NOT interchangeable with the Ensembl primary assembly. Ensembl names the same
# scaffold "KI270727.1" where the analysis set calls it "chr15_KI270727v1_random", and
# Ensembl uses "MT" where the analysis set uses "chrM". Aligning to Ensembl and then
# comparing to the vendor VCF would silently mismatch on every non-primary contig. We
# therefore rebuild the vendor's exact reference so our calls are directly comparable and
# any disagreement is biological rather than an artefact of reference choice.
#
# Note on decoys: hs38d1 exists to absorb reads from sequence absent in GRCh38. Omitting it
# pushes those reads onto real chromosomes as low-MAPQ noise. Including it is what the
# vendor did, so we do it too.
set -euo pipefail
OUT=${1:?usage: $0 OUTDIR}
mkdir -p "$OUT"; cd "$OUT"

URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
RAW=GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
FINAL=GRCh38_no_alt_plus_hs38d1_no_chr.fa

if [ ! -e "$FINAL" ]; then
  [ -e "$RAW" ] || { echo "[dl] analysis set (~900 MB)"; curl -fL --retry 5 --retry-all-errors -o "$RAW" "$URL"; }
  echo "[strip] removing 'chr' prefix to match the vendor's contig naming"
  # Only the leading "chr" of a header token is stripped: chr15 -> 15,
  # chr15_KI270727v1_random -> 15_KI270727v1_random, chrM -> M, chrUn_* -> Un_*.
  # Decoy contigs (chrUn_JTFH*) follow the same rule. Sequence lines are untouched.
  gzcat "$RAW" | sed -E '/^>/ s/^>chr/>/' > "$FINAL"
fi

echo "[index] samtools faidx"
samtools faidx "$FINAL"
echo "[index] GATK sequence dictionary"
${GATK:-gatk} CreateSequenceDictionary -R "$FINAL" 2>/dev/null || true
echo "[index] bwa index — slowest step, ~1 h for 3.1 Gb"
[ -e "${FINAL}.bwt" ] || bwa index "$FINAL"

echo "[verify] contig names/lengths vs the vendor VCF"
echo "  expected: 15 = 101991189 ; 15_KI270727v1_random = 448248"
awk '$1=="15"||$1=="15_KI270727v1_random"{printf "  actual:   %s = %s\n",$1,$2}' "${FINAL}.fai"
echo "[done] $OUT/$FINAL"
