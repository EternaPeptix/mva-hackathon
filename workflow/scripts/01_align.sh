#!/usr/bin/env bash
# 01_align.sh — FASTQ -> analysis-ready CRAM (GRCh38 no-alt + decoys), lane-resumable.
#
# Rationale: the challenge ships both a vendor VCF and the raw FASTQs. We re-align and
# re-call from FASTQ rather than trusting the vendor VCF, so that every reported variant
# is supported by evidence we generated ourselves and can inspect at read level.
#
# The vendor VCF header declares:
#   GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_hs38d1_maskedGRC_exclusions_v2_no_chr
# i.e. GRCh38 no-alt + decoys, contig names lacking the "chr" prefix. We align to the SAME
# reference so our calls are directly comparable to the vendor's, and normalise naming only
# at the reporting step.
#
# RESUMABILITY. This script is safe to run repeatedly and safe to run while the FASTQs are
# still downloading. It aligns any lane that (a) has both mates present and complete, and
# (b) has no finished BAM yet; it skips the rest and reports them. Merge/markdup/QC only run
# once every lane advertised in the run manifest is aligned. That means CPU can start work on
# early lanes while later lanes are still arriving, without ever producing a partial CRAM that
# looks complete. A lane BAM is written to a .tmp name and renamed only on success, so an
# interrupted run can never leave a truncated BAM that a later run would mistake for done.
set -euo pipefail
ORIGINAL_ARGS=("$@")
HERE=$(cd "$(dirname "$0")" && pwd)

usage() {
  cat >&2 <<'EOF'
usage: 01_align.sh -r REF.fa -o OUTDIR -s SAMPLE -f FASTQ_DIR [-t THREADS] [-n NLANES]

  -n NLANES  number of lanes the finished dataset should have (default: 4, the challenge
             dataset). A partial download must never silently redefine the expected total.
EOF
  exit 2
}

THREADS=8
NLANES=4
while getopts "r:o:s:f:t:n:" o; do case $o in
  r) REF=$OPTARG;; o) OUT=$OPTARG;; s) SAMPLE=$OPTARG;; f) FQDIR=$OPTARG;;
  t) THREADS=$OPTARG;; n) NLANES=$OPTARG;; *) usage;;
esac; done
: "${REF:?}" "${OUT:?}" "${SAMPLE:?}" "${FQDIR:?}"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] && [[ "$NLANES" =~ ^[1-9][0-9]*$ ]] || usage
mkdir -p "$OUT/logs"
if [ "${MVA_ALIGN_LOCK_HELD:-}" != 1 ]; then
  exec "${MVA_PYTHON:-python3}" "$HERE/with_lock.py" "$OUT/.alignment.lock" \
    env MVA_ALIGN_LOCK_HELD=1 bash "$0" "${ORIGINAL_ARGS[@]}"
fi
BWA=${BWA:-bwa}
command -v "$BWA"   >/dev/null || { echo "bwa not found"      >&2; exit 1; }
command -v samtools >/dev/null || { echo "samtools not found" >&2; exit 1; }
[ -e "${REF}.bwt" ] || { echo "bwa index missing for $REF — run: bwa index $REF" >&2; exit 1; }

# Sort memory is per-thread. samtools sort -@N -m M can reserve N*M, so cap it rather than
# letting a 16-core box try to reserve 32 GB in the middle of a markdup stream.
SORT_THREADS=$(( THREADS > 2 ? 2 : 1 ))
SORT_MEM=${SORT_MEM:-2G}

# A gzip member that is still being written fails `gzip -t`. That is exactly the test we want
# before feeding a file to bwa: it distinguishes "downloaded" from "downloading".
complete() { [ -s "$1" ] && gzip -t "$1" 2>/dev/null; }

pending=(); done_lanes=()
for R1 in "$FQDIR"/*_R1_001.fastq.gz; do
  [ -e "$R1" ] || { echo "no FASTQs in $FQDIR" >&2; exit 1; }
  R2=${R1/_R1_001/_R2_001}
  LANE=$(basename "$R1" | sed -E 's/.*_(L[0-9]{3})_R1_001\.fastq\.gz/\1/')
  BAM="$OUT/${SAMPLE}.${LANE}.bam"

  if [ -e "$BAM" ]; then
    samtools quickcheck -v "$BAM"
    if [ ! -s "${BAM}.bai" ]; then
      samtools index -@1 "$BAM" "${BAM}.bai.new"
      mv "${BAM}.bai.new" "${BAM}.bai"
    fi
    echo "[skip] $LANE already aligned"; done_lanes+=("$BAM"); continue
  fi
  if ! complete "$R1" || ! complete "$R2"; then
    echo "[wait] $LANE incomplete or missing mate — deferring"; pending+=("$LANE"); continue
  fi

  # PU (platform unit) = flowcell.lane, the field GATK BQSR actually groups on.
  #
  # Derive the flowcell and lane from the READ HEADER, not the filename. The filename's second
  # underscore-field is the sample ID (EX2312012), not the flowcell — taking it from there
  # silently mislabels every read group, and BQSR would then be grouping on a constant. The
  # Illumina header fields are instrument:run:flowcell:lane:tile:x:y.
  # The acquisition-specific example is omitted from this release copy.
  # Use flowcell.lane as the platform unit. A sample barcode can carry
  # sequencing error and is not used to define the platform unit.
  # Read the first line WITHOUT `| head -1`: head exits after one line, gzip dies of SIGPIPE,
  # and under `set -o pipefail` that surfaces as exit 141 and kills the run. Process
  # substitution plus a single `read` gets the same line without a checked pipeline.
  IFS= read -r HDR < <(gzip -cd "$R1" 2>/dev/null) || true
  FLOWCELL=$(printf '%s' "$HDR" | awk -F: '{print $3}')
  LANENUM=$( printf '%s' "$HDR" | awk -F: '{print $4}')
  [ -n "$FLOWCELL" ] && [ -n "$LANENUM" ] || { echo "cannot parse header of $R1" >&2; exit 1; }
  # Cross-check the header against the filename; disagreement means the file is mislabelled.
  [ "L00${LANENUM}" = "$LANE" ] || { echo "lane mismatch: file says $LANE, reads say L00${LANENUM}" >&2; exit 1; }
  RG="@RG\tID:${FLOWCELL}.${LANENUM}\tSM:${SAMPLE}\tLB:${SAMPLE}_lib1\tPL:ILLUMINA\tPU:${FLOWCELL}.${LANENUM}"
  echo "[align] $LANE  RG=${FLOWCELL}.${LANENUM}  threads=$THREADS"

  # bwa, not bwa-mem2: bwa-mem2 requires x86 AVX2 and will not run on Apple Silicon at all.
  # Plain bwa is the portable original, gives equivalent alignments, and is what the GATK Best
  # Practices are calibrated against.
  #   -Y  soft-clip supplementary alignments, keeping full read sequence visible for the
  #       SV/breakpoint checks downstream.
  #   -K 100000000  fix the batch size so alignment is deterministic regardless of thread
  #       count. Without it bwa's output depends on thread scheduling and the run is not
  #       reproducible — a silent reproducibility hole that is easy to miss.
  "$BWA" mem -t "$THREADS" -K 100000000 -Y -R "$RG" "$REF" "$R1" "$R2" \
      2> "$OUT/logs/bwa.${LANE}.log" \
    | samtools sort -@1 -m "$SORT_MEM" -o "${BAM}.tmp" -
  samtools quickcheck -v "${BAM}.tmp"
  samtools index -@1 "${BAM}.tmp"
  mv "${BAM}.tmp" "$BAM"; mv "${BAM}.tmp.bai" "${BAM}.bai"
  done_lanes+=("$BAM")
done

EXPECTED=$NLANES
echo "[status] ${#done_lanes[@]} lane BAM(s) ready; expecting $EXPECTED"

if [ "${#done_lanes[@]}" -lt "$EXPECTED" ]; then
  echo "[hold] deferring merge/markdup until all $EXPECTED lanes are aligned."
  [ "${#pending[@]}" -gt 0 ] && echo "[hold] pending lanes: ${pending[*]}"
  echo "[hold] re-run this script when the remaining FASTQs finish."
  exit 0
fi

if [ "${#done_lanes[@]}" -ne "$EXPECTED" ]; then
  echo "[FAIL] unexpected additional lanes; refusing finalisation" >&2; exit 1
fi
FINAL="$OUT/${SAMPLE}.cram"
FINGERPRINT=$("${MVA_PYTHON:-python3}" "$HERE/input_fingerprint.py" "$REF" "${done_lanes[@]}")
if [ -e "$FINAL" ]; then
  if [ -s "${FINAL}.crai" ] && [ -s "$OUT/logs/finalization.complete" ] && \
      [ "$(cat "$OUT/logs/finalization.complete")" = "$EXPECTED $FINGERPRINT" ]; then
    samtools quickcheck -v "$FINAL"
    echo "[done] existing verified CRAM retained"; exit 0
  fi
  echo "[FAIL] existing CRAM has no matching completion record; preserved for inspection" >&2
  exit 1
fi

echo "[merge] ${#done_lanes[@]} lanes"
MERGED="$OUT/${SAMPLE}.merged.$(date -u +%Y%m%dT%H%M%SZ).$$.bam"
samtools merge -@ "$SORT_THREADS" -o "$MERGED" "${done_lanes[@]}"
samtools quickcheck -v "$MERGED"
samtools index -@1 "$MERGED"

# Duplicate marking must see all lanes at once: a PCR duplicate pair can be split across lanes,
# and marking lane-by-lane would miss exactly those.
echo "[markdup]"
samtools collate -@1 -O -u "$MERGED" \
  | samtools fixmate -@1 -m -u - - \
  | samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM" -u - \
  | samtools markdup -@ "$SORT_THREADS" -f "$OUT/logs/markdup.stats" \
      --reference "$REF" -O CRAM - "${FINAL}.tmp"
samtools quickcheck -v "${FINAL}.tmp"
samtools index "${FINAL}.tmp" "${FINAL}.tmp.crai"

echo "[qc]"
samtools flagstat -@1 --input-fmt-option "reference=$REF" "${FINAL}.tmp" > "$OUT/logs/flagstat.txt"
samtools stats -@1 --input-fmt-option "reference=$REF" -r "$REF" "${FINAL}.tmp" > "$OUT/logs/samtools_stats.txt"
samtools coverage --input-fmt-option "reference=$REF" "${FINAL}.tmp" > "$OUT/logs/coverage_by_contig.txt"

# Keep the earlier lane BAMs and merge for independent validation and recovery.
mv "${FINAL}.tmp" "$FINAL"
mv "${FINAL}.tmp.crai" "${FINAL}.crai"
printf '%s %s\n' "$EXPECTED" "$FINGERPRINT" > "$OUT/logs/finalization.complete"
echo "[done] $FINAL; lane BAMs and merge retained"
