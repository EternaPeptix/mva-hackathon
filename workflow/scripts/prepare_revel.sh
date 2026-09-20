#!/usr/bin/env bash
# prepare_revel.sh — build a GRCh38-indexed REVEL file for the VEP REVEL plugin.
#
# REVEL v1.3 ships one row per possible missense change, sorted by *GRCh37* position, with the
# GRCh38 coordinate as a secondary column that is "." wherever the variant does not lift over.
# Two consequences, both of which silently break a naive `tabix -s 1 -b 3`:
#   1. rows with grch38_pos == "." are unparseable as positions;
#   2. the file is not in GRCh38 sort order, and tabix requires coordinate-sorted input.
# Both must be handled or the index build fails (or, worse, succeeds and returns wrong rows).
# This follows the procedure documented by Ensembl for the GRCh38 build of the plugin file.
set -euo pipefail

DIR=${1:-/Users/jeweled/mva-hackathon/ref/vep_plugins}
cd "$DIR"

SRC=revel_with_transcript_ids
OUT=revel_grch38.tsv.gz
[ -s "$SRC" ] || { echo "missing $SRC (unzip revel-v1.3_all_chromosomes.zip first)" >&2; exit 1; }

echo "[revel] filtering non-liftover rows, re-sorting on GRCh38 coordinate..."
{
  # Header must begin with '#' so tabix treats it as a comment line.
  head -1 "$SRC" | tr ',' '\t' | sed '1s/^/#/'
  tail -n +2 "$SRC" | tr ',' '\t' | awk -F'\t' '$3 != "." && $3 != ""' | sort -T . -k1,1 -k3,3n
} | bgzip -@ 6 -c > "$OUT"

echo "[revel] indexing..."
tabix -f -s 1 -b 3 -e 3 "$OUT"

# Prove the index actually resolves the variant we care about rather than trusting exit status.
echo "[revel] spot-check at BUB1B 15:40220612 (expect the N1002K row):"
tabix "$OUT" 15:40220612-40220612 || true
echo "[revel] done: $(ls -lh "$OUT" | awk '{print $5}')"
