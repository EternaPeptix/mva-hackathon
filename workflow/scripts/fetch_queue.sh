#!/usr/bin/env bash
# Fetch the challenge FASTQs and large reference resources serially, with verification.
# The Studio ranking controller can prepare Exomiser first without waiting on FASTQs.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=${1:-$(cd "$HERE/../../.." && pwd)}
mkdir -p "$ROOT/.locks"
if [ "${MVA_FETCH_QUEUE_LOCK_HELD:-}" != 1 ]; then
  exec "${MVA_PYTHON:-python3}" "$HERE/with_lock.py" "$ROOT/.locks/fetch-queue.oslock" \
    env MVA_FETCH_QUEUE_LOCK_HELD=1 bash "$0" "$ROOT"
fi

echo "[queue] verify or resume every FASTQ against publisher hashes"
bash "$HERE/fetch_fastqs.sh" "$ROOT/data"

echo "[queue] Exomiser phenotype database"
bash "$HERE/fetch_resource.sh" \
  "https://data.monarchinitiative.org/exomiser/data/2512_phenotype.zip" "$ROOT/exomiser" "" 13132847504

echo "[queue] VEP 115 GRCh38 cache"
bash "$HERE/fetch_resource.sh" \
  "https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz" \
  "$ROOT/vep_cache"
echo "[done] queue verified successfully"
