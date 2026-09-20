#!/usr/bin/env bash
# fetch_fastqs.sh — fetch the challenge FASTQs, one at a time, and verify every byte.
#
# Why not `hf download`. We used it first and it cost us three restarts and ~18 GB of wasted
# transfer. Two failure modes, both silent:
#
#   1. It parallelises across files. Four 10 GB streams over one link do not go faster than one;
#      they go the same speed and fail together. Measured here: ~157 KB/s aggregate either way.
#   2. Its partial file carries a per-session random suffix. On restart it does not find the
#      previous partial, so it begins again at byte 0 and orphans what was already fetched. We
#      caught this with two generations of `.incomplete` blobs on disk for the same content hash
#      — 18.3 GB of correct bytes that the tool could no longer see.
#
# The verified downloader uses a deterministic .part name and an OS-managed lock. A serial
# queue preserves interrupted transfers and checks publisher hashes before accepting files.
#
# The digests below come from the dataset's own LFS metadata (the `lfs.oid` field of the repo
# tree API), so they are the publisher's hashes rather than ours. Publishing them here means a
# reviewer can confirm they hold byte-identical inputs WITHOUT anyone redistributing controlled
# patient data — a hash is not the data. Verification is also what makes resuming from a
# salvaged partial safe: if the salvaged prefix were wrong, the sha256 would reject the file.
set -euo pipefail

DEST=${1:-/Users/jeweled/mva-hackathon/data}
HERE=$(cd "$(dirname "$0")" && pwd)
REPO="datasets/SageBio/mva-hackathon-2026-data"


# name <TAB> sha256 <TAB> size — the full lane set, so this script is also the manifest.
MANIFEST=$(cat <<'EOF'
WGS_EX2312012_HGWCNDSX7_S16_L001_R1_001.fastq.gz	5ec244c0648552f2b23ae5b1a0b1350e8e464eb5cccd40178a109c44c174939f	10467073561
WGS_EX2312012_HGWCNDSX7_S16_L001_R2_001.fastq.gz	b80f2cc9f29c229a2bfe557988a4b93c74494fc0d9108b3752bc04cf57239cab	11035873766
WGS_EX2312012_HGWCNDSX7_S16_L002_R1_001.fastq.gz	0f60e090d9873a18a6914d374fa070a20502aa546d05b51ec520bd9e06ae7e71	10375801767
WGS_EX2312012_HGWCNDSX7_S16_L002_R2_001.fastq.gz	992ee50db8a89709b7adb327daabcacbf58261503f8c439bb18a3c682e43cd91	10911120034
WGS_EX2312012_HGWCNDSX7_S16_L003_R1_001.fastq.gz	bbeb51349d4411a661c5922e44650b357d968551ef85735ebb39e6c6071e5d46	10195738985
WGS_EX2312012_HGWCNDSX7_S16_L003_R2_001.fastq.gz	f3f08f04cd17b089ceb8a8598634f3eb28bc299dc1c00cbfa8e24ba2adb9017b	10621354344
WGS_EX2312012_HGWCNDSX7_S16_L004_R1_001.fastq.gz	ebc3a6853b9afb5b839ba1e4ff7966eb4f716fba5fee5a39c526d1468079f438	10306364082
WGS_EX2312012_HGWCNDSX7_S16_L004_R2_001.fastq.gz	b16eeae1f84c57fdb889349c25d625eb96a9cb265059004081630131e1037d8b	10755107565
EOF
)


while IFS=$'\t' read -r name sha size; do
  [ -n "$name" ] || continue
  echo "[queue] verifying or resuming $name"
  bash "$HERE/fetch_resource.sh" \
    "https://huggingface.co/$REPO/resolve/main/$name" "$DEST" "$sha" "$size" \
    || { echo "[FAIL] $name — files preserved for inspection or resumption" >&2; exit 1; }
done <<< "$MANIFEST"

echo "[done] all lanes present"
ls -l "$DEST"/*.fastq.gz | awk '{printf "  %.2f GB  %s\n", $5/1073741824, $9}'
