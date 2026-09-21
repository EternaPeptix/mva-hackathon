#!/bin/sh
# Replay existing computation; no network or patient genomic files are needed.
set -eu
if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
  echo "Usage: bash workflow/track2-model/run.sh NEW_OUTPUT_DIRECTORY [PINNED_GRCH38_FASTA]" >&2
  exit 2
fi
model_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
model_python=${MVA_MODEL_PYTHON:-python3}
command -v "$model_python" >/dev/null
mkdir -- "$1"
model_output=$(CDPATH= cd -- "$1" && pwd)
if [ "$#" -eq 2 ]; then
  "$model_python" "$model_dir/alleles.py" --source-directory "$model_dir/public-transcript" \
    --reference "$2" --output "$model_output/allele-inputs.json"
else
  cp "$model_dir/allele-inputs.json" "$model_output/allele-inputs.json"
fi
"$model_python" "$model_dir/population.py" --alleles "$model_output/allele-inputs.json" \
  --output "$model_output/population-results.json"
echo "Completed assumption sensitivity. Results are not predicted drug responses: $model_output"
