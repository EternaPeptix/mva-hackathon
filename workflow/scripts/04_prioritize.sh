#!/usr/bin/env bash
# Genome-wide Exomiser ranking using explicit input paths and database versions.
# LIRICAL is a separate, still-unimplemented validation stage.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
exec "${MVA_PYTHON:-python3}" "$HERE/prioritize.py" "$@"
