#!/usr/bin/env bash
# Resume a single resource with an OS-level lock and publisher/integrity verification.
# Usage: fetch_resource.sh URL DESTDIR [SHA1_OR_SHA256] [EXPECTED_BYTES]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
exec "${MVA_PYTHON:-python3}" "$HERE/fetch_resource.py" "$@"
