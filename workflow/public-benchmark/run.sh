#!/usr/bin/env bash
set -euo pipefail
# Only system Python 3.9+ is required; all bioinformatics dependencies are isolated.
# Do not inherit user JVM flags into the isolated replay runtime.
unset JAVA_TOOL_OPTIONS _JAVA_OPTIONS JDK_JAVA_OPTIONS JAVA_OPTS JAVA_OPTIONS JVM_OPTS
# Linux JDK ZIP memory mapping failed on the tested QEMU/virtiofs mount.
# This changes jar I/O only, not GATK/RTG algorithms or benchmark settings.
if [[ "$(uname -s)" == Linux ]]; then
  export JDK_JAVA_OPTIONS='-Dsun.zip.disableMemoryMapping=true'
fi
exec python3 "$(cd -- "$(dirname -- "$0")" && pwd)/replay.py" "$@"
