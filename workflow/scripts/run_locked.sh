#!/usr/bin/env bash
# run_locked.sh — run a command at most once, even across shells and sessions.
#
# Written after a real incident in this project: two `bwa index` processes were started
# minutes apart against the same reference and both wrote to the same output prefix.
# bwa does not lock its outputs, so the .bwt/.sa files were interleaved garbage that
# would have produced silently wrong alignments. A first kill attempt missed one PID
# because the duplicate had been reparented, and the corruption resumed.
#
# macOS ships no flock(1), so we use mkdir as the atomic primitive: mkdir on an existing
# directory fails atomically on every POSIX filesystem, which is exactly a test-and-set.
#
# Usage: run_locked.sh <lock-name> <command...>
#   run_locked.sh bwa-index bwa index /path/to/ref.fa
#
# Exit codes: 0 = command ran and succeeded, or a live holder already owns the lock;
#             non-zero = the command itself failed.
set -euo pipefail

LOCKROOT="${LOCKROOT:-/Users/jeweled/mva-hackathon/.locks}"
mkdir -p "$LOCKROOT"

name="${1:?usage: run_locked.sh <lock-name> <command...>}"; shift
[ "$#" -gt 0 ] || { echo "run_locked: no command given" >&2; exit 64; }

lockdir="$LOCKROOT/$name.lock"
pidfile="$lockdir/pid"

if ! mkdir "$lockdir" 2>/dev/null; then
  # Lock exists. Distinguish a live holder from a stale lock left by a killed process.
  if [ -f "$pidfile" ] && kill -0 "$(cat "$pidfile" 2>/dev/null)" 2>/dev/null; then
    echo "run_locked: '$name' already running as PID $(cat "$pidfile"); refusing to start a second copy."
    exit 0
  fi
  echo "run_locked: stale lock for '$name' (holder is gone); reclaiming."
  rm -rf "$lockdir"
  mkdir "$lockdir"
fi

echo $$ > "$pidfile"
# Release the lock however we exit, including on SIGTERM/SIGINT.
trap 'rm -rf "$lockdir"' EXIT INT TERM

echo "run_locked: [$name] starting: $*"
"$@"
rc=$?
echo "run_locked: [$name] finished with status $rc"
exit $rc
