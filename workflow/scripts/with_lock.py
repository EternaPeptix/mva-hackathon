#!/usr/bin/env python3
"""Hold a POSIX lock for the lifetime of one foreground command."""
import fcntl
from pathlib import Path
import subprocess
import sys

if len(sys.argv) < 3:
    raise SystemExit('usage: with_lock.py LOCK_FILE COMMAND [ARG ...]')
path = Path(sys.argv[1])
path.parent.mkdir(parents=True, exist_ok=True)
with path.open('a+') as lock:
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit('another process holds this workflow lock')
    # The child keeps the same open-file lock if this wrapper is interrupted or killed.
    raise SystemExit(subprocess.run(sys.argv[2:], pass_fds=(lock.fileno(),)).returncode)
