#!/usr/bin/env python3
"""Detect changed input files when resuming a completed workflow (metadata, not content hashes)."""
import hashlib
import json
from pathlib import Path
import sys

records = []
for name in sys.argv[1:]:
    path = Path(name).resolve()
    info = path.stat()
    records.append({'path': str(path), 'size': info.st_size, 'mtime_ns': info.st_mtime_ns})
print(hashlib.sha256(json.dumps(records, sort_keys=True).encode()).hexdigest())
