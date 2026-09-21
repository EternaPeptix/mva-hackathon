#!/usr/bin/env bash
set -euo pipefail
if [[ $# != 1 || -e "$1" ]]; then
  echo "Usage: bash run.sh NEW_OUTPUT_DIRECTORY (must not exist)" >&2
  exit 2
fi
src="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$1/inputs" "$1/results"
out="$(cd "$1" && pwd)"
python="${MVA_RNA_PYTHON:-python3}"
rscript="${MVA_RNA_RSCRIPT:-Rscript}"
"$python" - "$out/inputs" "$src/sources.json" <<'PY'
import hashlib, json, pathlib, sys, urllib.request
dest = pathlib.Path(sys.argv[1])
expected = {x['file']: x['sha256'] for x in json.loads(pathlib.Path(sys.argv[2]).read_text())}
urls = [
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE134nnn/GSE134781/suppl/GSE134781_counts.txt.gz',
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE134nnn/GSE134781/soft/GSE134781_family.soft.gz',
    'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Mm/mh.all.v2025.1.Mm.symbols.gmt',
]
manifest = []
for url in urls:
    data = urllib.request.urlopen(url, timeout=90).read()
    name = url.rsplit('/', 1)[1]
    digest = hashlib.sha256(data).hexdigest()
    if digest != expected[name]:
        raise SystemExit('Public input changed: ' + name + '; review before replay')
    (dest / name).write_bytes(data)
    manifest.append(dict(file=name, url=url, bytes=len(data), sha256=digest))
(dest / 'sources.json').write_text(json.dumps(manifest, indent=2) + '\n')
PY
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
"$rscript" "$src/analysis.R" "$out/inputs" "$out/results"
