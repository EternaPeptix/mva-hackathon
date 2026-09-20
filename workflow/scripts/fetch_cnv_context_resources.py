#!/usr/bin/env python3
"""Snapshot public GRCh38 SV population and dosage sources, preserving identities.

Publisher GCS generation and size are required for gnomAD. Whole gzip streams
are verified by the existing downloader. SHA256 records local content identity;
published MD5 is checked when present; CRC32C is recorded without rechecking.
"""
import argparse
import base64
import datetime
import hashlib
import json
import os
from pathlib import Path
import time
import urllib.request
from fetch_resource import fetch

GNOMAD = [
    ("gnomad.v4.1.sv.sites.vcf.gz", 1743544245, "1713312040472933", "crc32c=7RKeww=="),
    ("gnomad.v4.1.sv.sites.vcf.gz.tbi", 511659, "1713312039906181", "crc32c=JcrmLQ=="),
]
CLINGEN = ["ClinGen_gene_curation_list_GRCh38.tsv", "ClinGen_region_curation_list_GRCh38.tsv", "README"]


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8*1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def head_metadata(message):
    # GCS sends separate x-goog-hash fields; dict(items()) loses the first.
    result = {k.lower():v for k, v in message.items()}
    result["x-goog-hash"] = ",".join(message.get_all("x-goog-hash", []))
    return result


def verify_publisher_md5(path, metadata):
    expected = [v.split("=", 1)[1].strip() for v in metadata.get("x-goog-hash", "").split(",")
                if v.strip().startswith("md5=")]
    if not expected:
        return None
    if len(expected) != 1:
        raise ValueError("Ambiguous publisher MD5")
    h = hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8*1024*1024), b""):
            h.update(block)
    if base64.b64encode(h.digest()).decode("ascii") != expected[0]:
        raise ValueError("Publisher MD5 mismatch; file preserved")
    return h.hexdigest()


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save(path, data):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    os.replace(tmp, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if (out / "summary.json").exists():
        raise ValueError("Completed source snapshot preserved; select a new directory")
    os.nice(10)
    started, tick = now(), time.monotonic()
    files = []
    sources = [("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/" + n + "?generation=" + g,
                n, size, g, crc) for n, size, g, crc in GNOMAD]
    # Small mutable daily files are copied first and pinned by the resulting SHA.
    sources = [("https://ftp.clinicalgenome.org/" + n, n, None, None, None) for n in CLINGEN] + sources
    save(out / "status.json", {"status":"running", "started_at":started, "pid":os.getpid()})
    try:
        for url, name, expected_size, generation, crc in sources:
            req = urllib.request.Request(url, method="HEAD", headers={"Accept-Encoding":"identity"})
            with urllib.request.urlopen(req, timeout=45) as response:
                headers = head_metadata(response.headers)
                if response.status != 200:
                    raise ValueError("Unexpected source HTTP status")
            if generation and (headers.get("x-goog-generation") != generation or
                    int(headers["content-length"]) != expected_size or
                    crc not in headers.get("x-goog-hash", "")):
                raise ValueError("gnomAD publisher object metadata differs")
            path = fetch(url, out, expected_size=expected_size)
            if name.endswith(".tsv"):
                lines = path.read_text().splitlines()
                if not lines or not any(line.startswith("#") and "Genomic Location" in line for line in lines):
                    raise ValueError("ClinGen TSV header not recognized")
            receipt = json.loads((out / (name + ".verified.json")).read_text())
            files.append({"filename":name, "path":str(path), "bytes":path.stat().st_size,
                "sha256":digest(path), "source_url":url, "head_metadata":headers,
                "download_receipt":receipt, "generation_verified":bool(generation),
                "publisher_crc32c_verified":False, "verified_publisher_md5":verify_publisher_md5(path, headers),
                "integrity_scope":"Complete gzip block checksums when applicable; content length, local SHA256, and publisher MD5 when present; recorded publisher CRC32C is not independently checked."})
            save(out / "progress.json", {"started_at":started, "files":files})
        result = {"status":"complete", "started_at":started, "completed_at":now(),
            "elapsed_seconds":time.monotonic()-tick, "files":files,
            "script":{"path":str(Path(__file__).resolve()), "sha256":digest(Path(__file__))},
            "limitations":["Daily ClinGen URLs may later change; use these snapshot hashes when reproducing.",
                "gnomAD SV overlap alone is not allele equivalence, causality, benignity or absence of disease."]}
        save(out / "summary.json", result)
        save(out / "status.json", {"status":"complete", "completed_at":result["completed_at"],
            "summary_sha256":digest(out / "summary.json")})
        print(json.dumps({"status":"complete", "files":len(files), "elapsed_seconds":result["elapsed_seconds"]}))
    except BaseException as error:
        save(out / "status.json", {"status":"failed", "failed_at":now(), "error":str(error), "type":type(error).__name__})
        raise


if __name__ == "__main__":
    main()
