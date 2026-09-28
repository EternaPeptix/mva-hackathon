#!/usr/bin/env python3
"""Reproduce the descriptive BUB1B probe readout from GEO GSE9008.

The deposited experiment has one vehicle array and three resveratrol arrays.
This script reports their signal ratios but does not calculate a significance
test because there is no replicated vehicle variance.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import statistics
from pathlib import Path


GEO_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE9nnn/GSE9008/"
    "matrix/GSE9008_series_matrix.txt.gz"
)
EXPECTED_MATRIX_SHA256 = (
    "4015268331cdaa0a6b913115cd54883200a671c7922f965b3d4aef212e28f86c"
)
TARGET_PROBE = "203755_at"
EXPECTED_GSM_ORDER = ["GSM228717", "GSM228718", "GSM228719", "GSM228720"]


def parse_matrix(path: Path) -> dict:
    raw = path.read_bytes()
    sha256 = hashlib.sha256(raw).hexdigest()
    if sha256 != EXPECTED_MATRIX_SHA256:
        raise SystemExit(
            f"GSE9008 matrix SHA-256 mismatch: expected {EXPECTED_MATRIX_SHA256}, "
            f"got {sha256}"
        )

    metadata: dict[str, list[str]] = {}
    probe_rows: list[list[str]] = []
    in_table = False
    with gzip.GzipFile(fileobj=io.BytesIO(raw), mode="rb") as compressed:
        text = io.TextIOWrapper(compressed, encoding="utf-8", errors="replace")
        reader = csv.reader(text, delimiter="\t", quotechar='"')
        for row in reader:
            if not row:
                continue
            first = row[0]
            if first.startswith("!") and first not in (
                "!series_matrix_table_begin",
                "!series_matrix_table_end",
            ):
                metadata[first] = row[1:]
            if first == "!series_matrix_table_begin":
                in_table = True
                continue
            if first == "!series_matrix_table_end":
                break
            if in_table and first == TARGET_PROBE:
                probe_rows.append(row)

    if metadata.get("!Series_geo_accession") != ["GSE9008"]:
        raise SystemExit("Input is not the expected GSE9008 series matrix")
    if metadata.get("!Series_pubmed_id") != ["18089832"]:
        raise SystemExit("Unexpected primary-source link")
    if metadata.get("!Series_platform_id") != ["GPL570"]:
        raise SystemExit("Unexpected platform; target-probe mapping needs review")
    gsm_ids = metadata.get("!Sample_geo_accession", [])
    titles = metadata.get("!Sample_title", [])
    if gsm_ids != EXPECTED_GSM_ORDER or len(titles) != len(gsm_ids):
        raise SystemExit("Unexpected sample order or count; refusing silent remapping")
    if len(probe_rows) != 1 or len(probe_rows[0]) != len(gsm_ids) + 1:
        raise SystemExit("Expected exactly one BUB1B probe row and four samples")

    values = [float(value) for value in probe_rows[0][1:]]
    control = values[0]
    if control <= 0:
        raise SystemExit("Control intensity must be positive to calculate a ratio")
    treated_values = values[1:]
    ratios = [value / control for value in treated_values]
    observations = [
        {
            "sample": gsm_ids[index + 1],
            "title": titles[index + 1],
            "mas5_signal": value,
            "treated_to_single_control_ratio": ratio,
            "relative_change_percent": (ratio - 1.0) * 100.0,
        }
        for index, (value, ratio) in enumerate(zip(treated_values, ratios))
    ]
    return {
        "geo_series": "GSE9008",
        "primary_paper_pmid": "18089832",
        "platform": "GPL570",
        "target_gene": "BUB1B",
        "target_probe": TARGET_PROBE,
        "matrix_url": GEO_MATRIX_URL,
        "matrix_sha256": sha256,
        "matrix_bytes": len(raw),
        "processing_reported_by_geo": "MAS5.0",
        "resveratrol_concentration_um": 25,
        "exposure_hours": 48,
        "control_sample": {
            "sample": gsm_ids[0],
            "title": titles[0],
            "mas5_signal": control,
        },
        "treated_array_count": len(treated_values),
        "treated_to_single_control_ratio_mean": statistics.mean(ratios),
        "treated_to_single_control_ratio_range": [min(ratios), max(ratios)],
        "observations": observations,
        "statistical_test": None,
        "interpretation": (
            "All three treated-array signals are descriptively below the single "
            "control; the design cannot estimate control variance or support an "
            "inferential treatment comparison."
        ),
        "not_measured": [
            "BUBR1 protein abundance",
            "BUBR1 checkpoint function",
            "chromosome segregation fidelity",
            "muscle output or renewal",
            "MVA1 cells or either nominated patient variant",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True, help="GSE9008 series matrix .txt.gz")
    parser.add_argument("--out", type=Path, required=True, help="Output JSON path")
    args = parser.parse_args()
    result = parse_matrix(args.matrix)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"result": "PASS", "out": str(args.out), "observations": len(result["observations"])}, indent=2))


if __name__ == "__main__":
    main()
