#!/usr/bin/env python3
"""Reproduce the descriptive BUB1B probe readout from GEO GSE15162.

This is not a differential-expression test: the GEO series matrix contains one
processed array for each treatment/timepoint, so replicate variance is absent.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path


GEO_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE15nnn/GSE15162/"
    "matrix/GSE15162_series_matrix.txt.gz"
)
EXPECTED_MATRIX_SHA256 = (
    "b85a836ead7ad65023dd1d53b5dc6ac304ea45e299884baee2d64f188f4b82b1"
)
TARGET_PROBE = "203755_at"
TARGET_GENE = "BUB1B"
EXPECTED_GSM_ORDER = [
    "GSM378688",
    "GSM378689",
    "GSM378690",
    "GSM378691",
    "GSM378692",
    "GSM378693",
    "GSM378694",
    "GSM378695",
]


def parse_matrix(path: Path) -> dict:
    raw = path.read_bytes()
    sha256 = hashlib.sha256(raw).hexdigest()
    if sha256 != EXPECTED_MATRIX_SHA256:
        raise SystemExit(
            f"GSE15162 matrix SHA-256 mismatch: expected {EXPECTED_MATRIX_SHA256}, "
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

    if metadata.get("!Series_geo_accession") != ["GSE15162"]:
        raise SystemExit("Input is not the expected GSE15162 series matrix")
    if metadata.get("!Series_platform_id") != ["GPL570"]:
        raise SystemExit("Unexpected platform; target-probe mapping needs review")
    gsm_ids = metadata.get("!Sample_geo_accession", [])
    titles = metadata.get("!Sample_title", [])
    if gsm_ids != EXPECTED_GSM_ORDER or len(titles) != len(gsm_ids):
        raise SystemExit("Unexpected sample order or count; refusing silent remapping")
    if len(probe_rows) != 1 or len(probe_rows[0]) != len(gsm_ids) + 1:
        raise SystemExit("Expected exactly one BUB1B probe row and eight samples")

    values = [float(value) for value in probe_rows[0][1:]]
    observations = []
    for index in range(0, len(values), 2):
        control_title = titles[index]
        treated_title = titles[index + 1]
        if "control" not in control_title.casefold() or "quercetin" not in treated_title.casefold():
            raise SystemExit("Sample titles do not alternate control/quercetin as expected")
        time_label = control_title.rsplit(" ", 1)[-1]
        if treated_title.rsplit(" ", 1)[-1] != time_label or not time_label.endswith("h"):
            raise SystemExit("Control and treated sample timepoints do not match")
        control = values[index]
        treated = values[index + 1]
        difference = treated - control
        observations.append(
            {
                "time_h": int(time_label[:-1]),
                "control_sample": gsm_ids[index],
                "quercetin_sample": gsm_ids[index + 1],
                "log2_expression_control": control,
                "log2_expression_quercetin": treated,
                "log2_difference_quercetin_minus_control": difference,
                "quercetin_to_control_linear_ratio": 2**difference,
                "statistical_test": None,
            }
        )

    return {
        "geo_series": "GSE15162",
        "primary_paper_pmid": "21782406",
        "platform": "GPL570",
        "target_gene": TARGET_GENE,
        "target_probe": TARGET_PROBE,
        "matrix_url": GEO_MATRIX_URL,
        "matrix_sha256": sha256,
        "matrix_bytes": len(raw),
        "processing_reported_by_geo": (
            "RMA normalization using raw data with GeneSpring GX v9.0"
        ),
        "sample_count": len(gsm_ids),
        "sample_design": "one processed array per treatment/timepoint in the deposited series matrix",
        "quercetin_concentration_um": 3,
        "observations": observations,
        "interpretation": (
            "Small, mixed-direction descriptive transcript shifts; no replicate-level "
            "variance or inferential test can be computed from this deposited matrix."
        ),
        "not_measured": [
            "BUBR1 protein abundance",
            "BUBR1 checkpoint function",
            "chromosome segregation fidelity",
            "muscle output or renewal",
            "MVA1 cells or either patient variant",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True, help="GSE15162 series matrix .txt.gz")
    parser.add_argument("--out", type=Path, required=True, help="Output JSON path")
    args = parser.parse_args()
    result = parse_matrix(args.matrix)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"result": "PASS", "out": str(args.out), "observations": len(result["observations"])}, indent=2))


if __name__ == "__main__":
    main()
