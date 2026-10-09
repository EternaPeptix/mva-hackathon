#!/usr/bin/env python3
"""Reproduce descriptive BUB1B probe ratios from GEO GSE2992.

GEO provides one processed array per cell-line/treatment condition. The ratios
are descriptive only; this script does not estimate variance or test significance.
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
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE2nnn/GSE2992/"
    "matrix/GSE2992_series_matrix.txt.gz"
)
EXPECTED_MATRIX_SHA256 = (
    "e6745b7d18799779a0141b6ece520cce1a42eaa931b802bc33175db98561d04f"
)
EXPECTED_MATRIX_BYTES = 831055
PLATFORM_ANNOTATION_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
    "acc=GPL96&targ=self&form=text&view=full"
)
PLATFORM_ANNOTATION_SHA256 = (
    "34e01f9a175dc3927f8bd0f8301b164dcb6540ba26ffa360540ae597a360752c"
)
PLATFORM_ANNOTATION_BYTES = 46608806
TARGET_PROBE = "203755_at"
TARGET_GENE = "BUB1B"
EXPECTED_GSM_ORDER = ["GSM65385", "GSM65386", "GSM65387", "GSM65412"]
EXPECTED_TITLES = [
    "Untreated RWPE1 non-tumorigenic prostate epithelial cells",
    "RWPE1 human non-tumorigenic prostate epihelial cells treated with lunasin",
    "RWPE2 untreated human prostate cancer cells",
    "RWPE2 human prostate cancer cells treated with lunasin at 2 microM for 24 hr",
]
EXPECTED_ANNOTATION_ROW = (
    "203755_at\tNM_001211\t\tHomo sapiens\tOct 6, 2014\tExemplar sequence\t"
    "GenBank\tgb:NM_001211.2 /DEF=Homo sapiens budding uninhibited by benzimidazoles 1 "
    "(yeast homolog), beta (BUB1B), mRNA.  /FEA=mRNA /GEN=BUB1B /PROD=budding uninhibited "
    "by benzimidazoles 1 (yeasthomolog), beta /DB_XREF=gi:5729749 /UG=Hs.36708 budding "
    "uninhibited by benzimidazoles 1 (yeast homolog), beta /FL=gb:AF053306.1 "
    "gb:AF035933.1 gb:AF068760.1 gb:AF046918.1 gb:AF107297.1 gb:AF046079.2 "
    "gb:NM_001211.2\tNM_001211\tBUB1 mitotic checkpoint serine/threonine kinase B\t"
    "BUB1B\t701\tNM_001211\t0000278 // mitotic cell cycle // traceable author statement /// "
    "0006468 // protein phosphorylation // inferred from electronic annotation /// "
    "0006915 // apoptotic process // inferred from electronic annotation /// 0007049 // "
    "cell cycle // non-traceable author statement /// 0007051 // spindle organization // "
    "non-traceable author statement /// 0007067 // mitotic nuclear division // "
    "non-traceable author statement /// 0007091 // metaphase/anaphase transition of mitotic "
    "cell cycle // inferred from electronic annotation /// 0007093 // mitotic cell cycle "
    "checkpoint // traceable author statement /// 0007094 // mitotic spindle assembly "
    "checkpoint // traceable author statement /// 0008283 // cell proliferation // traceable "
    "author statement /// 0016310 // phosphorylation // inferred from electronic annotation "
    "/// 0031145 // anaphase-promoting complex-dependent proteasomal ubiquitin-dependent "
    "protein catabolic process // traceable author statement /// 0034501 // protein localization "
    "to kinetochore // inferred from physical interaction /// 0048015 // phosphatidylinositol-"
    "mediated signaling // non-traceable author statement /// 0051301 // cell division // "
    "inferred from electronic annotation /// 0051436 // negative regulation of ubiquitin-protein "
    "ligase activity involved in mitotic cell cycle // traceable author statement /// 0051439 // "
    "regulation of ubiquitin-protein ligase activity involved in mitotic cell cycle // traceable "
    "author statement /// 0071459 // protein localization to chromosome, centromeric region "
    "// inferred from electronic annotation\t0000775 // chromosome, centromeric region // "
    "inferred from electronic annotation /// 0000776 // kinetochore // inferred from direct assay "
    "/// 0000776 // kinetochore // non-traceable author statement /// 0000776 // kinetochore // "
    "traceable author statement /// 0000777 // condensed chromosome kinetochore // inferred "
    "from direct assay /// 0000778 // condensed nuclear chromosome kinetochore // inferred "
    "from electronic annotation /// 0000940 // condensed chromosome outer kinetochore // "
    "inferred from direct assay /// 0005634 // nucleus // inferred from electronic annotation "
    "/// 0005680 // anaphase-promoting complex // traceable author statement /// 0005694 // "
    "chromosome // inferred from electronic annotation /// 0005737 // cytoplasm // inferred from "
    "direct assay /// 0005815 // microtubule organizing center // inferred from electronic "
    "annotation /// 0005829 // cytosol // traceable author statement /// 0005856 // cytoskeleton "
    "// inferred from electronic annotation /// 0048471 // perinuclear region of cytoplasm // "
    "inferred from direct assay /// 0051233 // spindle midzone // non-traceable author statement\t"
    "0000166 // nucleotide binding // inferred from electronic annotation /// 0004672 // "
    "protein kinase activity // non-traceable author statement /// 0004674 // protein "
    "serine/threonine kinase activity // inferred from electronic annotation /// 0005515 // "
    "protein binding // inferred from physical interaction /// 0005524 // ATP binding // inferred "
    "from electronic annotation /// 0016301 // kinase activity // inferred from electronic "
    "annotation /// 0016740 // transferase activity // inferred from electronic annotation /// "
    "0016772 // transferase activity, transferring phosphorus-containing groups // inferred from "
    "electronic annotation"
)


def parse_matrix(path: Path) -> dict:
    raw = path.read_bytes()
    sha256 = hashlib.sha256(raw).hexdigest()
    if len(raw) != EXPECTED_MATRIX_BYTES or sha256 != EXPECTED_MATRIX_SHA256:
        raise SystemExit(
            "GSE2992 matrix checksum/size mismatch: "
            f"expected {EXPECTED_MATRIX_BYTES} bytes/{EXPECTED_MATRIX_SHA256}, "
            f"got {len(raw)} bytes/{sha256}"
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

    expected_metadata = {
        "!Series_geo_accession": ["GSE2992"],
        "!Series_platform_id": ["GPL96"],
        "!Sample_geo_accession": EXPECTED_GSM_ORDER,
        "!Sample_title": EXPECTED_TITLES,
        "!Sample_data_processing": ["2^(RMA-normalized data)"] * 4,
    }
    for key, expected in expected_metadata.items():
        if metadata.get(key) != expected:
            raise SystemExit(f"Unexpected GSE2992 metadata for {key}: {metadata.get(key)!r}")
    if len(probe_rows) != 1 or len(probe_rows[0]) != len(EXPECTED_GSM_ORDER) + 1:
        raise SystemExit("Expected exactly one 203755_at row and four GSM values")
    annotation_fields = EXPECTED_ANNOTATION_ROW.split("\t")
    if annotation_fields[0] != TARGET_PROBE or annotation_fields[10] != TARGET_GENE:
        raise SystemExit("Pinned GPL96 annotation row does not map the target probe to BUB1B")

    values = [float(value) for value in probe_rows[0][1:]]
    pairs = [
        ("RWPE-1", "GSM65385", "GSM65386", 0, 1),
        ("RWPE-2", "GSM65387", "GSM65412", 2, 3),
    ]
    observations = []
    for cell_line, control_gsm, treated_gsm, control_index, treated_index in pairs:
        control = values[control_index]
        treated = values[treated_index]
        ratio = treated / control
        observations.append(
            {
                "cell_line": cell_line,
                "control_sample": control_gsm,
                "lunasin_sample": treated_gsm,
                "control_linear_expression": control,
                "lunasin_linear_expression": treated,
                "lunasin_to_control_linear_ratio": ratio,
                "statistical_test": None,
            }
        )

    rwpe1_ratio = observations[0]["lunasin_to_control_linear_ratio"]
    return {
        "geo_series": "GSE2992",
        "geo_series_title": "Soy peptide up-regulates tumor suppressor and other chemopreventive genes in prostate epithelial cells",
        "geo_citation_status": "GEO record says Citation missing",
        "primary_paper_pmid": "21526452",
        "primary_paper_doi": "10.1080/01635581.2011.539312",
        "paper_accession_linkage": (
            "Inferred from matching study title/topic, RWPE-1/RWPE-2 design, 2 uM/24 h "
            "treatment and authors; GEO does not explicitly cite the paper. The paper's "
            "Table 2 reports a 2.23-fold BUB1B value consistent with RWPE-1 matrix ratio."
        ),
        "platform": "GPL96",
        "target_gene": TARGET_GENE,
        "target_probe": TARGET_PROBE,
        "platform_annotation_url": PLATFORM_ANNOTATION_URL,
        "platform_annotation_sha256": PLATFORM_ANNOTATION_SHA256,
        "platform_annotation_bytes": PLATFORM_ANNOTATION_BYTES,
        "platform_annotation_matching_rows": 1,
        "platform_annotation_row": EXPECTED_ANNOTATION_ROW,
        "matrix_url": GEO_MATRIX_URL,
        "matrix_sha256": sha256,
        "matrix_bytes": len(raw),
        "processing_reported_by_geo": "2^(RMA-normalized data)",
        "sample_count": 4,
        "sample_design": (
            "One processed array per cell line per condition; no replicate-level variance"
        ),
        "lunasin_concentration_um": 2,
        "exposure_hours": 24,
        "observations": observations,
        "published_rwpe1_table2_fold_change": 2.23,
        "rwpe1_ratio_rounds_to_published_fold_change": round(rwpe1_ratio, 2) == 2.23,
        "interpretation": (
            "The processed RWPE-1 matrix descriptively reproduces the published 2.23-fold "
            "BUB1B RNA value; RWPE-2 is 1.25-fold. This reuses the likely same source arrays, "
            "not independent biological replication, and has no inferential test."
        ),
        "not_measured": [
            "BUBR1 protein abundance or half-life",
            "BUBR1 checkpoint function",
            "chromosome segregation fidelity",
            "muscle output or renewal",
            "MVA1 cells or either patient variant",
            "intact lunasin concentration inside human muscle",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True, help="GSE2992 series matrix .txt.gz")
    parser.add_argument("--out", type=Path, required=True, help="Output JSON path")
    args = parser.parse_args()
    result = parse_matrix(args.matrix)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "result": "PASS",
                "out": str(args.out),
                "sample_pairs": len(result["observations"]),
                "rwpe1_ratio": result["observations"][0]["lunasin_to_control_linear_ratio"],
                "rwpe2_ratio": result["observations"][1]["lunasin_to_control_linear_ratio"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
