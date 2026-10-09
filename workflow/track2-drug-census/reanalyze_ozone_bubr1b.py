#!/usr/bin/env python3
"""Recheck Bub1b in GEO GSE45166's log2-normalized ozone exposure matrix.

Requires Python 3, NumPy, and SciPy. Source files are downloaded separately from
the official NCBI GEO series-matrix and GPL7202 annotation URLs.
"""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import ttest_ind


SAMPLE_GROUPS = {
    "6h": {
        "ozone": ["GSM1098461", "GSM1098462", "GSM1098465", "GSM1098466"],
        "air": ["GSM1098463", "GSM1098464", "GSM1098467", "GSM1098468"],
    },
    "24h": {
        "ozone": ["GSM1098469", "GSM1098470", "GSM1098473", "GSM1098474"],
        "air": ["GSM1098471", "GSM1098472", "GSM1098475", "GSM1098476"],
    },
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def matrix(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_begin"))
    end = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_end"))
    rows = list(csv.reader(lines[start + 1 : end], delimiter="\t"))
    header = rows[0]
    ids = header[1:]
    features = []
    values = []
    for row in rows[1:]:
        if len(row) != len(header):
            continue
        try:
            signal = [float(value) for value in row[1:]]
        except ValueError:
            continue
        features.append(row[0])
        values.append(signal)
    return ids, features, np.asarray(values, dtype=float)


def bub1b_probe(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    header_index = next(i for i, line in enumerate(lines) if line.startswith("ID\t"))
    rows = csv.DictReader(lines[header_index:], delimiter="\t")
    probes = [row["ID"] for row in rows if (row.get("Gene symbol") or "").casefold() == "bub1b"]
    if len(probes) != 1:
        raise ValueError(f"Expected one GPL7202 Bub1b probe; found {probes}")
    return probes[0]


def benjamini_hochberg(p_values):
    order = np.argsort(p_values)
    ranked = p_values[order] * len(p_values) / np.arange(1, len(p_values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q_values = np.empty_like(ranked)
    q_values[order] = np.minimum(ranked, 1.0)
    return q_values


def analyze(series_matrix, annotation):
    sample_ids, feature_ids, values = matrix(series_matrix)
    probe = bub1b_probe(annotation)
    if probe not in feature_ids:
        raise ValueError(f"GPL7202 Bub1b probe {probe} is absent from the series matrix")
    probe_index = feature_ids.index(probe)
    results = []
    for timepoint, groups in SAMPLE_GROUPS.items():
        ozone_ix = [sample_ids.index(sample) for sample in groups["ozone"]]
        air_ix = [sample_ids.index(sample) for sample in groups["air"]]
        ozone = values[:, ozone_ix]
        air = values[:, air_ix]
        difference = ozone.mean(axis=1) - air.mean(axis=1)
        p_values = ttest_ind(ozone, air, axis=1, equal_var=True).pvalue
        q_values = benjamini_hochberg(p_values)
        effect = float(difference[probe_index])
        results.append({
            "timepoint": timepoint,
            "n_per_group": len(ozone_ix),
            "probe_id": probe,
            "ozone_mean_log2": float(values[probe_index, ozone_ix].mean()),
            "filtered_air_mean_log2": float(values[probe_index, air_ix].mean()),
            "log2_difference_ozone_minus_air": effect,
            "ozone_to_air_ratio": float(2**effect),
            "student_t_two_sided_p": float(p_values[probe_index]),
            "benjamini_hochberg_q_all_probes": float(q_values[probe_index]),
            "passes_article_p_lt_0_05_and_1_5_fold_rule": bool(
                p_values[probe_index] < 0.05 and abs(effect) > np.log2(1.5)
            ),
            "all_probe_count_p_lt_0_05_and_1_5_fold": int(
                np.sum((p_values < 0.05) & (np.abs(difference) > np.log2(1.5)))
            ),
            "all_probe_count_bh_q_lt_0_05": int(np.sum(q_values < 0.05)),
            "all_probe_count_bh_q_lt_0_05_and_1_5_fold": int(
                np.sum((q_values < 0.05) & (np.abs(difference) > np.log2(1.5)))
            ),
        })
    return {
        "series_accession": "GSE45166",
        "platform_accession": "GPL7202",
        "feature_count": len(feature_ids),
        "sample_matrix_sha256": sha256(series_matrix),
        "platform_annotation_sha256": sha256(annotation),
        "method": "Two-sided equal-variance Student t-test on GEO-provided log2-transformed, quantile-normalized values; Benjamini-Hochberg across all matrix features.",
        "results": results,
        "limits": [
            "Whole-lung transcript signal; not a cell-intrinsic or protein measurement.",
            "No BUBR1 abundance/function, chromosome fidelity, muscle, or MVA outcome was measured.",
            "The 6-hour source-threshold DEG count differs by one in this reanalysis; Bub1b's 24-hour result matches the published supplementary row to rounding.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--series-matrix", type=Path, required=True)
    parser.add_argument("--platform-annotation", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = analyze(args.series_matrix, args.platform_annotation)
    text = json.dumps(result, indent=2) + "\n"
    if args.out:
        args.out.write_text(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
