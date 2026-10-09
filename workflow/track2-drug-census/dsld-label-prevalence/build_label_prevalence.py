#!/usr/bin/env python3
"""Count exact DietarySupplementFacts ingredient strings across NIH DSLD labels."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sqlite3
import sys
import tempfile
import zipfile
from collections import Counter
from pathlib import Path


HERE = Path(__file__).resolve().parent
MANIFEST_PATH = HERE / "source-manifest.json"
CANDIDATES_PATH = HERE / "candidate-terms.tsv"
EXPECTED_MARKET_COUNTS = {"On Market": 121_959, "Off Market": 92_821}
EXPECTED_FACT_ROWS = 2_020_128
EXPECTED_ON_MARKET_FACT_LABELS = 121_949


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def digest_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def csv_rows(archive: zipfile.ZipFile, name: str):
    with archive.open(name) as raw:
        wrapper = io.TextIOWrapper(raw, encoding="utf-8-sig", newline="")
        yield from csv.DictReader(wrapper)


def source_shards(archive: zipfile.ZipFile, prefix: str) -> list[str]:
    return sorted(
        name for name in archive.namelist()
        if name.rsplit("/", 1)[-1].startswith(prefix)
        and name.endswith(".csv")
    )


def load_candidates() -> list[dict[str, str]]:
    with CANDIDATES_PATH.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--zip", required=True, type=Path, help="NIH DSLD CSV export ZIP")
    parser.add_argument("--out-dir", type=Path, default=HERE)
    args = parser.parse_args()

    source = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    actual_sha = digest_file(args.zip)
    require(actual_sha == source["sha256"], "Input ZIP SHA-256 differs from pinned NIH export")
    require(args.zip.stat().st_size == source["zip_bytes"], "Input ZIP size differs from pinned NIH export")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    market_counts: Counter[str] = Counter()
    total_labels = 0
    product_id_count = 0
    fact_rows = 0
    fact_id_count = 0
    on_market_fact_labels = 0
    ingredient_count_rows = 0
    candidate_terms = load_candidates()

    with tempfile.NamedTemporaryFile(prefix="dsld-label-prevalence-", suffix=".sqlite") as temp_db:
        db = sqlite3.connect(temp_db.name)
        db.execute("PRAGMA journal_mode=OFF")
        db.execute("PRAGMA synchronous=OFF")
        db.execute("PRAGMA temp_store=MEMORY")
        db.executescript(
            """
            CREATE TABLE product (
                dsld_id TEXT PRIMARY KEY,
                market_status TEXT NOT NULL,
                brand_name TEXT
            );
            CREATE TABLE fact (
                dsld_id TEXT NOT NULL,
                ingredient TEXT NOT NULL,
                PRIMARY KEY (dsld_id, ingredient)
            );
            CREATE INDEX fact_ingredient_idx ON fact(ingredient);
            """
        )

        with zipfile.ZipFile(args.zip) as archive:
            overview_shards = source_shards(archive, "ProductOverview_")
            fact_shards = source_shards(archive, "DietarySupplementFacts_")
            readme_shards = source_shards(archive, "ReadMe_")
            require(len(overview_shards) == 8, f"Expected 8 product shards, got {len(overview_shards)}")
            require(len(fact_shards) == 8, f"Expected 8 facts shards, got {len(fact_shards)}")
            require(len(readme_shards) == 8, f"Expected 8 README shards, got {len(readme_shards)}")

            readmes = [archive.read(name).decode("utf-8", "replace") for name in readme_shards]
            for readme in readmes:
                require("Labels Downloaded,214780" in readme, "Source README label total mismatch")
                require("Date of download,2026-01-10T" in readme, "Source README date mismatch")
                require("Version 9.4.0 | January 2026" in readme, "Source README version mismatch")

            product_batch = []
            for shard in overview_shards:
                for row in csv_rows(archive, shard):
                    dsld_id = row.get("DSLD ID", "").strip()
                    if not dsld_id:
                        continue
                    status = row.get("Market Status", "").strip()
                    require(status in EXPECTED_MARKET_COUNTS, f"Unexpected market status: {status!r}")
                    product_batch.append((dsld_id, status, row.get("Brand Name", "").strip()))
                    market_counts[status] += 1
                    total_labels += 1
                    if len(product_batch) >= 10_000:
                        db.executemany("INSERT INTO product VALUES (?, ?, ?)", product_batch)
                        product_batch.clear()
            if product_batch:
                db.executemany("INSERT INTO product VALUES (?, ?, ?)", product_batch)
            product_id_count = db.execute("SELECT COUNT(*) FROM product").fetchone()[0]
            require(total_labels == source["download_readme_label_count"], "Product-row count mismatch")
            require(product_id_count == total_labels, "Duplicate or missing product label IDs")
            require(dict(market_counts) == EXPECTED_MARKET_COUNTS, f"Market status totals mismatch: {market_counts}")

            fact_batch = []
            for shard in fact_shards:
                for row in csv_rows(archive, shard):
                    dsld_id = row.get("DSLD ID", "").strip()
                    ingredient = row.get("Ingredient", "").strip()
                    if not dsld_id or not ingredient:
                        continue
                    fact_batch.append((dsld_id, ingredient))
                    fact_rows += 1
                    if len(fact_batch) >= 10_000:
                        db.executemany("INSERT OR IGNORE INTO fact VALUES (?, ?)", fact_batch)
                        fact_batch.clear()
            if fact_batch:
                db.executemany("INSERT OR IGNORE INTO fact VALUES (?, ?)", fact_batch)
            require(fact_rows == EXPECTED_FACT_ROWS, f"Facts row count mismatch: {fact_rows}")
            require(
                db.execute(
                    "SELECT COUNT(*) FROM fact f LEFT JOIN product p USING(dsld_id) WHERE p.dsld_id IS NULL"
                ).fetchone()[0] == 0,
                "Facts rows reference unknown label IDs",
            )
            fact_id_count = db.execute("SELECT COUNT(DISTINCT dsld_id) FROM fact").fetchone()[0]
            on_market_fact_labels = db.execute(
                "SELECT COUNT(DISTINCT f.dsld_id) FROM fact f JOIN product p USING(dsld_id) WHERE p.market_status='On Market'"
            ).fetchone()[0]
            require(on_market_fact_labels == EXPECTED_ON_MARKET_FACT_LABELS, "On-market facts label count mismatch")

        output_rows = db.execute(
            """
            SELECT f.ingredient,
                   COUNT(*) AS all_market_labels,
                   SUM(CASE WHEN p.market_status='On Market' THEN 1 ELSE 0 END) AS on_market_labels,
                   SUM(CASE WHEN p.market_status='Off Market' THEN 1 ELSE 0 END) AS off_market_labels,
                   COUNT(DISTINCT CASE WHEN p.market_status='On Market' AND p.brand_name != '' THEN p.brand_name END) AS on_market_brands
            FROM fact f JOIN product p USING(dsld_id)
            GROUP BY f.ingredient
            ORDER BY on_market_labels DESC, f.ingredient COLLATE NOCASE ASC, f.ingredient ASC
            """
        ).fetchall()
        by_literal = {row[0]: row for row in output_rows}
        ingredient_count_rows = len(output_rows)

        prevalence_path = args.out_dir / "ingredient-label-prevalence.tsv"
        with prevalence_path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow([
                "ingredient_literal", "all_market_labels", "on_market_labels",
                "off_market_labels", "distinct_on_market_brands",
            ])
            for row in output_rows:
                writer.writerow(row)

        candidate_path = args.out_dir / "candidate-prevalence.tsv"
        missing_candidates = []
        with candidate_path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow([
                "ingredient_literal", "candidate_or_comparator", "role_or_relationship",
                "all_market_labels", "on_market_labels", "off_market_labels",
                "distinct_on_market_brands",
            ])
            for item in candidate_terms:
                row = by_literal.get(item["ingredient_literal"])
                if row is None:
                    missing_candidates.append(item["ingredient_literal"])
                    values = (0, 0, 0, 0)
                else:
                    values = row[1:]
                writer.writerow([
                    item["ingredient_literal"], item["candidate_or_comparator"],
                    item["role_or_relationship"], *values,
                ])

        summary = {
            "schema_version": 1,
            "source_manifest": "source-manifest.json",
            "source_zip_sha256": actual_sha,
            "source_version": source["source_version"],
            "source_export_timestamp": source["source_export_timestamp"],
            "retrieved_date": source["retrieved_date"],
            "label_counts": {
                "all_labels": total_labels,
                "on_market": market_counts["On Market"],
                "off_market": market_counts["Off Market"],
                "unique_product_ids": product_id_count,
                "on_market_label_records_with_supplement_facts": on_market_fact_labels,
                "all_label_records_with_supplement_facts": fact_id_count,
            },
            "ingredient_facts": {
                "source_rows": fact_rows,
                "distinct_label_ingredient_pairs": db.execute("SELECT COUNT(*) FROM fact").fetchone()[0],
                "distinct_exact_ingredient_strings": ingredient_count_rows,
                "source_strings_with_on_market_labels": sum(1 for row in output_rows if row[2] > 0),
            },
            "candidate_terms": len(candidate_terms),
            "candidate_literals_with_zero_on_market_labels": missing_candidates,
            "counts_are": "Distinct DSLD label IDs containing each exact Ingredient string in DietarySupplementFacts.",
            "counts_are_not": [
                "sales, purchases, product availability in current retail stock, consumers, or exposures",
                "verified chemical identities for broad label names",
                "ingredient amounts, dose feasibility, safety, BUB1B efficacy, or treatment benefit",
            ],
        }
        (args.out_dir / "summary.json").write_text(
            json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )

        db.close()

    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, sqlite3.Error, zipfile.BadZipFile) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
