#!/usr/bin/env python3
"""Check current-year FDA novel approvals against older catalog aliases.

This supplements, rather than replaces, the source-bounded drug census. A name
match is a review cross-reference, not an authoritative product-identity merge.
"""
import argparse
import csv
import datetime
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import urllib.request

URL = "https://www.fda.gov/drugs/novel-drug-approvals-fda/novel-drug-approvals-2026"


def norm(name):
    # Explicit isotope typography only; never fuzzy-match different ingredients.
    return re.sub(r"\bf[\s-]*18\b", "f18", name.casefold().strip())


EXTRA_ALIASES = {
    "daraxonrasib": {
        "aliases": ["RMC-6236", "RMC 6236", "Rasonque"],
        "source": "https://www.fda.gov/news-events/press-announcements/fda-permits-expanded-access-investigational-pancreatic-cancer-drug",
    },
    "milsaperidone": {
        "aliases": ["Bysanti", "P-88-8991", "(S)-Hydroxy iloperidone"],
        "source": "https://precision.fda.gov/ginas/app/ui/substances/7SV1ZOG031",
    },
}


class Tables(HTMLParser):
    def __init__(self):
        super().__init__()
        self.tables = []
        self.table = None
        self.row = None
        self.cell = None

    def handle_starttag(self, tag, attrs):
        if tag == "table": self.table = []
        if tag == "tr": self.row = []
        if tag in ("td", "th"): self.cell = []

    def handle_data(self, text):
        if self.cell is not None: self.cell.append(text)

    def handle_endtag(self, tag):
        if tag in ("td", "th") and self.cell is not None:
            self.row.append(" ".join("".join(self.cell).split()))
            self.cell = None
        if tag == "tr" and self.row is not None and self.table is not None:
            self.table.append(self.row)
            self.row = None
        if tag == "table" and self.table is not None:
            self.tables.append(self.table)
            self.table = None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--html", type=Path, help="Saved public FDA page, otherwise fetch live")
    p.add_argument("--catalog-dir", type=Path, default=Path(__file__).resolve().parent)
    p.add_argument("--corpus", type=Path, help="Complete saved Europe PMC gene corpus for fresh-name screening")
    args = p.parse_args()
    raw = args.html.read_text() if args.html else urllib.request.urlopen(URL, timeout=60).read().decode()
    parser = Tables()
    parser.feed(raw)
    tables = [t for t in parser.tables if t and "Active Ingredient" in t[0]]
    if len(tables) != 1:
        raise RuntimeError("FDA table schema changed")
    table = tables[0]
    header = table[0]
    index = {}
    with (args.catalog_dir / "aliases.tsv").open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            index.setdefault(norm(row["alias"]), set()).add(row["catalog_id"])
    corpus = json.loads(args.corpus.read_text()) if args.corpus else None
    if corpus and (corpus["hitCount"] != corpus["retrievedCount"] or
                   corpus["retrievedCount"] != len(corpus["results"])):
        raise RuntimeError("Incomplete gene literature corpus")
    records = []
    today = datetime.datetime.now(datetime.timezone.utc).date()
    for cells in table[1:]:
        if len(cells) != len(header):
            raise RuntimeError("Incomplete FDA row")
        row = dict(zip(header, cells))
        name = row["Active Ingredient"]
        date = datetime.datetime.strptime(row["Approval Date"], "%m/%d/%Y").date()
        if date > today:
            raise RuntimeError("Future approval date in FDA table")
        matches = index.get(norm(name), set()) | index.get(norm(row["Drug Name"]), set())
        method = "ingredient or brand alias, isotope typography normalized"
        if not matches:
            base = re.sub(r"-[a-z]{4}$", "", name.casefold())
            matches = index.get(base, set())
            method = "biologic suffix removed; manual identity review required" if matches else "unmatched; retain FDA entry"
        components = []
        if " and " in name:
            components = [{"name": n, "catalog_matches": sorted(index.get(norm(n), set()))}
                          for n in name.split(" and ")]
            method = "separate component cross-references; not a combination-product identity merge"
        extra = EXTRA_ALIASES.get(name.casefold(), {})
        aliases = [name, row["Drug Name"], re.sub(r"-[a-z]{4}$", "", name)] + extra.get("aliases", [])
        pattern = re.compile(r"(?<!\w)(?:" + "|".join(re.escape(a) for a in aliases) + r")(?!\w)", re.I)
        hits = []
        if corpus:
            for article in corpus["results"]:
                haystack = article.get("title", "") + " " + article.get("abstractText", "")
                if pattern.search(haystack):
                    hits.append({"source": article.get("source"), "id": article.get("id"),
                                 "pmid": article.get("pmid"), "doi": article.get("doi")})
        use = next(v for k, v in row.items() if k.startswith("FDA-approved use"))
        triage = "OTHER_INDICATION_NO_BUB1B_EFFICACY_INFERRED"
        if any(w in use.lower() for w in ("positron emission", "in conjunction with mri")):
            triage = "DIAGNOSTIC_PRODUCT_NOT_RESCUE_EVIDENCE"
        if name.startswith("apitegromab"):
            triage = "DETAILED_MUSCLE_COMPARATOR_REVIEW"
        if name in ("brepocitinib", "gedatolisib") or name.startswith("garetosmab"):
            triage = "ADJACENT_PATHWAY_REVIEW_REQUIRED"
        records.append({"brand": row["Drug Name"], "active_ingredient": name,
                        "fda_approval_date": date.isoformat(), "catalog_matches": sorted(matches),
                        "match_method": method, "component_cross_references": components,
                        "search_aliases": aliases, "extra_alias_source": extra.get("source"),
                        "indication_triage": triage,
                        "literature_screen_state": "TITLE_ABSTRACT_ALIAS_SCREEN" if corpus else "NOT_RUN",
                        "literature_hits": hits, "efficacy_in_bub1b": "UNKNOWN"})
    numbers = [int(dict(zip(header, r))["No."].rstrip(".")) for r in table[1:]]
    if sorted(numbers) != list(range(1, len(records) + 1)):
        raise RuntimeError("FDA numbered list is incomplete")
    output = {"source_url": URL, "generated_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "source_mode": "saved HTML snapshot" if args.html else "live retrieval",
              "source_sha256": hashlib.sha256(raw.encode()).hexdigest(),
              "literature_corpus_sha256": hashlib.sha256(args.corpus.read_bytes()).hexdigest() if corpus else None,
              "scope": "CDER 2026 novel approvals only. Not all FDA products or worldwide approvals. Approval in another indication is not BUB1B efficacy.",
              "count": len(records), "records": records}
    (args.catalog_dir / "approval-overrides.json").write_text(json.dumps(output, indent=2) + "\n")
    print("FDA novel approvals", len(records), "matched", sum(bool(r["catalog_matches"]) for r in records),
          "unmatched", sum(not r["catalog_matches"] for r in records))


if __name__ == "__main__":
    main()
