#!/usr/bin/env python3
"""Read public target indexes; absent annotations do not imply drug inefficacy."""
import argparse
import datetime
import json
from pathlib import Path
import urllib.request


def query(endpoint, text):
    request = urllib.request.Request(
        endpoint, data=json.dumps({"query": text}).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "MVA-public-research/1.0"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        result = json.load(response)
    if result.get("errors") or "data" not in result:
        raise RuntimeError(json.dumps(result))
    return {"endpoint": endpoint, "query": text, "response": result["data"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    records = []
    # BUB1 is a distinct target; EGFR is an annotation-positive service control.
    for symbol, hgnc, ensembl in [
        ("BUB1B", "hgnc:1149", "ENSG00000156970"),
        ("BUB1", "hgnc:1148", "ENSG00000169679"),
        ("EGFR", "hgnc:3236", "ENSG00000146648"),
    ]:
        dgidb = query("https://dgidb.org/api/graphql", '''{
          gene(conceptId:"%s") { name conceptId interactions {
            drug {name conceptId approved} publications {pmid}
            sources {sourceDbName sourceDbVersion license licenseLink}
          }}
        }''' % hgnc)
        opentargets = query("https://api.platform.opentargets.org/api/v4/graphql", '''{
          meta {dataVersion {year month iteration}}
          target(ensemblId:"%s") { id approvedSymbol
            drugAndClinicalCandidates {count rows {maxClinicalStage drug {id name}}}
            chemicalProbes {id targetFromSourceId drugId mechanismOfAction isHighQuality urls {url niceName}}
          }
        }''' % ensembl)
        gene = dgidb["response"]["gene"]
        target = opentargets["response"]["target"]
        if not gene or gene["name"] != symbol or not target or target["approvedSymbol"] != symbol:
            raise RuntimeError("Target identity did not resolve: " + symbol)
        records.append({"symbol": symbol, "dgidb": dgidb, "opentargets": opentargets})
        print(symbol, "DGIdb interactions", len(gene["interactions"]),
              "Open Targets clinical records", target["drugAndClinicalCandidates"]["count"],
              "probes", len(target["chemicalProbes"]))
    positive = records[-1]
    if not positive["dgidb"]["response"]["gene"]["interactions"]:
        raise RuntimeError("DGIdb positive annotation control failed")
    if not positive["opentargets"]["response"]["target"]["drugAndClinicalCandidates"]["count"]:
        raise RuntimeError("Open Targets positive annotation control failed")
    payload = {
        "retrieved_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "scope": "Public index annotations, not an efficacy screen. BUB1 is not BUB1B. Empty results do not exclude uncurated or indirect rescue mechanisms.",
        "records": [
            {
                "symbol": item["symbol"],
                "dgidb": {
                    "endpoint": item["dgidb"]["endpoint"],
                    "query": item["dgidb"]["query"],
                    "resolved_concept": item["dgidb"]["response"]["gene"]["conceptId"],
                    "interaction_count": len(item["dgidb"]["response"]["gene"]["interactions"]),
                },
                "opentargets": {
                    "endpoint": item["opentargets"]["endpoint"],
                    "query": item["opentargets"]["query"],
                    "data_version": item["opentargets"]["response"]["meta"]["dataVersion"],
                    "resolved_id": item["opentargets"]["response"]["target"]["id"],
                    "clinical_record_count": item["opentargets"]["response"]["target"]["drugAndClinicalCandidates"]["count"],
                    "chemical_probe_count": len(item["opentargets"]["response"]["target"]["chemicalProbes"]),
                },
            } for item in records
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
