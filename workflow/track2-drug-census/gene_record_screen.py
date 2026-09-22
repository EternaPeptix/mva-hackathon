#!/usr/bin/env python3
"""Reproduce the bounded, catalog-independent priority screen of 2026-09-21.

Regex matches are screening priorities, never efficacy/adjudication labels.
The fixed read ledger records actual reviewer actions; rerunning does not cause
an automated reader to become a manual reviewer. No clinical catalog is loaded.
"""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
import re

REVIEWED_CORPUS_SHA256 = "1df717bce779851e8f8f02654cc476de320cf368ad3cbc10dad1499b5a809e90"

BENEFIT = r"rescu|restor|protect|amelior|improv|alleviat|revers|attenuat|prevent|supplement"
INTERVENTION = r"treat|administr|supplement|compound|chemical|pharmacolog|inhibitor|activator|agonist|antagonist|extract|peptide|antioxidant|diet|vitamin|hormon"
SENSITIVITY_TITLE = r"oocyt|ageing|aging|senescen|muscle|proger|aneuploid|pharmacolog|drug|treat|antioxidant|diet|vitamin|hormon|supplement|compound|extract"

# All these 21 saved abstracts were actually read completely in this pass.
READ_LEDGER = {
    "42286606": "BDNF/MSC secretome/7,8-DHF: contextual chromosome and developmental benefit; not BUB1B-mutant rescue",
    "27598149": "Phyllanthus emblica extract: normal colon cell division/genomic outcomes; mixture not a defined drug",
    "31480299": "Melatonin: cryopreservation oocyte protection; BUBR1 mediation unproved",
    "39506684": "Folic acid: bovine cell apoptosis protection; BUB1B network/marker relationship",
    "30693443": "EUK-207: oxidative cochlear injury and SAMP8 hearing protection; not BUB1B-mutant treatment",
    "27284927": "LiCl: Axin1-depleted oocyte rescue; no FRZB or BUB1B-deficiency rescue",
    "22886373": "Pifithrin-alpha/NAC: chromium stress effects; release from arrest not proof of genomic repair",
    "40470792": "BUBR1 expression restoration in neonatal HIE; no named drug rescue established by abstract",
    "37201634": "Aristolochic acid I is an ovarian injury agent; not therapeutic rescue",
    "36539842": "HtrA2 overexpression genetic rescue in aged oocytes; not a drug result",
    "32078082": "Putrescine/polyamine deprivation in arginase-deficient CHO cells; culture requirement, not MVA pharmacological rescue",
    "19846757": "NAMPT overexpression protected endothelial polyploidy via SIRT1; genetic result, not niacin or SIRT2 drug result",
    "40750410": "Propargylglycine blocks engineered MEGL-induced cancer cell death; not desired normal-cell/MVA rescue",
    "37782908": "PAF(C-16) reverses BUB1B knockdown effects on cancer progression; malignant behavior not a therapeutic rescue",
    "39845481": "Ascorbic acid/NAC: Sod1-null oocyte chromosome segregation benefit; different causal gene",
    "21550059": "Aged vascular smooth muscle and ROS/BUBR1 signaling; kinase/oxidase inhibitors block induction, not restoration",
    "26883501": "BubR1-H/H blood-brain barrier impairment; descriptive no administered rescue drug",
    "28229533": "NF-kB genetic reduction benefits ERCC1 model; BubR1 mice establish MDSC phenotype only",
    "41579174": "Gentamicin nephrotoxicity transcriptomics and drug database suggestions; no new protective drug experiment",
    "42000973": "MLH3 depletion and oocyte/endometrial pathology; no drug rescue",
    "29067790": "SIRT2/BUBR1 K243R overexpression oocyte rescue; genetic intervention, not pharmacological",
}
FULLTEXT_PASSAGES = {
    "42286606": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13501573/fullTextXML",
    "27598149": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5037716/fullTextXML",
    "31480299": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6770451/fullTextXML",
    "39506684": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11542246/fullTextXML",
    "30693443": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6614136/fullTextXML",
    "27284927": "https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0157197",
    "39845481": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11751902/fullTextXML",
}
ABSTRACT_LIVE = {"22886373": "https://link.springer.com/article/10.1007/s11010-012-1423-7"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    blob = args.input.read_bytes()
    if hashlib.sha256(blob).hexdigest() != REVIEWED_CORPUS_SHA256:
        raise ValueError("Input differs from the manually reviewed snapshot; fresh screening cannot inherit its read ledger")
    corpus = json.loads(blob)
    records = corpus["results"]
    if len({(r["source"], str(r["id"])) for r in records}) != len(records):
        raise ValueError("Duplicate source/record identity")
    if not set(READ_LEDGER) <= {str(r["id"]) for r in records if r["source"] == "MED"}:
        raise ValueError("Manual read ledger contains absent records")
    benefit, intervention, sensitivity = [re.compile(p, re.I) for p in (BENEFIT, INTERVENTION, SENSITIVITY_TITLE)]
    rows = []
    for rec in records:
        rid, source = str(rec["id"]), rec["source"]
        title, abstract = rec.get("title", ""), rec.get("abstractText", "")
        combined = title + " " + abstract
        b = bool(benefit.search(combined))
        p = bool(intervention.search(combined))
        s = not b and bool(sensitivity.search(title))
        missing = not abstract.strip()
        read = source == "MED" and rid in READ_LEDGER
        verified_fulltext = source == "MED" and rid in FULLTEXT_PASSAGES
        verified_abstract = source == "MED" and rid in ABSTRACT_LIVE
        priority = "benefit_and_intervention" if b and p else "benefit_only" if b else "sensitivity_title" if s else "not_prioritized_by_filters"
        rows.append({
            "source": source, "record_id": rid, "pmid": rec.get("pmid", ""),
            "title": title, "abstract_available": int(not missing),
            "publication_types": " | ".join(rec.get("pubTypeList", {}).get("pubType", [])),
            "machine_priority": priority, "benefit_match": int(b),
            "generic_intervention_match": int(p), "sensitivity_title_match": int(s),
            "selected_title_screen": int(b or s),
            "manual_title_seen": int(b or s or missing),
            "manual_full_abstract_read": int(read),
            "verification_level": "primary_methods_results_passages" if verified_fulltext else "primary_abstract_live" if verified_abstract else "saved_abstract_only" if read else "not_primary_adjudicated",
            "primary_source_url": FULLTEXT_PASSAGES.get(rid, ABSTRACT_LIVE.get(rid, "")) if source == "MED" else "",
            "manual_appraisal": READ_LEDGER.get(rid, "") if source == "MED" else "",
        })
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "gene-record-screen.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    selected = [{"source": r["source"], "id": r["record_id"]} for r in rows if r["selected_title_screen"]]
    summary = {
        "date": "2026-09-21", "source_file": args.input.name,
        "source_sha256": hashlib.sha256(blob).hexdigest(), "query": corpus["query"],
        "filters": {"benefit": BENEFIT, "generic_intervention": INTERVENTION, "sensitivity_title_nonbenefit_only": SENSITIVITY_TITLE},
        "records": len(rows), "unique_source_ids": len({(r["source"], r["record_id"]) for r in rows}),
        "abstracts_available": sum(r["abstract_available"] for r in rows),
        "abstracts_missing": sum(not r["abstract_available"] for r in rows),
        "selected_title_records": len(selected),
        "manual_title_seen_including_all_missing_abstract_titles": sum(r["manual_title_seen"] for r in rows),
        "manual_full_abstracts_read": sum(r["manual_full_abstract_read"] for r in rows),
        "primary_methods_results_passages_verified": len(FULLTEXT_PASSAGES),
        "additional_primary_abstract_live_verified": len(ABSTRACT_LIVE),
        "machine_priority_counts": dict(Counter(r["machine_priority"] for r in rows)),
        "source_counts": dict(Counter(r["source"] for r in rows)),
        "review_tagged_records": sum(any("review" in t.lower() for t in rec.get("pubTypeList", {}).get("pubType", [])) for rec in records),
        "limitations": "Manual title screen is not manual abstract adjudication. Twenty-one abstracts and seven primary methods/results passage sets plus one live primary abstract were read. Remaining records are unadjudicated, not negative. Corpus can omit full-text-only interventions, missing abstracts, other genes, unpublished/non-indexed work. No drug catalog gate and no clinical efficacy imputation.",
    }
    (args.out_dir / "gene-record-screen-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
