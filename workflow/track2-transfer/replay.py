#!/usr/bin/env python3
"""Public-reference SNV consequences followed by separately curated route gates.

Python standard library only. No patient reads, VCF, phenotype file or model
fitting is used. Supports complete, contiguous RefSeq CDSs and coding SNVs;
unsupported inputs fail rather than silently receive a consequence.
"""

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import re
import sys
import urllib.request
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
CODE = dict(zip(("".join(x) for x in itertools.product("TCAG", repeat=3)),
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
AA3 = dict(zip("ACDEFGHIKLMNPQRSTVWY*", ["Ala", "Cys", "Asp", "Glu", "Phe",
               "Gly", "His", "Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln",
               "Arg", "Ser", "Thr", "Val", "Trp", "Tyr", "Ter"]))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_sources(manifest, refresh=False):
    loaded = {}
    for name, source in manifest.items():
        cached = (ROOT / "cache" / name).read_bytes()
        require(digest(cached) == source["sha256"], f"Cache identity mismatch: {name}")
        if refresh:
            with urllib.request.urlopen(source["url"], timeout=45) as response:
                fresh = response.read()
            # A refresh verifies live sources; it never replaces the pinned record.
            # Ensembl randomizes JSON key order between otherwise identical
            # responses. Ignore serialization order, never biological content.
            unchanged = (json.loads(fresh) == json.loads(cached)) if name.endswith(".json") else fresh == cached
            require(unchanged,
                    f"Live record changed: {name}. Inspect/re-pin deliberately; cache preserved.")
        loaded[name] = cached.decode("utf-8")
    return loaded


def refseq(record):
    accession = re.search(r"^VERSION\s+(\S+)", record, re.M).group(1)
    matches = re.findall(r"^     CDS             (\d+)\.\.(\d+)$", record, re.M)
    require(len(matches) == 1, f"Unsupported CDS location in {accession}")
    start, end = map(int, matches[0])
    origin = record.split("\nORIGIN", 1)[1].split("//", 1)[0]
    sequence = re.sub(r"[\s\d]", "", origin).upper()
    require(set(sequence) <= set("ACGT"), "Ambiguous reference sequence")
    cds = sequence[start - 1:end]
    require(len(cds) == end - start + 1, "CDS outside sequence")
    exons = [tuple(map(int, x)) for x in re.findall(
        r"^     exon            (\d+)\.\.(\d+)$", record, re.M)]
    require(exons and exons[0][0] == 1 and exons[-1][1] == len(sequence),
            "Incomplete exon annotation")
    require(all(a[1] + 1 == b[0] for a, b in zip(exons, exons[1:])),
            "Non-contiguous transcript exons")
    annotated = re.search(r'/translation="([A-Z\s]+)"', record).group(1)
    annotated = re.sub(r"\s", "", annotated)
    require(translate(cds) == annotated, "Computed translation differs from RefSeq annotation")
    return {"accession": accession, "cds": cds, "cds_start": start,
            "exons": exons, "reference_protein": annotated}


def translate(cds):
    require(len(cds) % 3 == 0 and cds.startswith("ATG"), "Unsupported CDS boundaries")
    protein = []
    for i in range(0, len(cds), 3):
        aa = CODE[cds[i:i + 3]]
        if aa == "*":
            return "".join(protein)
        protein.append(aa)
    raise ValueError("Missing terminal stop")


def apply_variant(reference, variant):
    cds = reference["cds"]
    match = re.fullmatch(r"c\.(\d+)([ACGT])>([ACGT])", variant)
    require(match is not None, f"Only single coding substitutions supported: {variant}")
    pos, old, new = int(match[1]), match[2], match[3]
    require(1 <= pos <= len(cds) and cds[pos - 1] == old, f"Reference mismatch: {variant}")
    require(old != new, "Identity substitution is not a variant")
    edited = cds[:pos - 1] + new + cds[pos:]
    codon_start = ((pos - 1) // 3) * 3
    codon, alt_codon = cds[codon_start:codon_start + 3], edited[codon_start:codon_start + 3]
    aa, alt_aa = CODE[codon], CODE[alt_codon]
    aa_pos = codon_start // 3 + 1
    require(aa != "*", "Stop-loss variants are outside this demo")
    consequence = "stop_gained" if alt_aa == "*" else "synonymous" if aa == alt_aa else "missense"
    protein = translate(edited)
    stop_end = reference["cds_start"] - 1 + (len(protein) + 1) * 3
    distance = reference["exons"][-2][1] - stop_end if len(reference["exons"]) > 1 else None
    return {"coding_hgvs": variant, "protein_hgvs": f"p.({AA3[aa]}{aa_pos}{AA3[alt_aa]})",
            "reference_codon": codon, "alternate_codon": alt_codon,
            "consequence": consequence, "protein_length_if_translated": len(protein),
            "full_length_sequence_possible": len(protein) == len(reference["reference_protein"]),
            "stop_end_to_final_junction_nt": distance,
            "positional_nmd_prediction": consequence == "stop_gained" and distance is not None and distance > 50,
            "nmd_limit": "Positional >50 nt rule only; not a measurement of RNA decay."}


def evaluate_route(route, mechanism, consequences, phase):
    evidence = mechanism["route_evidence"].get(route["id"], {})
    gates = []
    for gate in route["gates"]:
        if gate["kind"] == "sequence":
            state = any(x["full_length_sequence_possible"] for x in consequences.values())
            basis = "Computed canonical products; protein abundance and activity remain unmeasured."
        else:
            observed = evidence.get(gate["key"], {})
            state = observed.get("state")
            basis = observed.get("basis", "No supporting assertion in this bounded evidence set.")
        gates.append({"gate": gate["key"], "state": state, "basis": basis})
    states = [x["state"] for x in gates]
    if False in states:
        status = "not_transferred"
    elif None in states:
        status = "evidence_gap"
    else:
        status = "conditional_research_route"
    return {"route": route["name"], "status": status, "gates": gates,
            "clinical_recommendation": False,
            "phase_constraint": "Resolve trans; confirmed cis retires the proposed pair." if phase == "unresolved" else "Homozygosity is stipulated for the synthetic example, not observed here.",
            "required_next_evidence": evidence.get("next", "Qualify the gene-specific mechanism and drug bridge before nomination.")}


def run(config, sources):
    results = []
    for case in config["cases"]:
        ref = refseq(sources[case["refseq_file"]])
        require(ref["accession"] == case["transcript"], "Transcript version mismatch")
        ensembl = json.loads(sources[case["ensembl_file"]])
        require(f'{ensembl["id"]}.{ensembl["version"]}' == case["ensembl_transcript"],
                "Ensembl transcript version changed")
        require(ensembl["seq"] == ref["cds"], "RefSeq and Ensembl CDS disagree")
        variants = {x: apply_variant(ref, x) for x in case["variants"]}
        for variant, consequence in variants.items():
            expected = case["external_consequence_check"][variant]
            require(consequence["protein_hgvs"] == expected, "Sequence differs from external consequence assertion")
        if case.get("clinvar_file"):
            archive = ET.fromstring(sources[case["clinvar_file"]]).find("VariationArchive")
            require(archive is not None and archive.attrib["VariationID"] == case["clinvar_variation_id"],
                    "ClinVar identity mismatch")
            require(all(x in archive.attrib["VariationName"] for x in [case["transcript"], *case["variants"]]),
                    "ClinVar variant/transcript disagreement")
        mechanism = config["mechanisms"][case["mechanism_id"]]
        routes = [evaluate_route(x, mechanism, variants, case["phase"]) for x in config["routes"]]
        results.append({"id": case["id"], "input_type": case["input_type"], "gene": case["gene"],
            "genotype": case["genotype"], "transcript": ref["accession"],
            "cds_length_nt": len(ref["cds"]), "reference_protein_length_aa": len(ref["reference_protein"]),
            "sequence_checks": "RefSeq translation and complete Ensembl CDS agree",
            "variants": variants, "mechanism": mechanism["summary"],
            "mechanism_evidence_type": "Human-curated primary-source assertions, not inferred from sequence",
            "source_ids": mechanism["source_ids"], "routes": routes,
            "clinical_lead": None})
    return {"scope": config["scope"], "cases": results}


def markdown(result, config):
    lines = ["# Public transfer result", "", result["scope"], "",
             "| Example | Input and computed consequence | Curated mechanism | Route result |",
             "| --- | --- | --- | --- |"]
    for case in result["cases"]:
        seq = "; ".join(f'{v["coding_hgvs"]}: {v["reference_codon"]}→{v["alternate_codon"]}, {v["protein_hgvs"]}, {v["protein_length_if_translated"]} aa if translated' for v in case["variants"].values())
        decisions = "; ".join(f'{r["route"]}: **{r["status"]}**' for r in case["routes"])
        lines.append(f'| {case["gene"]}: {case["input_type"]} | {case["genotype"]}; {seq} | {case["mechanism"]} | {decisions} |')
    for case in result["cases"]:
        lines += ["", f'## {case["gene"]}: evidence gaps', "",
                  f'{case["transcript"]}: {case["cds_length_nt"]} nt CDS including termination; {case["reference_protein_length_aa"]} reference amino acids. {case["sequence_checks"]}.', ""]
        for r in case["routes"]:
            lines += [f'**{r["route"]} — {r["status"]}.** {r["required_next_evidence"]}', ""]
            for g in r["gates"]:
                lines.append(f'- {g["gate"]}: {str(g["state"]).lower()} — {g["basis"]}')
            lines += [""]
        lines += [r["phase_constraint"], "", "No qualified clinical drug lead is produced."]
    lines += ["", "## Evidence sources", ""]
    lines += [f'- **{key}:** [{s["title"]}]({s["url"]}) — {s["scope"]}' for key, s in config["citations"].items()]
    lines += ["", "Only sequence consequences and gate execution were recomputed. This is two worked examples, including one synthetic public genotype; it is not two-patient biological validation, a drug efficacy benchmark or an automated literature review.", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true", help="Verify live source content against pinned cache; never overwrite cache")
    parser.add_argument("--config", type=Path, default=ROOT / "cases.json")
    parser.add_argument("--out", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    manifest = json.loads((ROOT / "sources.json").read_text())
    sources = read_sources(manifest, args.refresh)
    result = run(config, sources)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "transfer.json").write_text(json.dumps(result, indent=2) + "\n")
    (args.out / "transfer.md").write_text(markdown(result, config))
    print(f'Completed {len(result["cases"])} public-reference examples; results: {args.out}')
    for case in result["cases"]:
        print(case["gene"] + ": " + "; ".join(r["route"] + "=" + r["status"] for r in case["routes"]))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, OSError) as error:
        sys.exit(f"Transfer replay failed: {error}")
