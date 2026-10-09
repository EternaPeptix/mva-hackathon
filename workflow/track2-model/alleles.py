#!/usr/bin/env python3
"""Reconstruct BUB1B canonical allele products from pinned public GRCh38 inputs.

This is sequence reconstruction, not a protein-function or treatment simulator.
Requires only Python 3.9+ and an uncompressed, already indexed reference FASTA.
The source directory is the retained public Ensembl transcript input directory.
No alignment, patient RNA, case VCF, or other raw case input is read.
"""

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time


TRANSCRIPT = "ENST00000287598"
VERSION = 11
REFERENCE_SHA256 = "187917c00fa7373a054f9768009d19e4f6e5fc87d3a6c125ce131164ca69fcd6"
SOURCES = {
    "exons.json": (
        "dec2ff9c23822d4ac3259f5153ea1d548b8ba264c4eddfdc49ab19b91bde2d80",
        "https://rest.ensembl.org/overlap/id/ENST00000287598?feature=exon;content-type=application/json",
    ),
    "cds.json": (
        "573be6f7288f6b59ff036cf31c36b17169152d196fcab5ffab137828fd79a8e1",
        "https://rest.ensembl.org/sequence/id/ENST00000287598?type=cds",
    ),
    "info.json": (
        "86e87c46ae1a424702676493d024d0508ce60f6c951b0b2803ce10269e09c20c",
        "https://rest.ensembl.org/info/data",
    ),
    "map-stop.json": (
        "3aabdfb29f2e443501e114daa50aec9bbfcca6f1d1b4e7041bfa0bc128df6b4b",
        "https://rest.ensembl.org/map/cds/ENST00000287598/2210..2210?content-type=application/json",
    ),
}
VARIANTS = {
    "L737Ter": {"position": 40209701, "coding_position": 2210,
                "protein_position": 737, "reference_codon": "TTA", "alternate_codon": "TGA"},
    "N1002Lys": {"position": 40220612, "coding_position": 3006,
                 "protein_position": 1002, "reference_codon": "AAT", "alternate_codon": "AAG"},
}
# Standard nuclear genetic code, codons ordered T/C/A/G at each position.
CODONS = dict(zip(
    ("".join(c) for c in itertools.product("TCAG", repeat=3)),
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sequence_sha256(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


class IndexedFasta:
    """Read an existing FASTA index; never create or modify reference files."""

    def __init__(self, path):
        self.path = path
        self.index_path = Path(str(path) + ".fai")
        self.index = {}
        for line in self.index_path.read_text().splitlines():
            fields = line.split("\t")
            require(len(fields) >= 5 and fields[0] not in self.index,
                    "malformed or duplicate FASTA index entry")
            self.index[fields[0]] = tuple(int(value) for value in fields[1:5])

    def fetch(self, contig, start, end):
        """Fetch zero-based, half-open genomic coordinates."""
        length, offset, bases, width = self.index[contig]
        require(0 <= start < end <= length and 0 < bases <= width,
                "invalid FASTA interval or wrapping metadata")
        first = offset + (start // bases) * width + start % bases
        last_base = end - 1
        last = offset + (last_base // bases) * width + last_base % bases
        with self.path.open("rb") as stream:
            stream.seek(first)
            sequence = stream.read(last - first + 1).replace(b"\n", b"").replace(b"\r", b"")
        require(len(sequence) == end - start, "FASTA index or sequence-length mismatch")
        return sequence.decode("ascii").upper()


def translate(cds):
    require(len(cds) % 3 == 0 and cds.startswith("ATG"), "unexpected complete CDS boundaries")
    protein = []
    for offset in range(0, len(cds), 3):
        codon = cds[offset:offset + 3]
        require(codon in CODONS, "noncanonical CDS base")
        amino_acid = CODONS[codon]
        if amino_acid == "*":
            return "".join(protein), offset // 3 + 1, codon
        protein.append(amino_acid)
    raise ValueError("CDS has no in-frame termination codon")


def reconstruct(source_directory, reference):
    source_data, source_identity = {}, {}
    for name, (expected_sha, url) in SOURCES.items():
        path = source_directory / name
        observed_sha = sha256(path)
        require(observed_sha == expected_sha, "public source identity changed: " + name)
        source_data[name] = json.loads(path.read_text())
        source_identity[name] = {"sha256": observed_sha, "bytes": path.stat().st_size, "url": url}

    cds_data = source_data["cds.json"]
    require(cds_data["id"] == TRANSCRIPT and cds_data["version"] == VERSION
            and cds_data["molecule"] == "dna", "unexpected transcript/CDS identity")
    require(source_data["info.json"]["releases"] == [116], "unexpected Ensembl release")
    mapping = source_data["map-stop.json"]["mappings"]
    require(len(mapping) == 1 and mapping[0]["start"] == mapping[0]["end"] == 40209701
            and mapping[0]["strand"] == 1 and mapping[0]["assembly_name"] == "GRCh38"
            and mapping[0]["seq_region_name"] == "15", "public stop-coordinate map disagrees")
    exons = sorted((e for e in source_data["exons.json"] if e["Parent"] == TRANSCRIPT),
                   key=lambda e: e["rank"])
    require(len(exons) == 23 and [e["rank"] for e in exons] == list(range(1, 24)),
            "missing, duplicate or unexpected transcript exons")
    require(all(e["strand"] == 1 and e["seq_region_name"] == "15"
                and e["assembly_name"] == "GRCh38" and e["start"] <= e["end"] for e in exons),
            "this reconstruction requires positive-strand GRCh38 chromosome 15")
    require(all(a["end"] < b["start"] for a, b in zip(exons, exons[1:])),
            "overlapping or unsorted exons")

    reference_sha = sha256(reference)
    require(reference_sha == REFERENCE_SHA256, "reference FASTA identity changed")
    fasta = IndexedFasta(reference)
    aliases = set(fasta.index) & {"15", "chr15"}
    require(len(aliases) == 1, "ambiguous or absent chromosome 15 alias")
    contig = next(iter(aliases))
    spliced = "".join(fasta.fetch(contig, e["start"] - 1, e["end"]) for e in exons)
    cds = cds_data["seq"].upper()
    require(len(cds) == 3153 and len(spliced) == 3669 and spliced.count(cds) == 1,
            "public CDS does not match the reconstructed reference transcript exactly once")
    cds_offset = spliced.index(cds)
    require(cds_offset == 152, "unexpected CDS start in spliced transcript")
    final_junction = sum(e["end"] - e["start"] + 1 for e in exons[:-1])
    findings = {}
    for name, variant in VARIANTS.items():
        preceding = 0
        for exon in exons:
            if exon["start"] <= variant["position"] <= exon["end"]:
                transcript_position = preceding + variant["position"] - exon["start"] + 1
                break
            preceding += exon["end"] - exon["start"] + 1
        else:
            raise ValueError("candidate outside transcript: " + name)
        c = transcript_position - cds_offset
        require(c == variant["coding_position"] and cds[c - 1] == "T"
                and fasta.fetch(contig, variant["position"] - 1, variant["position"]) == "T",
                "candidate coordinate/reference mismatch: " + name)
        start = (c - 1) // 3 * 3
        edited_codon = list(cds[start:start + 3])
        edited_codon[(c - 1) % 3] = "G"
        require(cds[start:start + 3] == variant["reference_codon"]
                and "".join(edited_codon) == variant["alternate_codon"],
                "candidate codon context mismatch: " + name)
        findings[name] = dict(variant, chrom="chr15", ref="T", alt="G",
                             transcript_position=transcript_position, exon_rank=exon["rank"])

    products = {}
    definitions = {"WT": [], "L737Ter": ["L737Ter"], "N1002Lys": ["N1002Lys"],
                   "L737Ter_N1002Lys_cis": ["L737Ter", "N1002Lys"]}
    for name, changes in definitions.items():
        edited = list(cds)
        for change in changes:
            edited[VARIANTS[change]["coding_position"] - 1] = "G"
        edited_cds = "".join(edited)
        protein, first_stop, stop_codon = translate(edited_cds)
        stop_end = cds_offset + first_stop * 3
        premature = first_stop < len(cds) // 3
        distance = final_junction - stop_end
        products[name] = {
            "variants": changes,
            "edited_cds_sha256": sequence_sha256(edited_cds),
            "protein_sequence": protein,
            "protein_sha256": sequence_sha256(protein),
            "protein_length_aa": len(protein),
            "first_stop_codon_aa_position": first_stop,
            "first_stop_codon_dna": stop_codon,
            "premature_stop": premature,
            "masked_downstream_variants": [v for v in changes
                                           if VARIANTS[v]["protein_position"] > first_stop],
            "positional_nmd": {
                "rule": "premature termination more than 50 nt upstream of final exon junction",
                "stop_codon_end_transcript_coordinate": stop_end,
                "final_junction_after_transcript_coordinate": final_junction,
                "spliced_bases_from_stop_end_to_final_junction": distance,
                "rule_positive": premature and distance > 50,
                "rna_survival_fraction": None,
                "interpretation": "positional prediction only; actual RNA decay is unmeasured",
            },
            "endogenous_product_detected": None,
            "relative_abundance": None,
            "functional_activity": None,
        }
    wt, trunc, miss, cis = (products[n] for n in definitions)
    require(wt["protein_length_aa"] == miss["protein_length_aa"] == 1050
            and trunc["protein_length_aa"] == cis["protein_length_aa"] == 736,
            "unexpected reconstructed protein lengths")
    require(trunc["protein_sequence"] == cis["protein_sequence"] == wt["protein_sequence"][:736],
            "stop-only/cis truncation invariant failed")
    differences = [i + 1 for i, (a, b) in enumerate(zip(wt["protein_sequence"], miss["protein_sequence"]))
                   if a != b]
    require(differences == [1002] and wt["protein_sequence"][1001] == "N"
            and miss["protein_sequence"][1001] == "K", "missense translation invariant failed")
    require(trunc["positional_nmd"]["spliced_bases_from_stop_end_to_final_junction"] == 746
            and cis["masked_downstream_variants"] == ["N1002Lys"], "phase/NMD invariant failed")
    return {
        "schema_version": 1,
        "scope": "Reference-based canonical sequence reconstruction. No patient RNA/protein, "
                 "folding, functional activity, quantitative NMD, drug response or clinical benefit is inferred.",
        "translation_assumptions": ["Standard nuclear code; initiation at the annotated CDS start",
                                    "Termination at the first in-frame stop; no readthrough, reinitiation or alternative splicing",
                                    "Products describe sequence if translated, not proof they exist endogenously"],
        "transcript": {"id": TRANSCRIPT + "." + str(VERSION), "assembly": "GRCh38", "strand": 1,
                       "ensembl_release": 116, "exon_count": len(exons), "spliced_length_nt": len(spliced),
                       "coding_length_nt_including_terminal_stop": len(cds), "cds_start_transcript_coordinate": cds_offset + 1,
                       "spliced_reference_sha256": sequence_sha256(spliced), "cds_sha256": sequence_sha256(cds),
                       "full_cds_matches_spliced_reference": True},
        "reference": {"sha256": reference_sha, "bytes": reference.stat().st_size,
                      "index_sha256": sha256(fasta.index_path), "contig": contig,
                      "identity_verified_by": "full-file SHA256 during this run"},
        "public_source_identity": source_identity,
        "variants": findings,
        "allele_products": products,
        "observed_phase": "unresolved",
        "phase_scenarios": {
            "reference": {"alleles": ["WT", "WT"]},
            "trans": {"alleles": ["L737Ter", "N1002Lys"]},
            "cis": {"alleles": ["L737Ter_N1002Lys_cis", "WT"],
                    "interpretation": "With only these two modeled variants, the second allele has WT coding sequence. "
                                      "Its expression, dosage compensation and functional sufficiency are not established."},
        },
        "phase_scenario_limit": "These are conditional two-variant configurations, not measured phase or an assertion "
                                "that other sequence, regulatory or structural changes are absent.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-directory", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), "output already exists; choose a new path")
    started = time.monotonic()
    result = reconstruct(args.source_directory, args.reference)
    result["execution"] = {"script_sha256": sha256(Path(__file__)),
                           "elapsed_seconds": round(time.monotonic() - started, 3),
                           "completed_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"output": str(args.output), "transcript": result["transcript"]["id"],
                      "protein_lengths_aa": {key: value["protein_length_aa"]
                                             for key, value in result["allele_products"].items()},
                      "phase": result["observed_phase"], "elapsed_seconds": result["execution"]["elapsed_seconds"]}))


if __name__ == "__main__":
    main()
