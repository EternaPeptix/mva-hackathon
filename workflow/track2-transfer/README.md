# Public transfer demonstration

**The workflow transfers; the drug nomination does not automatically transfer.**
We replay the existing BUB1B hypothesis and a synthetic homozygous TRIP13
Arg354Ter genotype using public reference sequences and curated primary evidence.
TRIP13 reaches a different decision: the BUBR1-stabilization route is not
supported, and no qualified drug lead emerges from the three reviewed routes.
This is a completed second worked example, **not validation in two patients**.

From a fresh clone, Python 3.9+ and its standard library are sufficient:

```sh
python3 workflow/track2-transfer/replay.py
```

No network, controlled data, package installation or genomic reference download
is needed. Read the generated [comparison table and gaps](results/transfer.md)
or [machine-readable result](results/transfer.json).

| Step | Recomputed or curated? | What is actually checked |
| --- | --- | --- |
| Selected coding variants → sequence | Computed | Apply coding SNVs to versioned RefSeq CDS; reject reference mismatches. This begins after variant selection, not at raw reads/VCF. |
| Sequence → consequence | Computed | Translate canonical products; report changed codon, amino acid, length and positional NMD rule. RefSeq protein annotation and full Ensembl CDS must agree. |
| Consequence → mechanism | Curated | Explicit source-backed assertions in `cases.json`; sequence alone cannot establish protein stability, checkpoint function or drug response. |
| Mechanism → route decision | Computed from curated assertions | Apply the same route gates to both examples. A failed gate prevents transfer; missing evidence remains missing. No clinical recommendation is generated. |

The synthetic second input is **NM_004237.4(TRIP13):c.1060C>T on both alleles**,
drawn from [Yost 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5493194/)
and [ClinVar VCV000431045.1](https://www.ncbi.nlm.nih.gov/clinvar/variation/431045/).
It reproduces CGA→TGA, p.(Arg354Ter), with 353 residues if translated
(reference: 432). The premature stop ends 141 nucleotides before the final
exon junction. The positional NMD rule agrees with the published RNA result;
this replay itself does not measure RNA. ClinVar is one literature-only
OMIM submission with no assertion criteria, not an independent expert review.
No individual from the paper is represented or given an invented phenotype.

The mechanistic distinction matters. The published TRIP13 experiments support
loss of its function and rescue by WT TRIP13. They do not show that increasing
BUBR1 substitutes for it. Thus a common syndrome/checkpoint phenotype is
insufficient to transfer niacin's proposed BUBR1-restoration rationale.
This is not evidence that every possible BUBR1 bypass is ineffective.
Arimoclomol is not transferred as a full-length missense-stabilization route;
everolimus lacks the required TRIP13-specific causal signaling evidence.

To check sources live without replacing the pinned cache:

```sh
python3 workflow/track2-transfer/replay.py --refresh
```

This compares current records with [source identities](sources.json), ignoring
only JSON key serialization order; cached byte identities are checked exactly.
A changed upstream record stops the refresh and leaves the cache intact; it
requires deliberate review/re-pinning. The five cached public records total
58,465 bytes. They contain reference sequences/annotations and the TRIP13
ClinVar record, with no case reads or case input hashes.

The small verification checks compare both reference providers, the published
TRIP13 consequence, reference-base rejection and behavior when a required
mechanistic assertion is removed:

```sh
python3 workflow/track2-transfer/test_replay.py
```

Add an example by supplying its versioned CDS records, coding SNVs and sourced
mechanism assertions in [cases.json](cases.json); the engine has no gene-name
branches. The source cache list is separate in `sources.json`. `--config` and
`--out` accept alternate locations. Current support is deliberately limited to
complete contiguous CDS annotations and nuclear-code coding SNVs. It does not
discover candidates, interpret CNVs/splice variants, resolve phase, infer
pathogenicity, rank all medicines or replace expert literature review.
The principal scaling bottleneck remains gene/allele-specific evidence
curation and experimental qualification. Two executions do not estimate
generalization accuracy, sensitivity or clinical benefit.
