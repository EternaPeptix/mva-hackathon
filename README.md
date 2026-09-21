# MVA: from genome evidence to a testable rescue experiment

Team **jeweled** focuses on the case's [recorded muscle atrophy and failure to thrive](workflow/configs/lirical-subject.json):
could a gene-informed medication hypothesis preserve muscle function and the cells that renew it?

- **[Track 1: candidate and evidence](reports/jeweled_track1_report_submission.md).**
  BUB1B **p.(Leu737Ter) + p.(Asn1002Lys)**, at subjective **EPCR 0.8**.
  Phase is unresolved and the missense remains VUS.
- **[Track 2: selected medications and proposal](reports/jeweled_track2_report_submission.md).**
  **Everolimus leads; niacin (nicotinic acid) is the main comparator**, tested separately.
  They test downstream mTORC1 normalization versus usable BUBR1 restoration, with the
  same goal: greater absolute muscle output, preserved renewal and no worsened chromosome fidelity.
- **[Submission index, methods and proposed pitch](reports/README.md).**
  **Track 1 model 1 scored 100 rank points and F-max 1.000 on 21 September.**
  [Submission record](reports/README.md#track-1-submission-record) distinguishes the
  confirmed score from the pending final report-persistence check. Track 2 submission
  and media/Vimeo remain on hold; no drug experiment has run.

## Start here

- [Reproduction guide and public replay](REPRODUCIBILITY.md)
- [Recorded environment](workflow/environments/20260919/runtime.json)
- [Benchmark methods and results](workflow/docs/PUBLIC-BENCHMARK-METHOD.md)
- [Evidence and component index](results/INDEX.md)

**Native** means an installed bioinformatics executable was run. **Synthetic controls**
(“fictional” in old files) use invented inputs with expected outcomes. A **receipt** records
commands, versions, input/output identities and results; it is not itself a biological finding.

## Pipeline

**Track 1 — rank a causal hypothesis**

1. Verify supplied reads; align four sequencing lanes to GRCh38.
2. Run GATK and bcftools; preserve original genotypes through normalization and annotation.
3. Review phenotype, population frequency, recessive models and alternative explanations.
4. Examine read support, phase, structural evidence and coverage limitations.
5. Rank the retained hypothesis with an explicit estimated probability of causality (EPCR).

**Track 2 — turn the mechanism into an experiment**

6. Qualify allele products, useful BUBR1 function and harmful mTORC1 excess in controlled models.
7. Compare everolimus with niacin separately at justified exposure, with genetic and tumor controls.
8. Measure absolute muscle output, renewal and fidelity; qualify arimoclomol and senolysis reserves.

## Quick start: public HG002 chromosome 20

```sh
bash workflow/public-benchmark/run.sh /absolute/path/to/new-output
```

The command downloads public inputs and runs both callers plus RTG vcfeval. See
[requirements and actual platform verification](REPRODUCIBILITY.md#public-benchmark-quick-start).
The full Mac Studio replay reproduced the primary counts in **134.3 minutes**, including
downloads and setup; controlled case data are separate.

## Validation summary

HG002 reads are labelled nominally 35×; measured mean depth is **36.97× inside the confident
regions**, with **99.25% of bases at least 20×** (MAPQ≥20/BQ≥20). Comparison uses GIAB v4.2.1,
its official chr20 confident BED and genotype-aware RTG Tools 3.13 vcfeval, with full-chromosome
context and raw calls. [Methods](workflow/docs/PUBLIC-BENCHMARK-METHOD.md) explain the definitions.

| Public check | Recorded result |
|---|---|
| GATK SNVs / indels, original split calls | Recall 99.3748% / 97.3902%; precision 98.9548% / 96.0168% |
| bcftools SNVs / indels, original split calls | Recall 99.3593% / 91.5115%; precision 98.9512% / 91.1827% |
| DELLY structural variants | 69/802 truth alleles matched (**8.6%**) |
| SvABA structural variants | 171/802 truth alleles matched (**21.3%**) |

Retaining native multiallelic genotypes improves the indel comparison; this later sensitivity
analysis is reported separately. [Structural results and miss analysis](results/INDEX.md#structural-variants)
retain the same truth denominator, including difficult insertions and repeat-associated alleles.

Three useful implementation findings changed the workflow: the pinned GATK reader rejects
CRAM 3.1; samtools' checksum misses a FLAG-preservation difference; and GATK uses the maximum
ploidy spanning an active region. Methods and control counts are in the component index.

## Findings and repurposing proposal

**Track 1:** BUB1B **p.(Leu737Ter) + p.(Asn1002Lys)** is our leading hypothesis, at subjective
**EPCR 0.8**. Both alleles are observed; trans is favored under the recessive model but unobserved.
The stop has a two-star ClinVar P/LP assertion and the missense remains VUS. The
[report](reports/jeweled_track1_report_submission.md) explains the two-allele census, ten-gene
differential, 288 bp phase blocks versus a 10,911 bp gap, and retracted mosaic signal.

**Track 2:** nominate **everolimus** as the lead drug hypothesis and **niacin (nicotinic acid)**
as the main comparator, tested separately. They address mTORC1 normalization versus restoration
of usable BUBR1. Everolimus has the clearer drug-to-target link; niacin connects human muscle
NAD engagement to a separate BUBR1-stability mechanism. Neither has demonstrated rescue of
this genotype. Arimoclomol and dasatinib remain conditional reserves. Our
[RNA reanalysis](workflow/track2-rna/README.md) adds no mTORC1 support and does not rebut the
published protein-activity result. The decision rests on literature and testability, not a
simulated drug ranking or clinical-benefit claim.

The [proposal](reports/jeweled_track2_report_submission.md) and
[evidence pipeline](reports/appendices/track2-evidence-pipeline.md) explain the
qualitative selection, exposure gate, absolute-output endpoints and rhabdomyosarcoma counter-screen.
The [unelected-route table](reports/jeweled_track2_report_submission.md#other-routes-considered-and-not-selected-as-the-lead) explains alternatives and weaknesses.
The [genome-informed model](workflow/track2-model/README.md) separates sequence reconstruction
from assumption-based population scenarios; it does not predict drug response. No drug experiment
has run, and neither chromosome correction nor clinical efficacy is claimed.

## Scope and limitations

**Low SV recall (8.6% DELLY; 21.3% SvABA) limits the causal search:** a negative short-read result
cannot exclude difficult structural alternatives. Annotation locates misses without proving their
cause. One public chromosome and synthetic tests do not establish case-wide clinical sensitivity.
Phase, endogenous allele function, rescue competence and drug benefit remain unresolved. EPCR is a research judgment.
Exploratory control development closed on 20 September; remaining work concerns reproduction,
documented defects and submission delivery before **24 October 2026**.

## Data handling, disclosure and license

The repository remains private; curated findings are included, raw reads and genome-scale case
payloads are excluded. [Release review](workflow/docs/PUBLIC-RELEASE-AUDIT.md) records decisions.
Claude Code/Codex used personal subscriptions with training disabled (user-confirmed). Narration
uses Cartesia Sonic 3.6/Lauren; Cartesia account settings remain unverified. See [LICENSE](LICENSE),
[required attribution and dataset citation status](reports/PROVENANCE.md#acknowledgement), and [CHANGELOG](CHANGELOG.md).
