# MVA: from genome evidence to a testable rescue experiment

Team **jeweled** asks whether a recessive BUB1B mechanism explains the supplied MVA case,
and whether restoring usable BUBR1 protein could improve mitotic fidelity. Track 1's CSV
and both research reports are prepared; submission is pending and Vimeo is on hold.

## Start here

- [Reproduction guide and public replay](REPRODUCIBILITY.md)
- [Recorded environment](workflow/environments/20260919/runtime.json)
- [Benchmark methods and results](workflow/docs/PUBLIC-BENCHMARK-METHOD.md)
- [Current reports, prediction, methods form and pitch](reports/README.md)
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

6. Test whether restoring residual BUBR1 abundance restores mitotic function.
7. Rank drug hypotheses; specify genetic controls, exposure and tumor counter-screens.
8. Track first-division errors and yield from each starting cell; define advancement criteria.

## Quick start: public HG002 chromosome 20

```sh
bash workflow/public-benchmark/run.sh /absolute/path/to/new-output
```

The command downloads public inputs and runs both callers plus RTG vcfeval. See
[requirements and actual platform verification](REPRODUCIBILITY.md#public-benchmark-quick-start).
The earlier full caller run took about 116 minutes; controlled case data are separate.

## Validation summary

HG002 reads are labelled nominally 35×. Comparison uses GIAB v4.2.1, its official chr20
confident BED and genotype-aware RTG Tools 3.13 vcfeval, with full-chromosome context and
raw calls. [Methods](workflow/docs/PUBLIC-BENCHMARK-METHOD.md) explain denominators and thresholds.

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

**Track 2:** predicted stop-allele loss plus a potentially rescuable missense product motivates
restoring **usable BUBR1 abundance**. The first experiment tests whether equalizing N1002K
protein abundance restores function. The resulting drug priorities are:

| Research rank | Candidate | Rationale and proposed test |
|---|---|---|
| 1 | **Arimoclomol** | Test proteostasis-mediated abundance rescue; compare vehicle/drug first-mitosis errors and first-anaphase yield without a scored visible error per starting cell in a qualified variant model, with genetic correction controls. |
| 2 | **Acetylcysteine (NAC)** | Test downstream stress protection if a measured oxidative defect exists; distinguish improved cell function/survival from restored segregation. |

The [Track 2 proposal](reports/jeweled_track2_report_submission.md) connects each drug to its
mechanism, evidence gaps, exposure gate and tumor counter-screen. No drug experiment has run.

## Scope and limitations

**Low SV recall (8.6% DELLY; 21.3% SvABA) limits the causal search:** a negative short-read result
cannot exclude difficult structural alternatives. Annotation locates misses without proving their
cause. One public chromosome and synthetic tests do not establish case-wide clinical sensitivity.
Phase, N1002K function and drug benefit remain unresolved. EPCR is a research judgment.
Exploratory control development closed on 20 September; remaining work concerns reproduction,
documented defects and submission delivery before **24 October 2026**.

## Data handling, disclosure and license

The repository remains private; curated findings are included, raw reads and genome-scale case
payloads are excluded. [Release review](workflow/docs/PUBLIC-RELEASE-AUDIT.md) records decisions.
Claude Code/Codex used personal subscriptions with training disabled (user-confirmed). Narration
uses Cartesia Sonic 3.6/Lauren; Cartesia account settings remain unverified. See [LICENSE](LICENSE),
[required attribution and dataset citation status](reports/PROVENANCE.md#acknowledgement), and [CHANGELOG](CHANGELOG.md).
