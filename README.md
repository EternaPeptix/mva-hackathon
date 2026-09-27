# MVA: from genome evidence to a testable rescue experiment

Team **jeweled** focuses on the case's [recorded muscle atrophy and failure to thrive](workflow/configs/lirical-subject.json):
could a gene-informed medication hypothesis preserve muscle function and the cells that renew it?

**[Review the current Track 2 package](reports/REVIEW.md):** report, route comparison,
revised pitch, editable slides and revision-specific reviewer assessments.

- **[Track 1: candidate and evidence](reports/jeweled_track1_report_submission.md).**
  BUB1B **p.(Leu737Ter) + p.(Asn1002Lys)**, at subjective **EPCR 0.8**.
  Phase is unresolved and the missense remains VUS.
- **[Track 2: selected medications and proposal](reports/jeweled_track2_report_submission.md).**
  **Niacin (nicotinic acid) and apitegromab are the provisional primary research comparison**, tested separately.
  They test usable BUBR1 restoration versus downstream muscle compensation; everolimus remains a conditional mechanism probe. They share the
  same goal: greater absolute muscle output, preserved renewal and no worsened chromosome fidelity.
- **[Submission index, methods and proposed pitch](reports/README.md).**
  **Track 1 model 1 scored 100 rank points and F-max 1.000 on 21 September.**
  [Submission record](reports/README.md#track-1-submission-record) distinguishes the
  confirmed score from the pending final report-persistence check. Track 2 has not been submitted. The new niacin/apitegromab Skylar video is **verified as V3 at the [existing Vimeo URL](https://vimeo.com/1228990028)**, with completed optimization and enabled English captions. No drug experiment has run.

## Start here

**Track 2 update:** the [broad drug search](workflow/track2-drug-census/decision-report.md) is incorporated into the report, methods, experimental plan and revised narration. Niacin tests restoration; apitegromab tests compensation. The census screens 15,357 entities plus 40 separate FDA product entries, with 247 curated route/compound appraisals; retrieval is not evidence of every drug’s efficacy. The new eight-slide movie measures 180 seconds: 166.64 seconds of speech and 13.36 seconds of transitions. Technical frame, audio and caption checks passed; no full independent listening review is claimed. Vimeo V3 replacement, new title/description and signed-out apitegromab-slide/caption correspondence at 01:19 were verified on 21 September 2026. Track 2 submission remains pending.

The [food-name review](workflow/track2-drug-census/food-substances/common-abstract-findings.md) retains 212 available abstracts across 213 source IDs. [Batch15](workflow/track2-drug-census/batch15-revalidation-findings.md) adds the final 532 saved strings: **3,598 accepted revalidations, 5,759 combined outer source/ID records and zero pending in the original selection**. This closes that queue, not the broader search or primary-source follow-up. The earlier 2,161 annotations were not independently re-audited here. Root checked 16 source/annotation pairs in the final batch; versions are not independent studies and the 143 initially empty fields retain separate recovery statuses.

The ledger remains **247 appraisals**, with only NMN further refined and the other 246 decision objects unchanged. A [separate earlier NMN trial](workflow/track2-drug-census/nmn-igarashi2022-qualification.md) adds exploratory motor findings alongside the [previously qualified human trials](workflow/track2-drug-census/nmn-human-functional-qualification.md). Wrong-supplement dispensing excluded 22 of 42 week-12 datasets; the registered primary muscle index was nonpositive. The qualified positives and counterevidence do not establish a dose optimum or MVA efficacy. The preceding [D3 assessment](workflow/track2-drug-census/vitamin-d3-functional-qualification.md) retains a reported within-arm chair-rise positive with unverified controlled contrast, a qualified daily-trial primary leg-power null, and monthly active-comparator function nulls/falls signals. Calcifediol combination results remain separate; no generalized functional benefit or comparative MVA efficacy is established.

The separate [NIH DSLD ingredient inventory](workflow/track2-drug-census/dsld-ingredients/README.md) preserves **6,467 source groups, 78,788 synonym cells and 85,255 total name/synonym cells**. Its **41,496 lexical matches, including 1,077 unflagged matches, provide no efficacy assessment** and confer no reading or appraisal credit. A separate first pass annotates six rows across four previously read abstracts, without new unique credit. Groups, label terms and matches are not unique drugs, verified chemical equivalents or additions to the main entity denominator. API version 9.5.0 (May 2026) and index-name suffix `2025-09-25n` are distinct source descriptors, not current-market or record-update guarantees.

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

6. Qualify allele products, useful BUBR1 function and a myostatin precursor/activation system in controlled muscle models.
7. Test niacin and apitegromab separately at justified exposure, with drug-appropriate, genetic and tumor controls.
8. Measure absolute muscle output, renewal and fidelity; retain everolimus as a conditional mTORC1 probe and qualify reserve routes.

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
Confirmed cis would retire this pair as the compound-heterozygous explanation.
The stop has a two-star ClinVar P/LP assertion and the missense remains VUS. The
[report](reports/jeweled_track1_report_submission.md) explains the two-allele census, ten-gene
differential, 288 bp phase blocks versus a 10,911 bp gap, and retracted mosaic signal.

**Track 2:** nominate **niacin (nicotinic acid)** for BUBR1 restoration and **apitegromab**
for downstream muscle compensation, tested separately. Human muscle NAD engagement
and separate SIRT2/BUBR1 experiments support testing niacin; they do not demonstrate
its restoration chain. Apitegromab has randomized motor-function evidence in SMA,
with age, dose and background-treatment limits, labeled fracture risk and intravenous
access burden. Neither primary candidate has demonstrated MVA rescue. The myostatin
assay must support actual-drug engagement; excess baseline myostatin is not mandatory
for compensation. Everolimus remains a conditional mTORC1 mechanism probe because
harmful pathway excess is unestablished here. Our [RNA reanalysis](workflow/track2-rna/README.md)
adds no mTORC1 support and does not settle the qualitative carrier-blot observation.
Arimoclomol and dasatinib remain conditional reserves.

The [proposal](reports/jeweled_track2_report_submission.md) and
[evidence pipeline](reports/appendices/track2-evidence-pipeline.md) explain the
qualitative selection, exposure gate, absolute-output endpoints and rhabdomyosarcoma counter-screen.
The [unelected-route table](reports/appendices/track2-route-decisions.md) explains alternatives and weaknesses.
The [Cell selection and tissue output explorer](workflow/track2-model/explorer/README.md)
makes the synthetic +16.8-point chromosome-normal fraction / −11.0% raw-output example
interactive offline. It displays retained aggregate populations, not individual lineages.
The [genome-informed model](workflow/track2-model/README.md) separates sequence reconstruction
from assumption-based population scenarios; it does not predict drug response. The
[public-source transfer replay](workflow/track2-transfer/README.md) executes the decision process
on BUB1B and a synthetic TRIP13 example, retaining gene-specific evidence gaps.
It is a second worked example, not validation in two patients. No drug experiment
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
uses Cartesia Sonic 3.6/Skylar on the Free plan; its training/retention settings remain unverified. See [LICENSE](LICENSE),
[required attribution and dataset citation status](reports/PROVENANCE.md#acknowledgement), and [CHANGELOG](CHANGELOG.md).
