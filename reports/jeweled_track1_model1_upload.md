> **Source revision.** Source: [cda0b0c12249b446adeff61cd1838f148ccb1768](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/reports/jeweled_track1_report_submission.md). Repository access is required.

# jeweled — Track 1: a recessive BUB1B hypothesis

**Research submission, revised 21 September 2026.**

Our primary hypothesis is the **BUB1B p.(Leu737Ter) plus p.(Asn1002Lys) pair**, at GRCh38 chr15:40209701 T>G and chr15:40220612 T>G. We assign the exact pair a **subjective causal probability (EPCR) of 0.8**. Trans is favored under the recessive disease model but has not been observed. Leu737Ter has a two-star ClinVar Pathogenic/Likely pathogenic assertion; Asn1002Lys remains a VUS. The [prediction CSV](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/reports/jeweled_evidence-reconciled-v2.csv) contains one primary row and no secondary finding.

## Why this is the leading explanation

The organizers describe a clinically solved MVA case. That supplied disease context is relevant prior evidence, although the public description does not confirm these exact two alleles or their segregation. Our review was informed by the previously nominated BUB1B hypothesis and is not a blinded discovery experiment.

Both GATK and bcftools recover the two heterozygous alleles. Filtered reference/alternate fragment counts are 21/26 and 15/12, with alternate observations on both strands. BUB1B coding bases have minimum mapped depth 24 in both the MQ20/BQ20 and MQ60/BQ20 reviews. A locus census finds only these two exon/CDS candidate overlaps and no additional screened splice-edge overlap. These results make an obvious missed third coding allele less likely.

The allele-specific autosomal inventory finds **two rare protein-altering BUB1B variants and no rare protein-altering or canonical-splice variant in ten comparison genes**: TRIP13, CEP57, CENATAC, MAD2L1BP, SLF2, SMC5, TRIM37, BUB1, BUB3 and MAD2L1. “Rare” here means retained population AF below 1%, a discovery filter rather than an ACMG threshold. Three comparison genes have common missense calls; none of their retained coding alleles has unresolved frequency. Both callers give the same result. Coding depth has no zero-depth base in these ten genes, although small low-depth intervals remain. This supports the differential within its declared variant classes; it is not exclusion of every possible mechanism.

The original vendor annotations have a maximum phase-block span of **288 bp**, against **10,911 bp between the candidates**. An earlier two-lane, 100 kb read-phasing analysis reaches 2,327 bp but does not connect them. The final four-lane bounded review also leaves the pair unphased, with a 6,769 bp intervening marker gap. Failure to bridge this distance is an information limit, not evidence that the variants are in cis or a reason to cap causal probability at one half. Parental segregation is the clearest next test. The [case-evidence appendix](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/reports/appendices/track1-case-evidence.md) connects these claims to selected numerical results and the authorized reconstruction workflow.

## Allele interpretation

The stop allele's external clinical assertion and our reconstruction are reported separately. [ClinVar VCV000533901.9](https://www.ncbi.nlm.nih.gov/clinvar/variation/533901/) combines one Labcorp Pathogenic assertion for MVA1 and one GeneDx Likely pathogenic assertion without conflict. It is not an expert-panel assertion. Relevant-transcript loss of function and nonsense-mediated decay are predicted; patient RNA has not been measured.

| Allele | Applied research criteria | Working category |
|---|---|---|
| Leu737Ter | PVS1 Very Strong + PM2 Supporting: 9 points | Likely pathogenic; external ClinVar **P/LP, two stars** |
| Asn1002Lys | PM2 Supporting + PM3 Supporting: 2 points | VUS |

We accept sufficient rarity of both alleles as an explicit general recessive-disease judgment: the stop's highest retained ancestry AF is approximately 9.98e-5 and the missense's is 8.99e-7. This is not a validated BUB1B-specific cutoff, and rarity alone is weak evidence. PM3 for the missense uses the independent Labcorp **Pathogenic** partner assertion: the unknown-phase observation contributes 0.5 case points, reaching Supporting. We do not use the missense to upgrade the stop, then reuse that upgrade. [ClinGen PM2 guidance](https://www.clinicalgenome.org/site/assets/files/5182/pm2_-_svi_recommendation_-_approved_sept2020.pdf); [PM3 guidance](https://www.clinicalgenome.org/site/assets/files/3717/svi_proposal_for_pm3_criterion_-_version_1.pdf).

PP4 is retained as a disclosed sensitivity because the supplied MVA phenotype is informative but a sufficiently supported, phenotype-matched BUB1B diagnostic-yield estimate is not established here. Genetic heterogeneity alone does not prohibit PP4. Adding Supporting PP4 would yield stop 10 points/Pathogenic and missense 3 points/VUS. The [interpretation appendix](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/reports/appendices/causal-interpretation.md) gives the framework, rarity judgments and remaining evidence gaps. Allele categories are not a numerical ceiling on the case-level hypothesis.

## Differential and causal probability

The broader review found no stronger explanation. Joint reconstruction changes the adjacent ADAMTS1 deletions from isolated frameshift labels to a net in-frame sequence; FAM136A's combined sequence is missense. Other leads retain mapping, phase or effect uncertainty, including an AFDN/KIF25/FRMD1 structural event without an established relevant pathogenic consequence. The apparent chromosome 20/21/22 mosaic-aneuploidy signal was **retracted** after centromeric/acrocentric masking removed the apparent chromosome-wide excess. This removes our claimed computational support; it does not refute the organizer-reported diagnosis.

The 0.8 estimate weighs the supplied MVA context, known BUB1B mechanism, credible stop allele, extremely rare second allele, locus census and completed differential. The exact-pair argument is narrower than gene-level support: N1002K is the only observed rare protein-altering partner in BUB1B, with credible read support and no additional emitted coding or screened splice-edge candidate. We favor this observed pair over a hypothesis requiring a separate undetected second hit, while recognizing that the reviewed variant classes do not cover every such hit and rarity does not establish missense dysfunction. Trans compatibility and a damaging N1002K effect remain assumptions of the proposed explanation.

The numerical EPCR expresses our weighting of this evidence. It is not uniquely determined by the data or proven better calibrated than a lower estimate. No new segregation or functional observation caused the revision; correcting an earlier overstatement can legitimately change a belief without a new adverse biological finding.

The residual uncertainty has identifiable explanations: Asn1002Lys could be incidental or insufficiently damaging, with an undetected BUB1B regulatory/structural second hit; the observed pair could be in cis; or another incompletely detected or unrecognized mechanism could explain the case. We do not invent numerical shares for these alternatives. EPCR is a stated research judgment, not a fitted posterior, prioritizer score, variant-classification probability or probability of treatment benefit. No empirical probability interval is claimed.

## Method, validation and limits

We verified the eight supplied FASTQ checksums, aligned four sequencing lanes to GRCh38 on Mac Studios, and ran GATK 4.6.1.0 and bcftools 1.24. Original genotypes were preserved through normalization and annotation. Autosomal coding inventories cover chromosomes 1–22; sex-chromosome, mitochondrial, copy-number and structural analyses have separate scopes. Caller agreement uses the same reads.

The completed rankings used the **supplied VCF**, not the later four-lane native calls. Exomiser 15.1.0 used hg38/phenotype resources 2512; LIRICAL 2.4.1 used HPO 2026-09-01 and shared Exomiser allele databases. These comparisons are not independent biological replication.

| Ranking comparison | BUB1B/MVA1 result | Interpretation |
|---|---|---|
| Exomiser, seven subject features, primary predictors | BUB1B AR **rank 2**, score **0.4877**; AD rank 1 | Rank unit is gene × inheritance; the AD result does not establish dominant MVA. |
| Exomiser, REVEL only, strict / permissive filter | AR **absent / rank 3**, score **0.1806** when retained | The exact pair is sensitive to predictor/filter choice. |
| LIRICAL, phenotype only, full / without rhabdomyosarcoma / rhabdomyosarcoma alone | MVA1 **rank 4 / 38 / 10** among 8,525 diseases | The cancer phenotype contributes materially to ranking. |
| LIRICAL, full global genotype-assisted JSON | MVA1 **rank 1/8,525**, model probability **0.9575393** | A model result, not EPCR, measured missense function or phase. |

N1002K has **REVEL 0.472** and **AlphaMissense 0.9229**. The adopted REVEL policy treats 0.472 as indeterminate (neither PP3 nor BP4); the alternative score is disclosed without stacking predictors or upgrading the classification. Gene prominence is more robust than the exact partner's computational support. [Complete sensitivities, resource identities and reconstruction map](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/reports/appendices/track1-case-evidence.md); [selected aggregate results](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/results/track1-selected-case-evidence.json).

Public HG002 chromosome-20 validation uses reads labelled nominally 35×, GIAB v4.2.1 confident regions and RTG vcfeval. Measured mean read depth is 35.75× across chr20 and 36.97× within the confident regions at MAPQ≥20/BQ≥20; 99.25% of confident-region bases have at least 20× depth. These descriptive filters do not reconstruct either caller's internal filtering. Original GATK SNV/indel recall was 99.3748%/97.3902%, with precision 98.9548%/96.0168%; preserving native multiallelic genotypes improves the later indel comparison. These benchmark results are not an estimate of case sensitivity. Public SV matching was low, so negative short-read evidence cannot exclude difficult insertions or repeat-associated events. The [methods form](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/reports/jeweled_methods_form_answers.md) separates original results from later representation analyses.

No parental, patient RNA/protein, functional or tumor measurements were available. The [secondary inventory](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/reports/appendices/secondary-findings.md) nominates no incidental finding; its [coverage review](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/reports/appendices/secondary-depth.md) retains local gaps. Read depth, negative inventories and the boundary-zero [nuclear-mixture estimate](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/reports/appendices/nuclear-mixture.md) do not establish complete clinical callability, purity or absence of mosaicism. The most consequential next evidence is segregation and a qualified N1002K abundance/function assay.

**AI use and data handling:** Anthropic Claude Code and OpenAI Codex, personal subscriptions, training on the user's content disabled (user-confirmed). Exact historical model/plan names are unverified. Computation used owned machines; candidate-level evidence summaries were reviewed with hosted assistants. The repository remains private. Curated findings are separated from restricted raw reads and genome-scale case payloads. This is a research hypothesis, not a clinical sign-out.

## Methods form — Track 1, model 1

### Team name

jeweled

### Model number (fill out this form once per model/approach submitted; up to 6)

1

### Please describe your model/approach in detail.

We ask whether the supplied MVA case is explained by a recessive BUB1B pair. Eight FASTQ checksums were verified, four sequencing lanes were aligned to GRCh38 on Mac Studios, and GATK 4.6.1.0 plus bcftools 1.24 were run independently on the same reads. Their completed whole-reference outputs were reconciled to the original call sections. Coding inventories cover all 22 autosomes; fixed-ploidy sex-chromosome, mitochondrial, copy-number and structural analyses have separate scopes.

We preserved original multiallelic genotypes, missing fields and annotation lineages. LIRICAL/Exomiser rankings used the supplied VCF; no ranking rerun on the final four-lane native callsets is claimed. Exomiser's BUB1B recessive entry ranked 2 (score 0.4877) with primary predictors, disappeared under REVEL-only strict filtering, and ranked 3 (0.1806) under REVEL-only permissive filtering. LIRICAL's MVA1 rank was 4/8,525 with seven phenotypes, 38 without rhabdomyosarcoma, 10 with rhabdomyosarcoma alone, and 1 with genotype integration (model probability 0.9575393, not EPCR). The exact pair is less robust to predictor choice than the gene's prominence. A separate final native allele-specific coding inventory used frequency filters without pathogenicity predictors. Regional read and joint-sequence review then tested alternative explanations. The [case-evidence appendix](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/reports/appendices/track1-case-evidence.md) supplies the complete comparison, selected aggregates and authorized execution map.

Both callers find the same two BUB1B protein-altering alleles and no rare protein-altering/canonical-splice allele in ten comparison genes. The current prediction is that BUB1B pair at subjective EPCR **0.8**, with phase unknown and Asn1002Lys still VUS. The [Track 1 report](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/reports/jeweled_track1_report_submission.md) and [interpretation appendix](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/reports/appendices/causal-interpretation.md) explain the inference. The previously nominated candidate and supplied solved-case context were known; this was not blinded discovery.

### Please state whether the submission file is the automated output of your computational approach (preferred), or if it has undergone downstream manual review and curation.

The CSV underwent manual review with AI assistance. The existing builder validates an explicit candidate decision, GRCh38 reference alleles, phase field and probability rationale. Selection and EPCR are research judgments, not calibrated automated outputs. The new one-row export replaces both the historical 0.95 and the later 0.35 judgments; earlier files remain archival.

### Please describe any downstream manual review in detail.

We separated allele presence, phase, allele classification and case-level causality. Leu737Ter has an external ClinVar Pathogenic/Likely pathogenic two-star assertion. Our declared points-based reconstruction applies PVS1 Very Strong plus PM2 Supporting (9 points, LP); Asn1002Lys receives PM2 Supporting plus PM3 Supporting (2 points, VUS). Both rarity judgments are explicitly accepted under general recessive-disease guidance. PM3 uses the independent Labcorp Pathogenic partner assertion and the unknown-phase allowance; no circular upgrade is made. PP4 is a disclosed sensitivity, not rejected merely because MVA is genetically heterogeneous. [Criteria and primary guidance](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/reports/appendices/causal-interpretation.md).

The vendor phase-block maximum was 288 bp versus a 10,911 bp intervariant gap. The earlier two-lane 100 kb review reached 2,327 bp; the final four-lane bounded review also left the pair unknown, with a 6,769 bp marker gap. This is missing information rather than evidence of cis. We favored trans under the recessive model without recording observed trans phase. Joint sequence review corrected misleading ADAMTS1 and FAM136A isolated annotations. An apparent chromosome-wide mosaic signal was retracted after masking localized bias. Remaining alternatives include a hidden second BUB1B hit, an incidental or insufficiently damaging missense allele, and incompletely detected mechanisms elsewhere.

EPCR 0.8 reflects the supplied disease context, gene mechanism, supported rare alleles and differential. It is not converted from ACMG points or prioritizer scores. No empirical confidence interval is claimed. No new functional or segregation observation caused this interpretation revision.

### Please state whether your approach only used publicly available data (preferred) or if proprietary data was also used. If proprietary data were used, where possible, please also include a submission based only on publicly available data.

We used the challenge-supplied case data plus public references and literature. No additional proprietary cohort, training dataset or laboratory dataset was used. Challenge case data are controlled access, not described as public.

### Please describe any public data used in detail (e.g., reference panels, population databases, literature).

| Resource | Use |
|---|---|
| GRCh38 with retained decoy/reference context | Alignment, calling, transcript reconstruction and reference-allele checks. |
| Ensembl release 116 GFF/transcripts | Consequence annotation and all-transcript CDS targets. |
| Exomiser 15.1.0 / release-2512 resources | Phenotype, inheritance, predictor and frequency ranking. |
| LIRICAL 2.4.1 / HPO 2026-09-01 / Ensembl transcript database | Phenotype-only and genotype-integrated comparisons; exact public-resource identities in the [selected evidence](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/results/track1-selected-case-evidence.json). |
| gnomAD v4.1, ClinGen and targeted ClinVar records | Exact-allele frequencies, gene validity and classification guidance. Missing observations are retained. |
| ClinVar GRCh38 snapshot dated 2026-09-13 | Frozen exact-allele and assertion-lineage review. |
| ACMG SF v3.3 / frozen ClinGen panel | Secondary inventory of 84 genes and reference-only coding targets. |
| Public nuclear marker resource | A fixed 37,267-marker model check, without a validated biological detection limit. |
| GIAB HG002 chr20 / v4.2.1 small-variant truth and confident BED | Public SNV/indel benchmarking. |
| GIAB HG002 structural truth / v5.0q, chr20 subset | Public structural benchmarking; [truth composition and residuals](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/reports/appendices/public-sv-residuals.md). |

HG002 reads are labelled nominally 35×. A subsequent descriptive read-depth measurement at MAPQ≥20/BQ≥20 found mean depth 35.75× across chr20 and 36.97× in the official confident regions; 99.25% of confident-region bases had at least 20× depth. Zero-depth positions were included; overlapping mates both count, deletion bases do not, and default unmapped/secondary/QC-fail/duplicate exclusions apply. This is not either caller's internal read-filter definition. RTG Tools 3.13 vcfeval compared diploid, raw calls in those confident regions, with full-chromosome context. Original GATK SNV/indel recall was 99.3748%/97.3902%, precision 98.9548%/96.0168%; bcftools recall was 99.3593%/91.5115%, precision 98.9512%/91.1827%. Later retained-multiallelic results are a separate representation sensitivity. We do not claim an additional hap.py cross-check. DELLY and SvABA matched 69/802 and 171/802 structural truth alleles; insertions and repeat-associated misses limit negative case interpretation.

N1002K scores are REVEL **0.472** and AlphaMissense **0.9229**. The adopted REVEL policy gives neither PP3 nor BP4; AlphaMissense remains a disclosed alternative, not stacked evidence. Historical provenance identifies REVEL v1.3 and the primary AlphaMissense hg38 release, without an exact AlphaMissense release date. LIRICAL's maximum-predictor rule and shared Exomiser resources do not provide independent confirmation. Scores were not multiplied into EPCR, and no new predictive model was trained. [Sensitivity details and resource-date limitations](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/reports/appendices/track1-case-evidence.md#ranking-and-predictor-sensitivity).

### Please describe any proprietary data used in detail.

The nonpublic inputs were organizer-supplied phenotype, variants and sequencing reads. No parental or tumor specimen data, patient-derived experimental observations or drug-response measurements were available.

### Is your approach able to output proposed pairs of compound heterozygous candidate variants, or only single candidate variants?

Pairs. The primary row contains heterozygous BUB1B p.(Leu737Ter) and p.(Asn1002Lys). Trans is the favored inference under the proposed recessive explanation, but phase remains explicitly unknown in the export. Phenotype or model ranking is not an observed segregation result.

### How did you handle secondary or incidental findings, if any? (See the finding_type column in your submission.)

The frozen 84-gene secondary review retained exact coding and noncoding matches, original genotypes, missing quality fields and assertion lineage. One called SCN5A allele per caller was individually Benign/Likely benign and selected only through a containing-haplotype assertion. It is not nominated. No exact called allele entered a clean aggregate P/LP stratum in the frozen resources. [Inventory](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/reports/appendices/secondary-findings.md).

Across the 393,880-base coding union, 99.30% had depth at least 20 at MQ20/BQ20 and 98.59% at MQ60/BQ20; the latter had 2,834 zero-depth bases in five genes. These are descriptive coverage results, not a clinical negative. [Coverage and local gaps](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/reports/appendices/secondary-depth.md).

The separate nuclear-mixture estimator returned a boundary value of 0.0 after 37,149/37,267 markers passed coverage gates. This does not establish purity, mosaicism sensitivity or variant phase. [Model limitations](https://github.com/EternaPeptix/mva-hackathon/blob/cda0b0c12249b446adeff61cd1838f148ccb1768/reports/appendices/nuclear-mixture.md).

### Please provide an estimate of run time and cost for your approach (if possible).

No validated end-to-end runtime or consolidated compute, electricity, subscription or labor cost was measured for the case analysis. A separate fresh public HG002 chr20 replay completed in 134.3 minutes, including downloads, setup, both callers, comparison and audit, with 5.31 GiB peak sampled process-group memory. It exactly reproduced the historical primary counts and audit partitions. The older 116-minute measurement covered calling only. These public timings are not the case's full analysis time or guarantees for other machines. Mac Studios performed the current genome work. Prospective resource ceilings and assay-capacity estimates are planning values.

### Provide a method abstract (up to 500 words). Please include strengths and limitations of your method as well as methodological detail.

We investigated a supplied MVA case using verified sequencing reads, GRCh38 alignment, two variant callers and allele-specific coding inventories. Phenotype/genotype ranking was supplemented by predictor sensitivities, population review, read phasing and joint-sequence reconstruction. Copy-number, structural, mitochondrial and secondary-finding analyses retain distinct scopes.

The leading hypothesis is BUB1B p.(Leu737Ter) plus p.(Asn1002Lys), at subjective EPCR 0.8. Both callers recover two rare BUB1B coding alleles and no rare protein-altering/canonical-splice alternative in ten comparison genes. Phase is unobserved: the 10,911 bp gap exceeds the original 288 bp phase-block maximum and later local blocks. Leu737Ter has a two-star external P/LP assertion; the declared research reconstruction is LP, while Asn1002Lys remains VUS. The supplied solved-case context informs the inference, so this is not blinded discovery.

Strengths include preservation of original genotypes, investigation of alternative sequence representations and correction of a misleading mosaic signal. Public HG002 chr20 benchmarking uses GIAB confident regions and RTG vcfeval, with separate SNV/indel results. Synthetic tests check implementation behavior. Neither establishes case-wide clinical sensitivity.

Limitations include one sample, absent parental and functional evidence, incomplete sensitivity to regulatory/structural/repeat mechanisms and no measured treatment response. Residual hypotheses are stated rather than hidden behind missing phase. The next decisive studies are segregation and a qualified N1002K abundance/function assay.

## Acknowledgement and data provenance

> This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.

Controlled case inputs were supplied through the [MVA Hackathon dataset](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data) and [Synapse project syn76251147](https://www.synapse.org/Synapse:syn76251147/wiki/642892). This report is submitted under the event CC BY 4.0 terms; underlying controlled case data are not redistributed. Repository license and third-party exclusions remain as documented.
