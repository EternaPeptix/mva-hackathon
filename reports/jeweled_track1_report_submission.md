# jeweled — Track 1: a recessive BUB1B hypothesis

**Research submission draft, revised 21 September 2026; not submitted.**

Our primary hypothesis is the **BUB1B p.(Leu737Ter) plus p.(Asn1002Lys) pair**, at GRCh38 chr15:40209701 T>G and chr15:40220612 T>G. We assign the exact pair a **subjective causal probability (EPCR) of 0.8**. Trans is favored under the recessive disease model but has not been observed. Leu737Ter has a two-star ClinVar Pathogenic/Likely pathogenic assertion; Asn1002Lys remains a VUS. The [prediction CSV](jeweled_evidence-reconciled-v2.csv) contains one primary row and no secondary finding.

## Why this is the leading explanation

The organizers describe a clinically solved MVA case. That supplied disease context is relevant prior evidence, although the public description does not confirm these exact two alleles or their segregation. Our review was informed by the previously nominated BUB1B hypothesis and is not a blinded discovery experiment.

Both GATK and bcftools recover the two heterozygous alleles. Filtered reference/alternate fragment counts are 21/26 and 15/12, with alternate observations on both strands. BUB1B coding bases have minimum mapped depth 24 in both the MQ20/BQ20 and MQ60/BQ20 reviews. A locus census finds only these two exon/CDS candidate overlaps and no additional screened splice-edge overlap. These results make an obvious missed third coding allele less likely.

The allele-specific autosomal inventory finds **two rare protein-altering BUB1B variants and no rare protein-altering or canonical-splice variant in ten comparison genes**: TRIP13, CEP57, CENATAC, MAD2L1BP, SLF2, SMC5, TRIM37, BUB1, BUB3 and MAD2L1. “Rare” here means retained population AF below 1%, a discovery filter rather than an ACMG threshold. Three comparison genes have common missense calls; none of their retained coding alleles has unresolved frequency. Both callers give the same result. Coding depth has no zero-depth base in these ten genes, although small low-depth intervals remain. This supports the differential within its declared variant classes; it is not exclusion of every possible mechanism.

The original vendor annotations have a maximum phase-block span of **288 bp**, against **10,911 bp between the candidates**. A later local read-phasing analysis reaches 2,327 bp but still does not connect them. Failure to bridge this distance is an information limit, not evidence that the variants are in cis or a reason to cap causal probability at one half. Parental segregation is the clearest next test.

## Allele interpretation

The stop allele's external clinical assertion and our reconstruction are reported separately. [ClinVar VCV000533901.9](https://www.ncbi.nlm.nih.gov/clinvar/variation/533901/) combines one Labcorp Pathogenic assertion for MVA1 and one GeneDx Likely pathogenic assertion without conflict. It is not an expert-panel assertion. Relevant-transcript loss of function and nonsense-mediated decay are predicted; patient RNA has not been measured.

| Allele | Applied research criteria | Working category |
|---|---|---|
| Leu737Ter | PVS1 Very Strong + PM2 Supporting: 9 points | Likely pathogenic; external ClinVar **P/LP, two stars** |
| Asn1002Lys | PM2 Supporting + PM3 Supporting: 2 points | VUS |

We accept sufficient rarity of both alleles as an explicit general recessive-disease judgment: the stop's highest retained ancestry AF is approximately 9.98e-5 and the missense's is 8.99e-7. This is not a validated BUB1B-specific cutoff, and rarity alone is weak evidence. PM3 for the missense uses the independent Labcorp **Pathogenic** partner assertion: the unknown-phase observation contributes 0.5 case points, reaching Supporting. We do not use the missense to upgrade the stop, then reuse that upgrade. [ClinGen PM2 guidance](https://www.clinicalgenome.org/site/assets/files/5182/pm2_-_svi_recommendation_-_approved_sept2020.pdf); [PM3 guidance](https://www.clinicalgenome.org/site/assets/files/3717/svi_proposal_for_pm3_criterion_-_version_1.pdf).

PP4 is retained as a disclosed sensitivity because the supplied MVA phenotype is informative but a sufficiently supported, phenotype-matched BUB1B diagnostic-yield estimate is not established here. Genetic heterogeneity alone does not prohibit PP4. Adding Supporting PP4 would yield stop 10 points/Pathogenic and missense 3 points/VUS. The [interpretation appendix](appendices/causal-interpretation.md) gives the framework, rarity judgments and remaining evidence gaps. Allele categories are not a numerical ceiling on the case-level hypothesis.

## Differential and causal probability

The broader review found no stronger explanation. Joint reconstruction changes the adjacent ADAMTS1 deletions from isolated frameshift labels to a net in-frame sequence; FAM136A's combined sequence is missense. Other leads retain mapping, phase or effect uncertainty, including an AFDN/KIF25/FRMD1 structural event without an established relevant pathogenic consequence. The apparent chromosome 20/21/22 mosaic-aneuploidy signal was **retracted** after centromeric/acrocentric masking removed the apparent chromosome-wide excess. This removes our claimed computational support; it does not refute the organizer-reported diagnosis.

The 0.8 estimate weighs the supplied MVA context, known BUB1B mechanism, credible stop allele, extremely rare second allele, locus census and completed differential. The exact-pair argument is narrower than gene-level support: N1002K is the only observed rare protein-altering partner in BUB1B, with credible read support and no additional emitted coding or screened splice-edge candidate. We favor this observed pair over a hypothesis requiring a separate undetected second hit, while recognizing that the reviewed variant classes do not cover every such hit and rarity does not establish missense dysfunction. Trans compatibility and a damaging N1002K effect remain assumptions of the proposed explanation.

The numerical EPCR expresses our weighting of this evidence. It is not uniquely determined by the data or proven better calibrated than a lower estimate. No new segregation or functional observation caused the revision; correcting an earlier overstatement can legitimately change a belief without a new adverse biological finding.

The residual uncertainty has identifiable explanations: Asn1002Lys could be incidental or insufficiently damaging, with an undetected BUB1B regulatory/structural second hit; the observed pair could be in cis; or another incompletely detected or unrecognized mechanism could explain the case. We do not invent numerical shares for these alternatives. EPCR is a stated research judgment, not a fitted posterior, prioritizer score, variant-classification probability or probability of treatment benefit. No empirical probability interval is claimed.

## Method, validation and limits

We verified the eight supplied FASTQ checksums, aligned four sequencing lanes to GRCh38 on Mac Studios, and ran GATK 4.6.1.0 and bcftools 1.24. Original genotypes were preserved through normalization and annotation. Autosomal coding inventories cover chromosomes 1–22; sex-chromosome, mitochondrial, copy-number and structural analyses have separate scopes. LIRICAL and Exomiser supplied ranking and sensitivity comparisons, not independent replication. Caller agreement uses the same reads.

Public HG002 chromosome-20 validation uses reads labelled nominally 35×, GIAB v4.2.1 confident regions and RTG vcfeval. Measured mean read depth is 35.75× across chr20 and 36.97× within the confident regions at MAPQ≥20/BQ≥20; 99.25% of confident-region bases have at least 20× depth. These descriptive filters do not reconstruct either caller's internal filtering. Original GATK SNV/indel recall was 99.3748%/97.3902%, with precision 98.9548%/96.0168%; preserving native multiallelic genotypes improves the later indel comparison. These benchmark results are not an estimate of case sensitivity. Public SV matching was low, so negative short-read evidence cannot exclude difficult insertions or repeat-associated events. The [methods form](jeweled_methods_form_answers.md) separates original results from later representation analyses.

No parental, patient RNA/protein, functional or tumor measurements were available. The [secondary inventory](appendices/secondary-findings.md) nominates no incidental finding; its [coverage review](appendices/secondary-depth.md) retains local gaps. Read depth, negative inventories and the boundary-zero [nuclear-mixture estimate](appendices/nuclear-mixture.md) do not establish complete clinical callability, purity or absence of mosaicism. The most consequential next evidence is segregation and a qualified N1002K abundance/function assay.

**AI use and data handling:** Anthropic Claude Code and OpenAI Codex, personal subscriptions, training on the user's content disabled (user-confirmed). Exact historical model/plan names are unverified. Computation used owned machines; candidate-level evidence summaries were reviewed with hosted assistants. The repository remains private. Curated findings are separated from restricted raw reads and genome-scale case payloads. No competition submission or clinical sign-out has occurred.
