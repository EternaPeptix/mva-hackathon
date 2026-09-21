# Methods description form — answers

**Team: jeweled. Research submission draft, revised 21 September 2026; not submitted.** The 22 question headings reproduce the original methods workbook. The [current submission index](README.md) identifies the active CSV, reports and pitch; dated earlier answers are archival.

The private repository is [EternaPeptix/mva-hackathon](https://github.com/EternaPeptix/mva-hackathon). It contains reproducible source, public/synthetic validation and selected report summaries. Raw case reads, genome-scale call payloads and internal working history remain excluded. Public HG002 replay is separate from authorized case analysis.

**AI use:** Anthropic Claude Code and OpenAI Codex, personal subscriptions, training on the user's content disabled (user-confirmed). Exact historical model/plan names are unverified. Computation used owned machines; candidate-level variants, evidence summaries and report text were reviewed with hosted assistants. This does not imply every historical targeted public-service query ran locally.

**Narration:** Cartesia Sonic 3.6 API, stock Lauren voice, approved by the user. Cartesia plan/tier and training/retention settings are unverified; Claude/Codex settings do not establish them. The service did not return an immutable model snapshot. The approved 180-second movie is attached to a private GitHub release; Vimeo upload remains pending.

# Track 1 (Variant Prediction) — Model 1

### Team name

jeweled

### Model number (fill out this form once per model/approach submitted; up to 6)

1

### Please describe your model/approach in detail.

We ask whether the supplied MVA case is explained by a recessive BUB1B pair. Eight FASTQ checksums were verified, four sequencing lanes were aligned to GRCh38 on Mac Studios, and GATK 4.6.1.0 plus bcftools 1.24 were run independently on the same reads. Their completed whole-reference outputs were reconciled to the original call sections. Coding inventories cover all 22 autosomes; fixed-ploidy sex-chromosome, mitochondrial, copy-number and structural analyses have separate scopes.

We preserved original multiallelic genotypes, missing fields and annotation lineages. LIRICAL/Exomiser rankings were compared with phenotype and predictor sensitivities. A separate allele-specific coding inventory used frequency filters without pathogenicity predictors. Regional read and joint-sequence review then tested alternative explanations.

Both callers find the same two BUB1B protein-altering alleles and no rare protein-altering/canonical-splice allele in ten comparison genes. The current prediction is that BUB1B pair at subjective EPCR **0.8**, with phase unknown and Asn1002Lys still VUS. The [Track 1 report](jeweled_track1_report_submission.md) and [interpretation appendix](appendices/causal-interpretation.md) explain the inference. The previously nominated candidate and supplied solved-case context were known; this was not blinded discovery.

### Please state whether the submission file is the automated output of your computational approach (preferred), or if it has undergone downstream manual review and curation.

The CSV underwent manual review with AI assistance. The existing builder validates an explicit candidate decision, GRCh38 reference alleles, phase field and probability rationale. Selection and EPCR are research judgments, not calibrated automated outputs. The new one-row export replaces both the historical 0.95 and the later 0.35 judgments; earlier files remain archival.

### Please describe any downstream manual review in detail.

We separated allele presence, phase, allele classification and case-level causality. Leu737Ter has an external ClinVar Pathogenic/Likely pathogenic two-star assertion. Our declared points-based reconstruction applies PVS1 Very Strong plus PM2 Supporting (9 points, LP); Asn1002Lys receives PM2 Supporting plus PM3 Supporting (2 points, VUS). Both rarity judgments are explicitly accepted under general recessive-disease guidance. PM3 uses the independent Labcorp Pathogenic partner assertion and the unknown-phase allowance; no circular upgrade is made. PP4 is a disclosed sensitivity, not rejected merely because MVA is genetically heterogeneous. [Criteria and primary guidance](appendices/causal-interpretation.md).

The vendor phase-block maximum was 288 bp versus a 10,911 bp intervariant gap; later local phasing still did not bridge the pair. This is missing information rather than evidence of cis. We favored trans under the recessive model without recording observed trans phase. Joint sequence review corrected misleading ADAMTS1 and FAM136A isolated annotations. An apparent chromosome-wide mosaic signal was retracted after masking localized bias. Remaining alternatives include a hidden second BUB1B hit, an incidental or insufficiently damaging missense allele, and incompletely detected mechanisms elsewhere.

EPCR 0.8 reflects the supplied disease context, gene mechanism, supported rare alleles and differential. It is not converted from ACMG points or prioritizer scores. No empirical confidence interval is claimed. No new functional or segregation observation caused this interpretation revision.

### Please state whether your approach only used publicly available data (preferred) or if proprietary data was also used. If proprietary data were used, where possible, please also include a submission based only on publicly available data.

We used the challenge-supplied case data plus public references and literature. No additional proprietary cohort, training dataset or laboratory dataset was used. Challenge case data are controlled access, not described as public.

### Please describe any public data used in detail (e.g., reference panels, population databases, literature).

| Resource | Use |
|---|---|
| GRCh38 with retained decoy/reference context | Alignment, calling, transcript reconstruction and reference-allele checks. |
| Ensembl release 116 GFF/transcripts | Consequence annotation and all-transcript CDS targets. |
| Exomiser 15.1.0 / release-2512 resources | Phenotype, inheritance, predictor and frequency ranking. |
| LIRICAL 2.4.1 and dated HPO/disease resources | Phenotype-only and genotype-integrated comparisons. |
| gnomAD v4.1, ClinGen and targeted ClinVar records | Exact-allele frequencies, gene validity and classification guidance. Missing observations are retained. |
| ClinVar GRCh38 snapshot dated 2026-09-13 | Frozen exact-allele and assertion-lineage review. |
| ACMG SF v3.3 / frozen ClinGen panel | Secondary inventory of 84 genes and reference-only coding targets. |
| Public nuclear marker resource | A fixed 37,267-marker model check, without a validated biological detection limit. |
| GIAB HG002 chr20 / v4.2.1 truth and confident BED | Public small-variant and structural benchmarking. |

HG002 reads are labelled nominally 35×. A subsequent descriptive read-depth measurement at MAPQ≥20/BQ≥20 found mean depth 35.75× across chr20 and 36.97× in the official confident regions; 99.25% of confident-region bases had at least 20× depth. Zero-depth positions were included; overlapping mates both count, deletion bases do not, and default unmapped/secondary/QC-fail/duplicate exclusions apply. This is not either caller's internal read-filter definition. RTG Tools 3.13 vcfeval compared diploid, raw calls in those confident regions, with full-chromosome context. Original GATK SNV/indel recall was 99.3748%/97.3902%, precision 98.9548%/96.0168%; bcftools recall was 99.3593%/91.5115%, precision 98.9512%/91.1827%. Later retained-multiallelic results are a separate representation sensitivity. We do not claim an additional hap.py cross-check. DELLY and SvABA matched 69/802 and 171/802 structural truth alleles; insertions and repeat-associated misses limit negative case interpretation.

REVEL and AlphaMissense informed the declared ranking/sensitivity policies; their scores were not multiplied into EPCR. No new predictive model was trained.

### Please describe any proprietary data used in detail.

The nonpublic inputs were organizer-supplied phenotype, variants and sequencing reads. No parental or tumor specimen data, patient-derived experimental observations or drug-response measurements were available.

### Is your approach able to output proposed pairs of compound heterozygous candidate variants, or only single candidate variants?

Pairs. The primary row contains heterozygous BUB1B p.(Leu737Ter) and p.(Asn1002Lys). Trans is the favored inference under the proposed recessive explanation, but phase remains explicitly unknown in the export. Phenotype or model ranking is not an observed segregation result.

### How did you handle secondary or incidental findings, if any? (See the finding_type column in your submission.)

The frozen 84-gene secondary review retained exact coding and noncoding matches, original genotypes, missing quality fields and assertion lineage. One called SCN5A allele per caller was individually Benign/Likely benign and selected only through a containing-haplotype assertion. It is not nominated. No exact called allele entered a clean aggregate P/LP stratum in the frozen resources. [Inventory](appendices/secondary-findings.md).

Across the 393,880-base coding union, 99.30% had depth at least 20 at MQ20/BQ20 and 98.59% at MQ60/BQ20; the latter had 2,834 zero-depth bases in five genes. These are descriptive coverage results, not a clinical negative. [Coverage and local gaps](appendices/secondary-depth.md).

The separate nuclear-mixture estimator returned a boundary value of 0.0 after 37,149/37,267 markers passed coverage gates. This does not establish purity, mosaicism sensitivity or variant phase. [Model limitations](appendices/nuclear-mixture.md).

### Please provide an estimate of run time and cost for your approach (if possible).

No validated end-to-end runtime or consolidated compute, electricity, subscription or labor cost was measured. Stage timings have different boundaries and cannot be summed reliably. The earlier public chr20 caller run took approximately 116 minutes; that is a benchmark-stage result, not this case's full analysis time. Mac Studios performed the current genome work. Prospective resource ceilings and assay-capacity estimates are planning values.

### Provide a method abstract (up to 500 words). Please include strengths and limitations of your method as well as methodological detail.

We investigated a supplied MVA case using verified sequencing reads, GRCh38 alignment, two variant callers and allele-specific coding inventories. Phenotype/genotype ranking was supplemented by predictor sensitivities, population review, read phasing and joint-sequence reconstruction. Copy-number, structural, mitochondrial and secondary-finding analyses retain distinct scopes.

The leading hypothesis is BUB1B p.(Leu737Ter) plus p.(Asn1002Lys), at subjective EPCR 0.8. Both callers recover two rare BUB1B coding alleles and no rare protein-altering/canonical-splice alternative in ten comparison genes. Phase is unobserved: the 10,911 bp gap exceeds the original 288 bp phase-block maximum and later local blocks. Leu737Ter has a two-star external P/LP assertion; the declared research reconstruction is LP, while Asn1002Lys remains VUS. The supplied solved-case context informs the inference, so this is not blinded discovery.

Strengths include preservation of original genotypes, investigation of alternative sequence representations and correction of a misleading mosaic signal. Public HG002 chr20 benchmarking uses GIAB confident regions and RTG vcfeval, with separate SNV/indel results. Synthetic tests check implementation behavior. Neither establishes case-wide clinical sensitivity.

Limitations include one sample, absent parental and functional evidence, incomplete sensitivity to regulatory/structural/repeat mechanisms and no measured treatment response. Residual hypotheses are stated rather than hidden behind missing phase. The next decisive studies are segregation and a qualified N1002K abundance/function assay.

# Track 2 (Drug Repositioning)

### Team name

jeweled

### Please describe your approach in detail (how you went from variant/mechanism to candidate medication(s)).

The lead therapeutic thesis is to restore the abundance of residual, functional BUBR1. The stop-plus-missense genotype motivates testing whether N1002K retains useful mitotic function when its abundance is restored. Selected BUBR1 mutants have shown abundance-related rescue in published expression experiments; N1002K has not.

Arimoclomol is the first experimental drug priority if a rescuable protein-quality defect is demonstrated. NAC is a downstream stress-protection comparator. The decisive genetic experiment precedes drug interpretation: if near-physiological N1002K abundance cannot restore function, an abundance-only strategy fails. Primary literature, exposure and tumor context inform this selection. [Track 2 report](jeweled_track2_report_submission.md).

Drug experiments require measured exposure, genetic correction controls and functional endpoints. First-division errors are interpreted with first-anaphase yield without a scored visible error per starting cell, so arrest or selection cannot masquerade as rescue. A tumor-protection counter-screen is a separate advancement gate. No drug-response experiment has been performed.

### Was candidate identification automated (e.g., via a drug-target database query or LLM-assisted search), or primarily manual literature review and expert curation?

Candidate identification and review were AI-assisted, using Claude Code and Codex for source searches, comparison, mechanistic critique and drafting. Selection is a curated research proposal. No validated automated drug-ranking model or independent clinical expert sign-off is claimed.

### Please describe any manual literature review or expert curation in detail.

We checked each bridge against its actual variant, protein, species, model, endpoint and exposure. Published BUBR1 expression rescue provides a rationale for abundance experiments. HSP90 dependence does not establish arimoclomol/HSP70 rescue; the drug-specific bridge remains a question. NAC evidence was separated by model and endpoint, and free reduced parent exposure distinguished from total NAC or cysteine/glutathione changes.

Checkpoint inhibition opposes the intended nonmalignant rescue mechanism and is excluded from that shortlist; this is not a blanket clinical contraindication. Readthrough is deferred on RNA, product/function and exposure grounds, not a claim that UGA can never insert leucine. A published reporter study reported Leu with Leu/Ile measurement ambiguity. The oncology review treats tumor protection as a testable risk, without claiming HSP70 inhibition proves arimoclomol would protect this tumor. The contribution is the genotype-to-protein-to-experiment decision, not demonstrated efficacy or a claim of global priority.

### Please state whether you relied solely on publicly available drug/target databases and literature (preferred), or whether proprietary data sources were also used.

Drug/mechanism selection used public literature and regulatory information, informed by the supplied case phenotype and genetic hypothesis. No proprietary drug screen, unpublished sponsor data or patient response measurements were used.

### Please describe the public data sources used in detail (e.g., DrugBank, ChEMBL, ClinicalTrials.gov, OMIM, PubMed literature search).

Primary literature includes BUBR1 stability/function experiments (Suijkerbuijk et al., PMID 20516114), arimoclomol Gaucher-cell experiments (Fog et al., PMID 30497978), mitotic/aneuploidy-model NAC studies, readthrough amino-acid studies (Xue et al., PMID 28575328), and HSP70/antioxidant oncology studies. Regulatory prescribing information supplies indication and exposure context. The [report](jeweled_track2_report_submission.md) links these primary studies and exposure/oncology appendices. Their findings are not treated as measured benefit for this genotype. We do not imply use of every database named as an example in the form.

### Please describe any proprietary data sources used in detail.

No proprietary drug/target dataset was used. Organizer-supplied information provided the case context. Laboratory/model access is unestablished, and no patient-derived drug-response data are available.

### How did you characterize the variant's mechanism (e.g., loss-of-function vs. gain-of-function, pathway disrupted, downstream biological consequence)?

Leu737Ter has predicted loss of function/nonsense-mediated decay in the relevant transcript. This is not measured endogenous RNA loss. N1002K remains VUS: its abundance, turnover, localization and mitotic function are unmeasured. Trans is favored under the recessive model but unobserved.

We therefore test whether residual missense protein is both unstable and functionally rescuable. More protein alone does not establish rescue. A second route asks whether downstream stress protection improves function while segregation errors persist. The causal EPCR of 0.8 is not a probability that either drug works.

### Please provide an estimate of time or effort spent on this analysis (if possible).

No consolidated person-hour or cost ledger was measured. Work comprised computational interpretation, primary-source review, exposure/oncology appraisal and report/pitch preparation. No laboratory execution has occurred. Assay-capacity scenarios are not measured throughput or power.

### Provide a method abstract (up to 500 words). Please include strengths and limitations of your approach as well as methodological detail.

We propose restoring usable BUBR1 protein as the lead therapeutic strategy for a BUB1B stop-plus-missense hypothesis. The key experiment tests whether restoring N1002K abundance restores mitotic function. This distinguishes a potentially rescuable dosage defect from a function-specific defect and determines whether a proteostasis drug deserves further testing.

Arimoclomol is the first experimental drug priority; NAC is a comparator for downstream stress protection. Neither has demonstrated efficacy for this genotype. Qualification of a consented model, genetic correction controls and exposure measurement precede drug interpretation. First-mitosis errors are assessed jointly with first-anaphase yield without a scored visible error per pretreatment starting cell, retaining death, non-entry, slippage and tracking loss. Marker changes or selective survival alone cannot establish functional rescue.

A parallel oncology counter-screen tests tumor protection and interactions with a confirmed active regimen. Credible adverse effects block advancement; imprecise results remain inconclusive. Meaningful benefit/adverse-effect margins and confirmation rules must be set before confirmatory experiments.

The approach connects allele mechanism to a practical drug-testing decision. Its main limitations are absent patient-derived functional/response evidence, unknown phase, uncertain useful exposure and unestablished laboratory access. It is a falsifiable repurposing proposal, not a treatment recommendation.
