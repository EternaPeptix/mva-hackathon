# Methods description form — answers

**Team: jeweled. Research submission draft, revised 21 September 2026; not submitted.** The 22 question headings reproduce the original methods workbook. The [current submission index](README.md) identifies the active CSV, reports and pitch; dated earlier answers are archival.

The private repository is [EternaPeptix/mva-hackathon](https://github.com/EternaPeptix/mva-hackathon). It contains reproducible source, public/synthetic validation and selected report summaries. Raw case reads, genome-scale call payloads and internal working history remain excluded. Public HG002 replay is separate from authorized case analysis.

**AI use:** Anthropic Claude Code and OpenAI Codex, personal subscriptions, training on the user's content disabled (user-confirmed). Exact historical model/plan names are unverified. Computation used owned machines; candidate-level variants, evidence summaries and report text were reviewed with hosted assistants. This does not imply every historical targeted public-service query ran locally.

**Narration:** Cartesia Sonic 3.6 API, stock Lauren voice, approved by the user. Cartesia plan/tier and training/retention settings are unverified; Claude/Codex settings do not establish them. The service did not return an immutable model snapshot. The approved 180-second movie is attached to a private GitHub release and predates the new Track 2 lead. The recording remains unchanged; media work and Vimeo are on hold.

# Track 1 (Variant Prediction) — Model 1

### Team name

jeweled

### Model number (fill out this form once per model/approach submitted; up to 6)

1

### Please describe your model/approach in detail.

We ask whether the supplied MVA case is explained by a recessive BUB1B pair. Eight FASTQ checksums were verified, four sequencing lanes were aligned to GRCh38 on Mac Studios, and GATK 4.6.1.0 plus bcftools 1.24 were run independently on the same reads. Their completed whole-reference outputs were reconciled to the original call sections. Coding inventories cover all 22 autosomes; fixed-ploidy sex-chromosome, mitochondrial, copy-number and structural analyses have separate scopes.

We preserved original multiallelic genotypes, missing fields and annotation lineages. LIRICAL/Exomiser rankings used the supplied VCF; no ranking rerun on the final four-lane native callsets is claimed. Exomiser's BUB1B recessive entry ranked 2 (score 0.4877) with primary predictors, disappeared under REVEL-only strict filtering, and ranked 3 (0.1806) under REVEL-only permissive filtering. LIRICAL's MVA1 rank was 4/8,525 with seven phenotypes, 38 without rhabdomyosarcoma, 10 with rhabdomyosarcoma alone, and 1 with genotype integration (model probability 0.9575393, not EPCR). The exact pair is less robust to predictor choice than the gene's prominence. A separate final native allele-specific coding inventory used frequency filters without pathogenicity predictors. Regional read and joint-sequence review then tested alternative explanations. The [case-evidence appendix](appendices/track1-case-evidence.md) supplies the complete comparison, selected aggregates and authorized execution map.

Both callers find the same two BUB1B protein-altering alleles and no rare protein-altering/canonical-splice allele in ten comparison genes. The current prediction is that BUB1B pair at subjective EPCR **0.8**, with phase unknown and Asn1002Lys still VUS. The [Track 1 report](jeweled_track1_report_submission.md) and [interpretation appendix](appendices/causal-interpretation.md) explain the inference. The previously nominated candidate and supplied solved-case context were known; this was not blinded discovery.

### Please state whether the submission file is the automated output of your computational approach (preferred), or if it has undergone downstream manual review and curation.

The CSV underwent manual review with AI assistance. The existing builder validates an explicit candidate decision, GRCh38 reference alleles, phase field and probability rationale. Selection and EPCR are research judgments, not calibrated automated outputs. The new one-row export replaces both the historical 0.95 and the later 0.35 judgments; earlier files remain archival.

### Please describe any downstream manual review in detail.

We separated allele presence, phase, allele classification and case-level causality. Leu737Ter has an external ClinVar Pathogenic/Likely pathogenic two-star assertion. Our declared points-based reconstruction applies PVS1 Very Strong plus PM2 Supporting (9 points, LP); Asn1002Lys receives PM2 Supporting plus PM3 Supporting (2 points, VUS). Both rarity judgments are explicitly accepted under general recessive-disease guidance. PM3 uses the independent Labcorp Pathogenic partner assertion and the unknown-phase allowance; no circular upgrade is made. PP4 is a disclosed sensitivity, not rejected merely because MVA is genetically heterogeneous. [Criteria and primary guidance](appendices/causal-interpretation.md).

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
| LIRICAL 2.4.1 / HPO 2026-09-01 / Ensembl transcript database | Phenotype-only and genotype-integrated comparisons; exact public-resource identities in the [selected evidence](../results/track1-selected-case-evidence.json). |
| gnomAD v4.1, ClinGen and targeted ClinVar records | Exact-allele frequencies, gene validity and classification guidance. Missing observations are retained. |
| ClinVar GRCh38 snapshot dated 2026-09-13 | Frozen exact-allele and assertion-lineage review. |
| ACMG SF v3.3 / frozen ClinGen panel | Secondary inventory of 84 genes and reference-only coding targets. |
| Public nuclear marker resource | A fixed 37,267-marker model check, without a validated biological detection limit. |
| GIAB HG002 chr20 / v4.2.1 small-variant truth and confident BED | Public SNV/indel benchmarking. |
| GIAB HG002 structural truth / v5.0q, chr20 subset | Public structural benchmarking; [truth composition and residuals](appendices/public-sv-residuals.md). |

HG002 reads are labelled nominally 35×. A subsequent descriptive read-depth measurement at MAPQ≥20/BQ≥20 found mean depth 35.75× across chr20 and 36.97× in the official confident regions; 99.25% of confident-region bases had at least 20× depth. Zero-depth positions were included; overlapping mates both count, deletion bases do not, and default unmapped/secondary/QC-fail/duplicate exclusions apply. This is not either caller's internal read-filter definition. RTG Tools 3.13 vcfeval compared diploid, raw calls in those confident regions, with full-chromosome context. Original GATK SNV/indel recall was 99.3748%/97.3902%, precision 98.9548%/96.0168%; bcftools recall was 99.3593%/91.5115%, precision 98.9512%/91.1827%. Later retained-multiallelic results are a separate representation sensitivity. We do not claim an additional hap.py cross-check. DELLY and SvABA matched 69/802 and 171/802 structural truth alleles; insertions and repeat-associated misses limit negative case interpretation.

N1002K scores are REVEL **0.472** and AlphaMissense **0.9229**. The adopted REVEL policy gives neither PP3 nor BP4; AlphaMissense remains a disclosed alternative, not stacked evidence. Historical provenance identifies REVEL v1.3 and the primary AlphaMissense hg38 release, without an exact AlphaMissense release date. LIRICAL's maximum-predictor rule and shared Exomiser resources do not provide independent confirmation. Scores were not multiplied into EPCR, and no new predictive model was trained. [Sensitivity details and resource-date limitations](appendices/track1-case-evidence.md#ranking-and-predictor-sensitivity).

### Please describe any proprietary data used in detail.

The nonpublic inputs were organizer-supplied phenotype, variants and sequencing reads. No parental or tumor specimen data, patient-derived experimental observations or drug-response measurements were available.

### Is your approach able to output proposed pairs of compound heterozygous candidate variants, or only single candidate variants?

Pairs. The primary row contains heterozygous BUB1B p.(Leu737Ter) and p.(Asn1002Lys). Trans is the favored inference under the proposed recessive explanation, but phase remains explicitly unknown in the export. Phenotype or model ranking is not an observed segregation result.

### How did you handle secondary or incidental findings, if any? (See the finding_type column in your submission.)

The frozen 84-gene secondary review retained exact coding and noncoding matches, original genotypes, missing quality fields and assertion lineage. One called SCN5A allele per caller was individually Benign/Likely benign and selected only through a containing-haplotype assertion. It is not nominated. No exact called allele entered a clean aggregate P/LP stratum in the frozen resources. [Inventory](appendices/secondary-findings.md).

Across the 393,880-base coding union, 99.30% had depth at least 20 at MQ20/BQ20 and 98.59% at MQ60/BQ20; the latter had 2,834 zero-depth bases in five genes. These are descriptive coverage results, not a clinical negative. [Coverage and local gaps](appendices/secondary-depth.md).

The separate nuclear-mixture estimator returned a boundary value of 0.0 after 37,149/37,267 markers passed coverage gates. This does not establish purity, mosaicism sensitivity or variant phase. [Model limitations](appendices/nuclear-mixture.md).

### Please provide an estimate of run time and cost for your approach (if possible).

No validated end-to-end runtime or consolidated compute, electricity, subscription or labor cost was measured for the case analysis. A separate fresh public HG002 chr20 replay completed in 134.3 minutes, including downloads, setup, both callers, comparison and audit, with 5.31 GiB peak sampled process-group memory. It exactly reproduced the historical primary counts and audit partitions. The older 116-minute measurement covered calling only. These public timings are not the case's full analysis time or guarantees for other machines. Mac Studios performed the current genome work. Prospective resource ceilings and assay-capacity estimates are planning values.

### Provide a method abstract (up to 500 words). Please include strengths and limitations of your method as well as methodological detail.

We investigated a supplied MVA case using verified sequencing reads, GRCh38 alignment, two variant callers and allele-specific coding inventories. Phenotype/genotype ranking was supplemented by predictor sensitivities, population review, read phasing and joint-sequence reconstruction. Copy-number, structural, mitochondrial and secondary-finding analyses retain distinct scopes.

The leading hypothesis is BUB1B p.(Leu737Ter) plus p.(Asn1002Lys), at subjective EPCR 0.8. Both callers recover two rare BUB1B coding alleles and no rare protein-altering/canonical-splice alternative in ten comparison genes. Phase is unobserved: the 10,911 bp gap exceeds the original 288 bp phase-block maximum and later local blocks. Leu737Ter has a two-star external P/LP assertion; the declared research reconstruction is LP, while Asn1002Lys remains VUS. The supplied solved-case context informs the inference, so this is not blinded discovery.

Strengths include preservation of original genotypes, investigation of alternative sequence representations and correction of a misleading mosaic signal. Public HG002 chr20 benchmarking uses GIAB confident regions and RTG vcfeval, with separate SNV/indel results. Synthetic tests check implementation behavior. Neither establishes case-wide clinical sensitivity.

Limitations include one sample, absent parental and functional evidence, incomplete sensitivity to regulatory/structural/repeat mechanisms and no measured treatment response. Residual hypotheses are stated rather than hidden behind missing phase. The next decisive studies are segregation and a qualified N1002K abundance/function assay.

# Track 2 (Drug Repositioning)

### Team name

jeweled

### Please describe your approach in detail (how you went from variant/mechanism to candidate medication(s)).

We nominate everolimus as the lead drug hypothesis and niacin (nicotinic acid) as the primary comparator, tested separately rather than in combination. They test mTORC1 normalization versus restoration of usable BUBR1. The BUB1B stop-plus-missense hypothesis has unresolved phase and endogenous protein consequences; genotype alone does not establish either drug-responsive mechanism.

Everolimus has an established drug-to-target link and indication-specific pediatric pharmacology. In a different BUB1B allele model, gastrocnemius showed increased p70 S6 kinase and 4EBP1 phosphorylation; the signaling blot used two mice per genotype. This supports a pathway hypothesis, not causal or pharmacological rescue. We narrowly prioritize the implementable medication hypothesis; prioritizing proximal causal restoration would put niacin first. The two-route set is more defensible than claiming a universal drug order.

Niacin has human muscle NAD engagement evidence. Separate NMN/SIRT2 experiments establish a BUBR1 stability mechanism, and genetic SIRT2 elevation improved selected outcomes in BubR1-hypomorphic mice. These are distinct interventions, not a demonstrated niacin-to-functional-BUBR1 chain. Qualify actual precursor utilization, SIRT2/BUBR1 effects and useful function. Arimoclomol is a restoration-specific reserve: HSP90 dependence does not establish its proposed chaperone rescue, and useful exposure is unresolved.

Dasatinib is a separate reserve requiring harmful senescent burden, selective susceptibility and improved absolute output. Genetic clearance supplies the strongest same-gene organismal functional precedent but does not establish dasatinib activity; D+Q evidence cannot be transferred to monotherapy and quercetin eligibility is unresolved.

Compare absolute muscle output, progenitor renewal and chromosome fidelity in qualified affected/corrected models. Everolimus need not repair the initiating mitotic defect to improve downstream function, but must not worsen fidelity or exhaust renewal. A claim of niacin-mediated BUBR1 restoration requires the corresponding protein/function evidence. The [report](jeweled_track2_report_submission.md) and [pipeline](appendices/track2-evidence-pipeline.md) specify source boundaries and prospective decisions.

### Was candidate identification automated (e.g., via a drug-target database query or LLM-assisted search), or primarily manual literature review and expert curation?

Candidate identification and review were AI-assisted, using Claude Code and Codex for source searches, comparison, mechanistic critique and drafting. Selection is a curated research judgment. No validated automated drug-ranking model, numerical therapeutic score or independent clinical expert sign-off is claimed.

### Please describe any manual literature review or expert curation in detail.

We compared actual drug-to-target evidence, disease causality, tissue relevance, useful exposure and pediatric feasibility using the same standard for each candidate. No independent clinical expert sign-off or quantitative efficacy ranking is claimed. Everolimus wins the implementation comparison; niacin has the stronger proximal restoration rationale. Neither has demonstrated rescue of the nominated genotype.

We retained adverse and null evidence without assigning it excessive weight. The public mouse RNA Hallmark null does not rebut measured substrate phosphorylation. Preliminary rapamycin nonbenefit in a different BubR1 model lacks treatment details; everolimus-associated BUBR1 loss occurred in late-passage nonmutant Wharton-jelly cells. These findings define functional and renewal checks, not proof of MVA harm.

Niacin is not NMN and neither inherits genetic SIRT2 outcomes. NAD normalization alone is insufficient; a pediatric mitochondrial-myopathy case progressed despite biochemical correction. SIRT2 has multiple substrates, so dependence on SIRT2 alone does not prove BUBR1 dependence. Human-myotube utilization experiments retain panel-specific statistics and exposure limitations.

Arimoclomol's other-protein effects do not establish BUBR1 rescue. Failure of additional N1002K at matched abundance weakens quantity-only rescue without excluding improved folding/function per molecule. Genetic INK-ATTAC benefit, marker-level D+Q patent evidence and Ercc1 drug results remain separate. The alternatives table records the strongest case, weakness and condition for reconsidering each route.

The sequence reconstruction and 1,536 assumed population runs inform conditional products and endpoint design. The nine-mouse RNA reanalysis found no Hallmark program significant after global correction across 150 tests. Neither computation ranks drugs or demonstrates biological rescue.

### Please state whether you relied solely on publicly available drug/target databases and literature (preferred), or whether proprietary data sources were also used.

Drug/mechanism selection used public literature and regulatory information, informed by the supplied case phenotype and genetic hypothesis. No proprietary drug screen, unpublished sponsor data or patient response measurements were used.

### Please describe the public data sources used in detail (e.g., DrugBank, ChEMBL, ClinicalTrials.gov, OMIM, PubMed literature search).

Primary sources include [Suijkerbuijk 2010, MVA expression/function](https://pubmed.ncbi.nlm.nih.gov/20516114/); [North 2014, NAD/SIRT2/BUBR1](https://pmc.ncbi.nlm.nih.gov/articles/PMC4194088/); [Pirinen 2020, human muscle niacin/NAD](https://pubmed.ncbi.nlm.nih.gov/32386566/); [Membrez 2024, NAPRT-dependent utilization in human myotubes](https://pmc.ncbi.nlm.nih.gov/articles/PMC10963276/); and [Aaltio 2024, pediatric niacin counterexample](https://pubmed.ncbi.nlm.nih.gov/39173541/). [Everolimus labeling](https://dailymed.nlm.nih.gov/dailymed/fda/fdaDrugXsl.cfm?setid=2150f73a-179b-4afc-b8ce-67c85cc72f04), [niacin labeling](https://dailymed.nlm.nih.gov/dailymed/fda/fdaDrugXsl.cfm?setid=dd60fe81-3d54-4a29-9502-601be89c575e) and [arimoclomol labeling](https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid=5feffc0e-453d-47fa-91dd-38d4952309bc&version=5) supply approved-use and pharmacology context, not MVA efficacy.

Alternatives use [Baker 2011, INK-ATTAC](https://www.nature.com/articles/nature10600), [D+Q patent Example 11](https://patents.google.com/patent/WO2015116735A1/en), [Zhu 2015, pharmacological senolysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC4531078/), [Sieben 2020, BubR1 muscle signaling](https://www.jci.org/articles/view/126863/version/2), [Baker 2013, rapamycin observation](https://pmc.ncbi.nlm.nih.gov/articles/PMC3785294/) and [Goutas 2023, autophagy/BUBR1](https://pmc.ncbi.nlm.nih.gov/articles/PMC10149375/). Additional primary sources and rejected routes appear in the [proposal table](jeweled_track2_report_submission.md#other-routes-considered-and-not-selected-as-the-lead).

The executed [RNA analysis](../workflow/track2-rna/README.md) uses GSE134781 and MSigDB 2025.1.Mm orthology-mapped Hallmarks, with complete results, sensitivities and environment. Public GRCh38/Ensembl sequence supports the [genome-informed analysis](../workflow/track2-model/README.md). No use of every example database named in the form is implied.

### Please describe any proprietary data sources used in detail.

No proprietary drug/target dataset was used. Organizer-supplied information provided the case context. Laboratory/model access is unestablished, and no patient-derived drug-response data are available.

### How did you characterize the variant's mechanism (e.g., loss-of-function vs. gain-of-function, pathway disrupted, downstream biological consequence)?

Leu737Ter predicts truncation and possible NMD; reference reconstruction does not measure endogenous decay. N1002K remains VUS, with unmeasured abundance, turnover, localization and mitotic function. Trans is favored under the recessive model but unobserved. Canonical cis translation masks the downstream substitution; this does not establish half-normal usable BUBR1.

Genetic correction and controlled actual-allele expression tests determine whether the nominated genotype produces a functional defect and whether surviving protein is useful, inert or interfering. A separate branch asks whether harmful mTORC1 excess contributes to muscle dysfunction; the different-allele mouse result cannot establish that state in this case.

If myogenic mTORC1 excess is established, a [conditional IR-trafficking follow-up](appendices/track2-evidence-pipeline.md#conditional-follow-up-insulin-receptor-trafficking) compares endogenous receptor uptake and signaling duration in affected/corrected cells. The [Choi 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5560052/) and [2019](https://www.nature.com/articles/s41467-019-09318-3) studies motivate this proposed mechanism; it is unestablished for N1002K and muscle disease. Receptor dependence, trafficking mediation and harmful downstream activity are separate questions. This adds no drug or completed experiment.

Niacin tests a NAD/SIRT2-to-usable-BUBR1 hypothesis; everolimus tests partial mTORC1 normalization. Restoration requires relevant protein/function evidence, whereas downstream benefit may occur without correcting the initiating segregation defect. Both must improve useful absolute output without unacceptable renewal loss or worsened fidelity. A benefit without the proposed causal dependence is retained but reclassified. Genomic EPCR 0.8 is not a probability of drug benefit.

Either functional cohort can establish meaningful benefit if the other cohort,
renewable progenitor reserve and fidelity meet predeclared preservation margins.
The four drug-by-benefit-cohort claims form one confirmation family, with one
fixed regimen per drug. Each claim requires all benefit/preservation components;
the largest component P value is tested using Holm correction across four claims.
Pilot selection and independent confirmation remain separate; the [analysis rule](appendices/track2-evidence-pipeline.md#confirmatory-benefit-rule)
defines the scope and prevents additional doses or time points from becoming
unadjusted opportunities for success.

### Please provide an estimate of time or effort spent on this analysis (if possible).

No consolidated person-hour or cost ledger was measured. Work comprised computational interpretation and sequence reconstruction, assumption-based population simulations on a Mac Studio, primary-source comparison, exposure/oncology appraisal and report/pitch preparation. Public mouse count reanalysis also ran on a Mac Studio. No laboratory execution has occurred. Assay-capacity scenarios and simulation runtimes are not measured laboratory throughput or power.

### Provide a method abstract (up to 500 words). Please include strengths and limitations of your approach as well as methodological detail.

We propose everolimus as the lead drug hypothesis and niacin (nicotinic acid) as the primary comparator for a BUB1B stop-plus-missense case, tested separately. The question is whether mTORC1 normalization or restoration of usable BUBR1 can preserve productive muscle function. Everolimus has established target pharmacology and a small same-gene muscle signaling precedent, but no relevant rescue result. Niacin has human muscle NAD engagement; separate NMN/SIRT2 experiments support BUBR1 stabilization without establishing niacin rescue. Prioritizing implementation favors everolimus; prioritizing proximal causality favors niacin. This selects a focused hypothesis set, not demonstrated drug superiority.

Model qualification establishes endogenous RNA/protein, phase where possible, allele function, signaling and reproducible tissue defects with genetic controls. Abundance-matched allele studies distinguish quantity from function per molecule. Niacin experiments qualify precursor utilization and NAD/SIRT2/BUBR1 dependence. Partial genetic mTORC1 normalization tests pathway causality before attributing drug benefit to that pathway. Drug experiments use measured, pharmacologically justified exposures; genetic results cannot be transferred to a medication.

The shared comparison measures absolute muscle output, productive progenitor yield and chromosome fidelity. Death, arrest, division time and tracking loss prevent selective loss or slowed cycling from masquerading as rescue. A downstream intervention need not repair the primary mitotic defect, but must not worsen fidelity or exhaust renewal. Mature-muscle function and progenitor replacement are measured separately; either can supply benefit if the other is preserved. Four joint benefit/preservation claims undergo multiplicity correction. Independent confirmation, withdrawal follow-up, prespecified margins and relevant tumor/regimen review follow useful functional results.

Arimoclomol is a restoration reserve requiring exposure-qualified mutant-function rescue. Dasatinib is a senolytic reserve requiring harmful burden, selective susceptibility and increased absolute output; it cannot inherit genetic clearance or D+Q combination results.

Reference-sequence reconstruction and assumed population scenarios inform interpretation and endpoints. Public nine-mouse RNA reanalysis found no globally significant Hallmark program; it does not rebut the published phosphorylation result or rank drugs. Strengths are orthogonal hypotheses and tests that separate target success from drug failure. Limitations include unknown phase and allele function, absent drug-response evidence, uncertain useful exposure and unestablished laboratory access. No clinical efficacy is claimed.
