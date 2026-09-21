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

We select usable BUBR1 restoration as the primary research route, with niacin (nicotinic acid) as the first exploratory drug nomination and arimoclomol as its prespecified challenger. The nominated BUB1B stop-plus-missense genotype has unresolved phase and endogenous protein consequences. Selected human MVA mutants have expression-rescue precedents; neither nominated drug has demonstrated N1002K rescue.

Niacin has human muscle NAD target-engagement evidence. Human-myotube NA/NAPRT evidence supports testing precursor utilization; panel-specific limitations are retained in the pipeline. Precursor utilization must be qualified in the actual model before interpreting nonresponse. Separate NMN/SIRT2 experiments establish a BUBR1 stability mechanism, and genetic SIRT2 elevation improved selected outcomes in BubR1-hypomorphic mice. These are distinct interventions, not a continuous demonstrated niacin-to-functional-BUBR1 chain. Arimoclomol has pediatric approved-use experience and other-protein processing evidence, but HSP90 dependence does not establish its HSP70-to-BUBR1 bridge. We prioritize the niacin experiment without claiming drug superiority.

Selection considers same-gene functional intervention evidence, actual drug implementation, exposure, approved-medication eligibility and a decisive functional test. Senescent-cell clearance ranks highest for same-gene organismal functional evidence; the gap is an eligible, selective drug implementation. D+Q eligibility is unresolved, and dasatinib alone cannot inherit combination evidence. Dasatinib remains a conditional reserve if relevant senescent populations impair function and respond selectively. Arimoclomol is chosen to challenge the same restoration question, not to declare it universally better than senolysis.

Everolimus is now a conditional pathway alternative. Public muscle RNA reanalysis supplies no positive mTORC1 Hallmark corroboration; a preliminary unpublished rapamycin nonbenefit observation and a separate late-passage-cell BUBR1 reduction with everolimus further weaken sole-lead status. None establishes MVA harm. The [report](jeweled_track2_report_submission.md) and [evidence pipeline](appendices/track2-evidence-pipeline.md) connect each source to its scope and proposed decisive experiment.

### Was candidate identification automated (e.g., via a drug-target database query or LLM-assisted search), or primarily manual literature review and expert curation?

Candidate identification and review were AI-assisted, using Claude Code and Codex for source searches, comparison, mechanistic critique and drafting. Selection is a curated research judgment. No validated automated drug-ranking model, numerical therapeutic score or independent clinical expert sign-off is claimed.

### Please describe any manual literature review or expert curation in detail.

We checked allele, species, tissue, developmental stage, intervention and endpoint, explicitly separating pathway evidence from named-drug rescue. Niacin is not NMN, and neither inherits genetic SIRT2 lifespan results. NAD normalization is not functional rescue; a pediatric mitochondrial-myopathy case progressed despite biochemical correction. SIRT2 also has multiple substrates. The niacin proposal therefore requires evidence tying functional benefit to usable BUBR1, not NAD alone.

Arimoclomol's other-protein results do not establish BUBR1 rescue. An abundance-matched N1002K failure weakens pure quantity rescue, but does not exclude improved folding/function per molecule. Conversely, a niacin/arimoclomol response difference cannot by itself distinguish these mechanisms. Checkpoint inhibition conflicts with the proposed nonmalignant fidelity objective; it is not described as a universal clinical contraindication.

Genetic INK-ATTAC clearance improved selected BubR1-hypomorph tissue outcomes. D+Q same-model patent evidence concerns markers, not functional rescue; peer-reviewed D+Q progeroid functional results concern Ercc1. Everolimus-associated BUBR1 reduction was measured in late-passage Wharton-jelly MSCs, not the patient genotype. Preliminary rapamycin nonbenefit lacks treatment details. All remain correctly bounded counterevidence.

The sequence analysis reconstructs conditional coding products. The separate population model's 1,536 deterministic runs use declared assumptions and arbitrary time horizons, not patient-calibrated drug effects. A secondary public mouse RNA analysis retained 12,989 genes across nine mice and found no Hallmark program significant after global correction across 150 tests. Neither computation ranks drugs or establishes biological rescue.

### Please state whether you relied solely on publicly available drug/target databases and literature (preferred), or whether proprietary data sources were also used.

Drug/mechanism selection used public literature and regulatory information, informed by the supplied case phenotype and genetic hypothesis. No proprietary drug screen, unpublished sponsor data or patient response measurements were used.

### Please describe the public data sources used in detail (e.g., DrugBank, ChEMBL, ClinicalTrials.gov, OMIM, PubMed literature search).

Primary sources include [Suijkerbuijk 2010, MVA expression/function](https://pubmed.ncbi.nlm.nih.gov/20516114/); [North 2014, NAD/SIRT2/BUBR1](https://pmc.ncbi.nlm.nih.gov/articles/PMC4194088/); [Pirinen 2020, human muscle niacin/NAD](https://pubmed.ncbi.nlm.nih.gov/32386566/); [Membrez 2024, NAPRT-dependent utilization in human myotubes](https://pmc.ncbi.nlm.nih.gov/articles/PMC10963276/); and [Aaltio 2024, pediatric niacin counterexample](https://pubmed.ncbi.nlm.nih.gov/39173541/). [Niacin labeling](https://dailymed.nlm.nih.gov/dailymed/fda/fdaDrugXsl.cfm?setid=dd60fe81-3d54-4a29-9502-601be89c575e) and [arimoclomol labeling](https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid=5feffc0e-453d-47fa-91dd-38d4952309bc&version=5) supply approved-use and pharmacology context, not MVA efficacy.

Alternatives use [Baker 2011, INK-ATTAC](https://www.nature.com/articles/nature10600), [D+Q patent Example 11](https://patents.google.com/patent/WO2015116735A1/en), [Zhu 2015, pharmacological senolysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC4531078/), [Sieben 2020, BubR1 muscle signaling](https://www.jci.org/articles/view/126863/version/2), [Baker 2013, rapamycin observation](https://pmc.ncbi.nlm.nih.gov/articles/PMC3785294/) and [Goutas 2023, autophagy/BUBR1](https://pmc.ncbi.nlm.nih.gov/articles/PMC10149375/). Additional primary sources and rejected routes appear in the [proposal table](jeweled_track2_report_submission.md#other-routes-considered-and-not-selected-as-the-lead).

The executed [RNA analysis](../workflow/track2-rna/README.md) uses GSE134781 and MSigDB 2025.1.Mm orthology-mapped Hallmarks, with complete results, sensitivities and environment. Public GRCh38/Ensembl sequence supports the [genome-informed analysis](../workflow/track2-model/README.md). No use of every example database named in the form is implied.

### Please describe any proprietary data sources used in detail.

No proprietary drug/target dataset was used. Organizer-supplied information provided the case context. Laboratory/model access is unestablished, and no patient-derived drug-response data are available.

### How did you characterize the variant's mechanism (e.g., loss-of-function vs. gain-of-function, pathway disrupted, downstream biological consequence)?

Leu737Ter predicts truncation and possible NMD; reference reconstruction does not measure endogenous decay. N1002K remains VUS, with unmeasured abundance, turnover, localization and mitotic function. Trans is favored under the recessive model but unobserved. Canonical cis translation masks the downstream substitution; this does not establish half-normal usable BUBR1.

The central hypothesis is that a surviving product has function that can be preserved or improved. Genetic correction and controlled WT/actual-allele expression tests qualify that hypothesis before assigning drug action. A surviving truncation can be helpful, inert or interfering; abundance alone cannot determine the sign. Joint fidelity and absolute-output measurements distinguish useful restoration from slower cycling or selective loss. If a drug helps without BUBR1 dependence, its mechanism must be reclassified. The genomic EPCR 0.8 is not a probability of drug benefit.

### Please provide an estimate of time or effort spent on this analysis (if possible).

No consolidated person-hour or cost ledger was measured. Work comprised computational interpretation and sequence reconstruction, assumption-based population simulations on a Mac Studio, primary-source comparison, exposure/oncology appraisal and report/pitch preparation. Public mouse count reanalysis also ran on a Mac Studio. No laboratory execution has occurred. Assay-capacity scenarios and simulation runtimes are not measured laboratory throughput or power.

### Provide a method abstract (up to 500 words). Please include strengths and limitations of your approach as well as methodological detail.

We propose testing usable BUBR1 restoration in a BUB1B stop-plus-missense hypothesis, prioritizing niacin (nicotinic acid) with arimoclomol as a prespecified challenger. Niacin has human muscle NAD engagement evidence; separate NMN/SIRT2 experiments support BUBR1 stabilization. Genetic SIRT2 results do not establish niacin rescue. Arimoclomol has pediatric approved-use and other-protein processing evidence, but no demonstrated BUBR1 implementation. Neither drug is a proven winner.

Model qualification establishes phase where possible, endogenous RNA/protein, allele function and reproducible segregation/myogenic defects with genetic controls. Abundance-matched WT and variant experiments distinguish quantity from function per molecule. Vehicle, niacin and arimoclomol are tested separately at measured, pharmacologically justified exposures. NMN and controlled SIRT2 perturbation interrogate the niacin mechanism without substituting for the nominated drug.

Joint functional gates require improved chromosome fidelity with preserved productive progenitor output, followed by absolute muscle function per starting cohort. Death, arrest, division time and tracking loss prevent selection or slowed cycling from masquerading as rescue. Mature-muscle function and progenitor replacement are measured separately. Independent confirmation, withdrawal follow-up, prespecified margins and relevant tumor/regimen assessment are required.

Senolysis supplies the strongest same-gene organismal functional precedent but lacks a validated eligible implementation. Dasatinib is a conditional reserve; D+Q combination results cannot establish monotherapy activity. Everolimus is a pathway-qualified alternative, weakened by absent positive RNA corroboration, underreported rapamycin nonbenefit and a different-cell-state BUBR1 depletion observation.

Reference-sequence reconstruction and assumed population scenarios inform interpretation and endpoint design. A public nine-mouse count reanalysis found no globally significant Hallmark program; it does not measure kinase activity or rank drugs. Strengths are explicit intervention boundaries and tests that can reject target or compound separately. Limitations are unknown phase and allele function, absent drug-response evidence, uncertain useful exposure and unestablished laboratory access. No chromosome correction or clinical efficacy is claimed.
