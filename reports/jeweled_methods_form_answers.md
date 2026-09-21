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

We selected everolimus as the experimental lead for normalizing mTORC1 if excess signaling drives a demonstrated muscle defect. The BUB1B stop-plus-missense hypothesis motivates a genotype-qualified model. A published mouse model carrying a different BubR1 truncation showed muscle mTORC1 activation, but did not test rapalog rescue. This supplies a pathway hypothesis, not evidence that the nominated human genotype has the same defect. Elevated signaling could be compensatory.

Selection is a qualitative comparison of mechanistic proximity, intervention evidence, drug-to-target support, exposure feasibility and decisive failure tests. Everolimus has an established mTORC1 target relationship, labeled TSC experience in younger children and a dispersible formulation. Sirolimus also has an oral solution and therapeutic monitoring. The within-class choice is pragmatic, not demonstrated MVA superiority. Rapamycin/sirolimus, rather than everolimus, is the actual compound in the cited fly chromosome-instability study. Arimoclomol remains an alternative if N1002K has a functionally rescuable abundance defect; genetic senescent-cell clearance has a stronger organism-level intervention precedent than the current senolytic drug bridge. The [evidence pipeline](appendices/track2-evidence-pipeline.md) and [unelected-route table](jeweled_track2_report_submission.md#other-routes-considered-and-not-selected-as-the-lead) record the comparison.

We also completed an independent secondary analysis of GSE134781 public mouse muscle counts on a Mac Studio: nine mice, 50 Hallmark programs and three planned contrasts, using TMM/voom and CAMERA with estimated inter-gene correlation. No program met global FDR <0.05; the mTORC1 X753–WT result was weakly downward relative to background (P=0.219; FDR=0.983). It supplies no positive transcriptional corroboration, does not measure kinase activity and does not rank alternatives. Full results and sensitivities are retained in the [public RNA analysis](../workflow/track2-rna/README.md). This limits our evidence claim and keeps direct pathway qualification decisive.

A qualified model must first establish a reproducible muscle defect and test whether partial genetic normalization of excess signaling improves function. Complete pathway ablation would not answer this question. Matched mature myotubes provide the functional test; separate progenitor differentiation/regeneration experiments test whether apparent benefit sacrifices tissue renewal. Everolimus experiments then require measured exposure, target engagement, functional improvement and preserved absolute viable-cell/muscle output and growth. A rhabdomyosarcoma counter-screen addresses tumor protection and interaction with any confirmed active regimen. No patient-derived or drug-response experiment has run. [Track 2 report](jeweled_track2_report_submission.md).

### Was candidate identification automated (e.g., via a drug-target database query or LLM-assisted search), or primarily manual literature review and expert curation?

Candidate identification and review were AI-assisted, using Claude Code and Codex for source searches, comparison, mechanistic critique and drafting. Selection is a curated research judgment. No validated automated drug-ranking model, numerical therapeutic score or independent clinical expert sign-off is claimed.

### Please describe any manual literature review or expert curation in detail.

We checked each bridge against the actual allele, protein, species, tissue, developmental stage, intervention and endpoint. The BubR1 truncation mouse supplies muscle-pathway evidence; it does not establish a drug response or harmful mTOR activation in this case. The fly rapamycin result concerns preservation of stem-cell numbers and does not establish restoration of progeny or brain size. Neither result is relabelled as everolimus efficacy. Normal muscle development and regeneration also require mTOR signaling, so lower phosphorylation or smaller cells cannot serve as the benefit endpoint.

Alternatives receive the same bar. Expression rescue of selected BUBR1 mutants does not show N1002K rescue, and HSP90 dependence does not establish arimoclomol/HSP70 rescue. Genetic removal of p16-positive cells is not a validated senolytic drug. NAC and elamipretide have overlapping downstream effects, so differences in drug response alone would not distinguish redox from membrane/bioenergetic mechanisms. Checkpoint inhibition opposes the intended nonmalignant rescue mechanism and is excluded from that shortlist; this is not a blanket clinical contraindication. Readthrough remains dependent on RNA availability, product identity/function and exposure, rather than an absolute claim that UGA cannot insert leucine.

The imported sequence analysis checks the nominated coding consequences against public reference sequence. Its separate population model explores how continuing missegregation, growth, death, senescence and differentiation affect absolute output and composition. The 1,536 deterministic runs span declared assumptions and two arbitrary time horizons; they are neither independent biological samples nor fitted drug-response predictions. Their role is to expose misleading surrogate endpoints and specify measurements, not choose a medication.

### Please state whether you relied solely on publicly available drug/target databases and literature (preferred), or whether proprietary data sources were also used.

Drug/mechanism selection used public literature and regulatory information, informed by the supplied case phenotype and genetic hypothesis. No proprietary drug screen, unpublished sponsor data or patient response measurements were used.

### Please describe the public data sources used in detail (e.g., DrugBank, ChEMBL, ClinicalTrials.gov, OMIM, PubMed literature search).

Primary sources include [Sieben et al., BubR1 truncation and muscle signaling](https://www.jci.org/articles/view/126863/version/2); the [rapamycin chromosome-instability fly study](https://doi.org/10.1038/s41467-026-70521-0); [Baker et al., genetic senescent-cell clearance](https://pmc.ncbi.nlm.nih.gov/articles/PMC3468323/); [Suijkerbuijk et al., BUBR1 expression/function](https://pubmed.ncbi.nlm.nih.gov/20516114/); and [Xu et al., PP2A-B56 recruitment in MVA fibroblasts](https://pmc.ncbi.nlm.nih.gov/articles/PMC3654266/). The latter are mechanistic comparators, not evidence for the selected drug.

Our executed [public-count analysis](../workflow/track2-rna/README.md) uses
[GSE134781](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134781) and the
MSigDB 2025.1.Mm orthology-mapped Hallmark collection. The source identities,
software environment, all gene/set results and sensitivity analyses are retained.
The [2026 randomized sirolimus exercise study](https://doi.org/10.1002/jcsm.70274)
is included as counterevidence to assuming that intermittent exposure preserves
function; it is not a study of MVA.

The [everolimus prescribing information](https://www.novartis.com/us-en/product/pi/pdf/afinitor.pdf) and [sirolimus prescribing information](https://labeling.pfizer.com/showlabeling.aspx?id=139) supply approved-use, formulation, pharmacology and safety context. Primary studies of [mTOR-dependent muscle regeneration](https://pmc.ncbi.nlm.nih.gov/articles/PMC2793064/) and [muscle Raptor ablation](https://doi.org/10.1016/j.cmet.2008.10.002) motivate separate mature-muscle and progenitor safeguards; complete genetic loss is not equated with partial drug inhibition. Public GRCh38 and Ensembl transcript sequence support the [genome-informed analysis](../workflow/track2-model/README.md). The governing [evidence pipeline](appendices/track2-evidence-pipeline.md) links claims to source scope and stopping rules. We do not imply use of every database named as an example in the form.

### Please describe any proprietary data sources used in detail.

No proprietary drug/target dataset was used. Organizer-supplied information provided the case context. Laboratory/model access is unestablished, and no patient-derived drug-response data are available.

### How did you characterize the variant's mechanism (e.g., loss-of-function vs. gain-of-function, pathway disrupted, downstream biological consequence)?

Leu737Ter predicts truncation and possible nonsense-mediated decay in the relevant transcript. Reference-sequence reconstruction does not measure endogenous RNA loss or protein abundance. N1002K remains VUS: its abundance, turnover, localization and mitotic function are unmeasured. Trans is favored under the recessive model but unobserved. These uncertainties remain in the model-qualification stage.

The selected downstream hypothesis is that a BUB1B-dependent muscle defect includes harmful excess mTORC1 signaling. Matched nutrient conditions and differentiation stage, genetic correction controls and partial pathway titration must distinguish that hypothesis from normal growth signaling or compensation. Only functional benefit at an appropriate degree of normalization supports the mechanism. Everolimus is not proposed to correct the chromosomes, and genomic EPCR 0.8 is not a probability of drug benefit.

### Please provide an estimate of time or effort spent on this analysis (if possible).

No consolidated person-hour or cost ledger was measured. Work comprised computational interpretation and sequence reconstruction, assumption-based population simulations on a Mac Studio, primary-source comparison, exposure/oncology appraisal and report/pitch preparation. Public mouse count reanalysis also ran on a Mac Studio. No laboratory execution has occurred. Assay-capacity scenarios and simulation runtimes are not measured laboratory throughput or power.

### Provide a method abstract (up to 500 words). Please include strengths and limitations of your approach as well as methodological detail.

We propose everolimus as a conditional experimental lead for muscle dysfunction associated with a BUB1B stop-plus-missense hypothesis. A different BubR1 truncation mouse model links the gene to excess muscle mTORC1 signaling, without demonstrating rapalog rescue. Everolimus supplies an established drug-to-target relationship and practical formulation/approved-use experience. Its selection is qualitative and does not establish superiority to sirolimus or benefit for this genotype.

Model qualification must establish the endogenous allele consequences, muscle dysfunction and signaling state. Genetic correction and partial pathway normalization test whether excess signaling drives the defect rather than compensating for it. Matched mature myotubes assess muscle function; separate progenitor differentiation/regeneration experiments detect a cost to growth and renewal. Complete pathway ablation would not test the proposed normalization mechanism.

A secondary analysis of the nine public mouse muscle samples did not corroborate an upregulated mTORC1 Hallmark program (global FDR 0.983 for X753–WT). This is not a kinase-activity measurement or independent biological replication; the unfavorable result limits the convergence claimed for the lead.

Drug studies require measured everolimus exposure, target engagement, improved function and preserved absolute viable output, alongside a rhabdomyosarcoma counter-screen and interactions with any confirmed active regimen. Lower pathway markers, selective survival or improved fractions alone cannot establish rescue. Continuous and pulse/washout exposure must be compared without presuming protective timing; sirolimus is the reference compound. Reference-sequence reconstruction and assumption-based population scenarios inform endpoint design without fitting a genotype-to-drug response model.

Alternatives include abundance rescue, senescent-cell clearance and downstream stress protection; each retains its own missing mechanistic or pharmacological bridge. Strengths are explicit source boundaries, falsifiable stage-specific experiments and attention to absolute output. Limitations include uncertain phase and missense function, absent case-specific mTOR evidence, unknown useful exposure and unestablished laboratory access. No drug experiment, chromosome correction or clinical efficacy is claimed.
