# Interpretation of the proposed BUB1B pair

**21 September 2026.** We assign a **subjective research EPCR of 0.8** to the exact GRCh38 pair **chr15:40209701 T>G, p.(Leu737Ter), and chr15:40220612 T>G, p.(Asn1002Lys)** as the explanation of the primary phenotype. This is a current evidence judgment, not a calibrated posterior, clinical diagnosis, classifier score or probability of treatment benefit. The two allele classifications and the pair's causal probability answer different questions.

## Clinical context and independence

The organizer describes a child living with MVA and a clinically solved answer: the [overview](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/about.py) states that the answer is clinically confirmed, and the [rules](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/rules.py) refer to NHS-validated causal variant(s). That supplied disease context informs the review. The [FAQ](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/faq.py) and overview explain partial credit conditionally for a compound-heterozygous pair; this scoring example does not unambiguously establish the hidden answer's allelic structure. We do not use it as independent evidence that this case must have a causal pair.

Those statements do **not** publicly confirm our exact two positions, supply their segregation results or disclose the underlying clinical laboratory evidence. We have not examined the private gold-standard answer or scored a real submission. Earlier candidate awareness prevents a blinded-discovery claim. Public filename hints are not used as genetic evidence.

## What the completed analysis supports

The two native caller views agree that both proposed BUB1B variants are heterozygous. Filtered REF/ALT fragment counts are 21/26 and 15/12, respectively, with alternate support on both strands. The annotated BUB1B coding union has minimum mapped depth 24 in both reviewed MQ20/BQ20 and MQ60/BQ20 views. These findings support allele presence; repeated analysis of the same reads is not independent biological replication. The [case-evidence map](track1-case-evidence.md) links the [selected numerical aggregates](../../results/track1-selected-case-evidence.json), original analysis scopes and reconstruction entry points.

The completed chromosome 1–22 inventory contains **exactly two BUB1B protein-altering/canonical-splice alleles**, both the proposed candidates. The separate BUB1B gene-plus-5 kb census finds no additional emitted exon/CDS or screened splice-edge overlap. The following comparison agrees between GATK and bcftools:

| Gene | Protein-altering/canonical-splice inventory alleles | Alleles with observed retained maximum AF below 1% |
|---|---:|---:|
| BUB1B | 2 | 2 |
| TRIP13 | 0 | 0 |
| CEP57 | 0 | 0 |
| CENATAC | 1 | 0 |
| MAD2L1BP | 0 | 0 |
| SLF2 | 1 | 0 |
| SMC5 | 1 | 0 |
| TRIM37 | 0 | 0 |
| BUB1 | 0 | 0 |
| BUB3 | 0 | 0 |
| MAD2L1 | 0 | 0 |

Thus **none of these ten comparison genes has a rare allele in the reviewed coding/canonical-splice inventory**. Their three retained coding alleles have population maximum frequencies of approximately 4.47%, 61.42% and 93.39%, respectively; frequency is not unresolved for these rows. The inventory includes stop/frameshift, canonical splice donor/acceptor, missense, in-frame indel, start-lost and stop-lost effects on supported protein-coding transcripts. It preserves native genotype, filter and missing quality fields, without a newly imposed PASS/GQ admission threshold. Its frequency view is the retained gnomAD-prefix AC/AN maximum, not a complete population census. The 1% comparison is a discovery view, not our PM2 cutoff.

This is meaningful negative evidence within the tested scope. It is not proof that every comparison gene is normal: broader splice effects, deep-intronic/regulatory variants, difficult structural events and low-fraction variants remain incompletely assessed. None of these ten genes has a zero-depth coding base in the reviewed depth views, but depth alone does not establish a reference genotype or clinical sensitivity. The historical count of “523 variants excluding the differential” is not reused; its corrected source count was 518, and the newer exact-allele inventory provides the stronger comparison.

Completed broader reviews have not supplied a better-supported causal explanation. ADAMTS1's adjacent isolated frameshift annotations combine into a supported net 21-base in-frame sequence; FAM136A's combined consequence is missense. FANCD2 retains locus ambiguity. The AFDN/KIF25/FRMD1 structural event is technically supported without an established explanatory effect. These findings weaken particular alternatives without proving all alternatives benign or closing the differential.

## Phase and the retracted mosaic signal

**Phase remains unobserved.** The vendor VCF's 304,329 annotated phase blocks have a maximum observed variant-position span of **288 bp**; the two targets are **10,911 bp apart**. That vendor annotation statistic is not a universal limit of short-read phasing. The earlier two-lane 100 kb WhatsHap review recovered local blocks up to **2,327 bp**, still without linking the targets. The final four-lane bounded review also leaves a 6,769 bp marker gap and an unknown pair relation in every run; its longest local block is 77 bp. These scopes must not be conflated.

These unsuccessful links do not demonstrate cis phase. Trans is the favored configuration **under the recessive causal model**, but has not been observed; parental segregation could falsify that model. Unknown phase is neither a measured 50:50 probability nor repeated negative evidence each time the same reads fail to bridge the gap.

The earlier apparent chromosome 20/21/22 mosaic signal was **retracted**: the skew localized to pericentromeric/acrocentric sequence and disappeared as a chromosome-wide excess after masking. We do not count it as support for MVA. Its rejection also does not refute the organizer's clinical diagnosis; this bulk short-read analysis does not provide a validated cell-by-cell aneuploidy exclusion.

## Allele classifications and the adopted evidence policy

For Leu737Ter, the external [ClinVar VCV000533901.9](https://www.ncbi.nlm.nih.gov/clinvar/variation/533901/) classification is **Pathogenic/Likely pathogenic, two stars**: Labcorp classifies it Pathogenic for MVA1 and GeneDx Likely pathogenic with condition not provided. This is nonconflicting multiple-submitter evidence, not unanimous Pathogenic or an expert-panel classification. Our narrower research reconstruction is reported alongside it, not substituted for it.

We use the [Tavtigian points framework](https://pmc.ncbi.nlm.nih.gov/articles/PMC8011844/): Supporting, Moderate, Strong and Very Strong contribute 1, 2, 4 and 8 points; LP is 6–9 and Pathogenic is at least 10. The adopted table is:

| Criterion | Leu737Ter | Asn1002Lys |
|---|---|---|
| PVS1 | Very Strong, 8: relevant-transcript premature stop with predicted NMD | Not applicable |
| PM2 | Supporting, 1: adopted sufficient-rarity judgment for recessive disease | Supporting, 1: adopted sufficient-rarity judgment for recessive disease |
| PM3 | Not applied: unphased VUS partner | Supporting, 1: unphased co-occurrence with an independently accepted Pathogenic partner |
| PP4 | Not applied; Supporting sensitivity below | Not applied; Supporting sensitivity below |
| **Applied total and working category** | **9, Likely pathogenic** | **2, VUS** |

**Rarity policy.** We judge **both** alleles sufficiently rare in this recessive-disease context under [general ClinGen PM2 guidance](https://www.clinicalgenome.org/site/assets/files/5182/pm2_-_svi_recommendation_-_approved_sept2020.pdf). This is a disclosed research judgment, not a BUB1B VCEP threshold or a calibrated population likelihood ratio. The stop's exome AC/AN is 115/1,461,846 (AF approximately 7.87e-5), largest observed ancestry AF 9.98214e-5 and FAF95 maximum 8.395e-5. The missense is a high-quality singleton, 1/1,461,878, with ancestry AF 8.99276e-7. A null genome response for the missense is not a verified zero count. Recessive alleles need not be absent from controls. The older Labcorp population-quality warning is retained as historical context; current empty filter/flag lists and quality summaries do not reproduce it. We do not use zero homozygotes as an additional pathogenicity criterion.

**Loss of function.** PVS1 follows the existing MANE/exon assessment and predicted NMD, with consideration of transcript relevance under [ClinGen PVS1 guidance](https://pubmed.ncbi.nlm.nih.gov/30192042/). No patient RNA decay, protein deficiency or alternative-transcript rescue was measured. The stop's evidence does not depend on calling N1002K pathogenic.

**Co-occurrence.** [ClinGen PM3 v1.0](https://www.clinicalgenome.org/site/assets/files/3717/svi_proposal_for_pm3_criterion_-_version_1.pdf) awards 0.5 **case points** for unknown phase with a Pathogenic partner, reaching Supporting; an LP partner would give only 0.25 and would not reach Supporting in this single case. We explicitly accept Labcorp's independently derived Pathogenic assertion for this partner decision. Both alleles meet our adopted rarity judgment, and the organizer supplies affected-case context. This is noncircular despite our separate research LP reconstruction. PM3 Supporting contributes one final classification point; it neither proves trans nor requires PP4 to be awarded. We do not feed the missense's PM3 evidence back into the stop's classification.

**Phenotype sensitivity.** MVA with rhabdomyosarcoma favors BUB1B, consistent with [Hanks et al.](https://pubmed.ncbi.nlm.nih.gov/15475955/) and the broader [BUB1B/TRIP13 tumor evidence](https://pmc.ncbi.nlm.nih.gov/articles/PMC5493194/). Current [ClinGen PP4 guidance](https://pmc.ncbi.nlm.nih.gov/articles/PMC10806742/) allows genetic heterogeneity using appropriately defined diagnostic yield. Our omission is therefore not a blanket “several genes” rule, an insistence on obtaining our own cytogenetics, or a claim that now-completed structural studies are unfinished. We have not established a robust, method-matched BUB1B yield for this precise phenotype sufficient to quantify PP4. If Supporting PP4 is accepted, the stop becomes **10 points, Pathogenic**, and the missense **3 points, VUS**. That is an explicit sensitivity, not the applied table.

No qualifying functional assay, new segregation evidence or measured N1002K effect is added. The retained REVEL policy is unchanged: **0.472** is indeterminate, supplying neither PP3 nor BP4. AlphaMissense **0.9229** is a disclosed alternative, not selected or multiplied to obtain a preferred classification. The [ranking sensitivities](track1-case-evidence.md#ranking-and-predictor-sensitivity) show that strict REVEL-only filtering removes the recessive pair even while BUB1B remains prominent. A VUS category does not set a numerical upper bound on this pair's probability of explaining a particular clinically solved case.

## Why EPCR is 0.8, and what could change it

The supplied solved-MVA context, supported alleles, stop-plus-rare-missense mechanism, exact ten-gene comparison and failure of reviewed alternatives to displace the pair together support the leading hypothesis. Gene-level support alone does not settle the exact partner. N1002K is the only observed rare protein-altering second BUB1B allele, has credible read support, and the locus census adds no other emitted coding or screened splice-edge candidate. We give this observed-pair explanation greater weight than one requiring an additional unseen hit. That is a qualitative judgment: the review does not quantify sensitivity to every noncoding or complex hit, and the missense may still be incidental or insufficiently damaging.

We retain a rounded subjective EPCR of 0.8 under this weighting, without converting rarity, prioritizer scores or ACMG points into a calibrated posterior. The earlier 0.35 and 0.95 were also judgments. A downward revision need not await a new cis, normal-function or competing-cause observation: removal of overstated support can itself justify reassessment. Our changed weighting is not empirical proof that 0.35 was wrong or that 0.8 is uniquely warranted. No new segregation/function result or probability confidence band is claimed.

Residual possibilities are described without invented percentage allocations:

- **A different second BUB1B hit:** N1002K could be incidental while a missed noncoding or complex allele supplies the causal partner, especially if the observed pair is cis. A [reported BUB1B intronic insertion](https://pubmed.ncbi.nlm.nih.gov/38102195/) illustrates this residual class; it is not evidence that such an insertion occurs here.
- **Insufficient missense effect:** N1002K may not produce the damaging effect required by the proposed mechanism. An informative validated functional result, including appropriate genetic rescue, would materially change this judgment.
- **Another or blended mechanism:** remaining structural, regulatory, low-fraction or incompletely captured chromosome-instability causes could explain the phenotype.

Parental phase and an informative N1002K functional assessment would be more discriminating than another same-read caller. The subjective EPCR remains separate from allele-classification points and does not establish that raising BUBR1 abundance, correcting a cellular surrogate or administering a proposed drug would benefit the child.
