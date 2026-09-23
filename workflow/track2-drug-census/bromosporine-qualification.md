# Bromosporine and anti-PD-1: primary-source qualification

Assessed 2026-09-23. Read-only repository context: `725447082e515c090c2139c89b23a398e4f25079`. This is an original evidence assessment, not an experimental replication or clinical recommendation. Raw source assets and derived inspection files remain in this external assessment directory.

**Disposition:** retain bromosporine as a broad bromodomain research compound with reported pancreatic-tumor monotherapy and combination activity. Keep the anti-PD-1 combination as a distinct preclinical cancer finding. Do not promote either to BUB1B restoration, muscle compensation, MVA benefit, or a selective catalytic PCAF inhibitor. The reported tumor benefit is substantial enough to preserve; uncertainty concerns mechanism, formal synergy and transfer to constitutional disease.

## What was actually read

| Source | Reading depth and limits |
|---|---|
| Zhang et al., Nature Communications 2026, PMID 41663394, DOI 10.1038/s41467-026-69311-5 | Full article XML and 17-page PDF obtained. Results, discussion, relevant experimental/statistical methods and main figure legends reviewed. Main Figures 4 and 5 inspected as original publisher images. Figures 1–3 provide pathway context; their image panels were not independently audited. |
| Zhang supplementary information | All 16 pages text-extracted and read, including all 10 figure legends and three supplementary tables. Figures S7, S8 and S10 visually inspected on PDF pages 8, 9 and 11; S8 legend continues on p10. Supplementary data-table sequences were inspected only for whether BUB1B was assayed, not sequence-validated. |
| Zhang source data | Original XLSX acquired unchanged. Selected values/labels for Figures 4, 5, S8 and S10 inspected; tumor-weight means below are simple arithmetic summaries of published values. No new hypothesis tests, curve fitting or synergy analysis performed. Embedded original blot images and all unrelated sheets were not exhaustively audited. |
| Zhang reporting summary / transparent peer review | Reporting-summary pp2–3 visually inspected for sample-size, exclusion, randomization and blinding declarations. Peer-review text searched and relevant PCAF, immune-cell and dose-schedule exchanges read. These are author/reviewer statements, not independent validation. |
| Picaud et al., Science Advances 2016, PMID 27757418, DOI 10.1126/sciadv.1600760 | Full text XML acquired; biochemical/selectivity, cellular engagement, cancer-response results and relevant methods read. Table 1 inspected numerically; Figures 1–2 visually inspected. Entire 32-page supplement acquired; targeted S2 selectivity pages 5–6 inspected visually and S3/S4/S15/Table S3 descriptions inspected as text. Not a complete audit of every transcriptomic supplement. |
| SGC bromosporine page | Full official HTML downloaded and identity, property, selectivity and cell-assay sections read. |
| Moustakim et al., 2017, PMC5412877 | Full text XML acquired; targeted domain-architecture and bromosporine-to-L-Moses context read. No transfer of L-Moses efficacy or PBMC tolerance to bromosporine. |
| Li et al., Molecular Cancer 2024, PMID 38711083 | Full text XML acquired; targeted intervention methods and H3K18la–TTK/BUB1B results/legends read to resolve the citation link. Figures and supplement were not independently reviewed; this is a citation-trace check, not completion of that paper's drug qualification. |
| Manufacturer identities / later PCAF claim | Official indexed Bio X Cell and MCE content consulted. Direct downloads returned 403; indexed excerpts retained. El-Shershaby et al. 2021, PMID 34126285, DOI 10.1016/j.bmc.2021.116266, was available only as indexed publisher abstract/introduction excerpts; its assay methods were not qualified. It receives no full-text or mechanism-completion credit. |

The main [Zhang article](https://www.nature.com/articles/s41467-026-69311-5), [supplement](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-69311-5/MediaObjects/41467_2026_69311_MOESM2_ESM.pdf), and [source data](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-69311-5/MediaObjects/41467_2026_69311_MOESM4_ESM.xlsx) anchor the experiment-level claims below. Download URLs and SHA-256 hashes are in [bromosporine-qualification.json](bromosporine-qualification.json) and the external retrieval manifest.

## Identity and target attribution

The pancreatic paper used **bromosporine, MCE HY-15815**. The official supplier identifies CAS **1619994-69-2**, formula **C17H20N6O4S**, molecular weight **404.44**. SGC reports the same formula/weight and InChIKey **UYBRROMMFMPJAN-UHFFFAOYSA-N**, with a triazolopyridazine scaffold. This is consistent with the parent bromosporine/BSP characterized by Picaud, not L-Moses, its inactive stereoisomer D-Moses, a biotin-tagged BSP derivative or an unspecified “BSP” abbreviation. The study's actual lot, purity certificate and analytical identity were not provided in the inspected sources. [SGC](https://www.thesgc.org/chemical-probes/bromosporine-bsp), [MCE identity](https://www.medchemexpress.com/Bromosporine.html).

The strongest primary characterization establishes **acetyl-lysine reader-domain binding with broad bromodomain coverage**:

| Recombinant bromodomain | Bromosporine ITC KD, nM (mean ± SD) |
|---|---:|
| PCAF / KAT2B | 4,762 ± 459 |
| BRD4(1) | 41.8 ± 2.8 |
| BRD4(2) | 39.7 ± 2.2 |
| BRD9 | 41.7 ± 3.8 |
| CECR2 | 8.0 ± 1.0 |
| TAF1(2) | 16.6 ± 2.7 |

These are **binding affinities**, not PCAF catalytic acetylation/lactylation IC50 values. Picaud also presents a 42-bromodomain BLI screen, acetyl-lysine-pocket crystal structures, cellular BRD4/BRD9 displacement and BET-like transcriptional responses in leukemia cells. The much weaker PCAF-domain affinity and multiple strong alternate targets make selective PCAF attribution untenable from bromosporine exposure alone. The 5–20 μM cancer-cell concentrations in Zhang overlap broad target pharmacology; biochemical affinity does not quantify intracellular occupancy. [Primary probe paper, Table 1 and Figs1–2](https://pmc.ncbi.nlm.nih.gov/articles/PMC5061470/).

PCAF has a bromodomain and a distinct HAT catalytic domain. Its bromodomain can influence chromatin recruitment without being the catalytic active site. The subsequent L-Moses paper explicitly develops a selective bromodomain tool to separate these functions. [Moustakim primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5412877/).

Zhang supplies real **PCAF-pathway evidence**: Pcaf knockdown lowers H3K18la/CXCL1, and Fig4F shows purified PCAF-dependent histone lactylation with appropriate component omissions. But **Fig4F does not contain bromosporine**. Its lactyl-CoA docking is also not a bromosporine-binding assay. Cellular H3K18la lowering, a second compound (embelin), and phenocopy by gene knockdown do not establish that bromosporine directly or selectively inhibits the catalytic PCAF reaction. No PCAF target-engagement assay in PDAC, drug-resistant PCAF rescue, selective bromodomain-versus-catalytic comparison or purified-PCAF bromosporine inhibition curve was identified.

An auxiliary unresolved claim must remain visible: the MCE page and the indexed 2021 El-Shershaby abstract give a bromosporine “PCAF IC50” of 2.10 μM. The available excerpt repeatedly frames the project as bromodomain inhibition, but the exact assay construct/readout was not accessible. It cannot be recoded as catalytic lactyltransferase inhibition, nor can this assessment claim that every possible catalytic effect has been disproved. The qualified conclusion is narrower: **broad bromodomain inhibition is established; selective catalytic PCAF inhibition and its responsibility for the PDAC phenotype remain unproven.** [2021 publisher record](https://www.sciencedirect.com/science/article/abs/pii/S0968089621002741).

## Actual administered agents, models and controls

All tumor experiments here used **male C57BL/6J mice aged 6–8 weeks**, not MVA models. Anti-PD-1 was **Bio X Cell BE0146**, mapped by the official manufacturer sheet to **RMP1-14, rat IgG2a κ, anti-mouse PD-1**. The control was **BE0089**, mapped to **2A3, rat IgG2a κ anti-trinitrophenol**. This is not pembrolizumab or nivolumab. Clone mapping is manufacturer-based; the paper supplies catalog numbers rather than the clone name. [BE0146](https://bioxcell.com/pub/media/tds/BE0146-TDS.pdf), [BE0089](https://bioxcell.com/pub/media/tds/BE0089-TDS.pdf).

| Experiment | Intervention and comparator | Endpoint/sample |
|---|---|---|
| Cancer cells, Fig4C/G/H | Bromosporine 0, 5, 20 μM in PANC1, BxPC3, KPC and PANC02 for H3K18la blots; PANC1/BxPC3 for CXCL1 RNA/protein. Embelin 0, 10, 20 μM is a separate corroborating cell experiment in Fig4D. Exact incubation duration for these bromosporine panels was not located. | H3K18la lowering; CXCL1 RNA/secretion lowering. CXCL1 assays n=3 independent experiments. No direct cell-viability normalization establishes that all transcriptional changes are separable from cytostasis/toxicity. |
| Orthotopic monotherapy, Fig4J–R | 2×10^6 KPC-luc cells/pancreas; bromosporine 30 mg/kg i.p. daily on days7–13, versus vehicle. | Tumor weight/immune/CXCL1 n=5 per group, tumor blots n=4. |
| Subcutaneous monotherapy, FigS8J–U | KPC-luc cells; 3×10^6 cells/right flank per model methods; bromosporine 30 mg/kg i.p. daily starting day9, 11 doses per methods, versus vehicle. | n=6 per group. The schematic marks through day20 with harvest day21, while body-weight plots run day9–19; exact schedule alignment remains a minor reporting gap. |
| Subcutaneous combination, Fig5A–G | 6×10^6 PANC02 cells/right flank. Four arms: DMSO+IgG2a; bromosporine+IgG2a; DMSO+anti-PD-1; bromosporine+anti-PD-1. Bromosporine 30 mg/kg i.p. 3×/week and antibody 200 μg/mouse i.p. 2×/week; two five-day cycles separated by two rest days. Fig5A: antibodies days13,16,20,23; small molecule days14,15,17,21,22,24; harvest day25. | n=6 tumors per arm. Fig5E caption says n=3 mice per arm; source workbook has nine values per arm and appears to use antibody-positive-cell percentages. Their sampling/aggregation unit requires clarification. |
| Orthotopic combination, Fig5H–Q | 2×10^6 KPC-luc cells/pancreas; same four arms/doses. Fig5H: antibodies days6,9,13,16; bromosporine days7,8,10,14,15,17; imaging day12; endpoint harvest day18. | 21 mice per arm: 6 for tissue/flow endpoints and 15 for survival. Source flux table contains 15 values per arm. |

Vehicle is labelled DMSO in combination panels; its full injected formulation/concentration was not located. Administration details are recorded for source qualification, not proposed use. The authors describe the two rest days as a precaution against potential combined-treatment toxicity in peer-review correspondence, not as proof that toxicity occurred.

## Benefit and null results

The following tumor-weight values are **descriptive means computed from the supplied source workbook**, in grams. They are not independent replication and the groups belong to different experiments.

| Model/figure | Vehicle/control | Bromosporine alone | Anti-PD-1 alone | Combination |
|---|---:|---:|---:|---:|
| Orthotopic monotherapy, Fig4L, n=5 | 0.980 | 0.518 | — | — |
| Subcutaneous monotherapy, FigS8M, n=6 | 0.610 | 0.259 | — | — |
| Subcutaneous combination study, Fig5D, n=6 | 0.465 | 0.286 | 0.286 | 0.180 |
| Orthotopic combination study, Fig5J, n=6 | 1.100 | 0.773 | 0.415 | 0.233 |

Locators: `Figure 4!A35:B39`; `Figure S8!B67:G68`; `Figure 5!C19:H22` and `C60:H63`. Published comparisons support lower monotherapy tumor weight (Fig4L p=0.0093; S8M p=0.0485) and additional combination benefit versus either single agent (Fig5D p=0.0365/0.0489; Fig5J p<0.0001 versus bromosporine and p=0.0005 versus anti-PD-1). Fig5P/Q also reports lower bioluminescence and longer survival in the combination arm. The same survival plot shows benefit in each monotherapy arm relative to vehicle; the additional combination benefit must remain attributed to the combination. These results are not evidence of a cure or a clinical effect.

**Formal synergy was not established by the reported analysis.** The abstract/discussion use “synergistic,” while the PANC02 Results paragraph calls tumor-growth inhibition “additive.” Figure5 uses pairwise two-tailed unpaired t-tests and log-rank survival tests. No specified additive null, factorial interaction test, Bliss/Loewe/HSA analysis, combination index or dose matrix was identified. Superiority to each single agent supports combination benefit; it does not by itself establish a supra-additive interaction. No independent synergy test was run here.

Immune and biochemical results remain endpoint-specific:

- Bromosporine monotherapy lowers tumor H3K18la, serum CXCL1 and tumor-infiltrating neutrophils, while increasing CD8 and granzyme-B-positive CD8 fractions. The orthotopic **PD-1-positive CD8 fraction is unchanged** (Fig4R p=0.4041).
- Combination treatment lowers tumor burden, neutrophil fractions and increases granzyme-B-positive CD8 cells. But **orthotopic total CD8 infiltration is not significantly higher than either monotherapy** (Fig5L p=0.6784 and 0.8616 for the displayed comparisons). **PD-1-positive CD8 is not different from anti-PD-1 alone** (Fig5N p=0.9871); bromosporine alone is also not different from vehicle there (p=0.1696). Calling a PD-1-positive fraction “exhausted” is the authors' marker-based interpretation, not a comprehensive exhaustion-function assay.
- Subcutaneous serum CXCL1 is lowered by bromosporine alone and by combination; adding anti-PD-1 does not significantly improve it over bromosporine alone (Fig5G p=0.2814). Anti-PD-1 alone does not alter serum CXCL1 versus control (p=0.9453). Subcutaneous neutrophil percentage is not further reduced versus bromosporine alone (Fig5E p=0.6229).
- H3K18ac declines less than H3K18la in treated cultured cells; H4K8ac is unchanged in the tested gradients. Tumor H3K18ac/H4K8ac and promoter enrichment for these acetylation marks do not show the claimed lactylation change (S8G–I). This supports context-specific differential readouts, not biochemical selectivity for one PCAF acyltransferase reaction.

## Normal-cell evidence, quality limits and MVA transfer

Body weight is reported unaffected in the monotherapy experiment. Subcutaneous combination S10 reports no significant group difference in body/kidney/spleen weights. These are limited gross-tolerance observations in tumor-bearing mice; no normal-organ functional assays, detailed systemic toxicity assessment, nonmalignant BUB1B-deficient cells, muscle physiology or chromosome-segregation rescue were identified.

S7B compares baseline histone marks in two “normal” cell lanes (N1/N2) and cancer lines. S7H's direct immune-cell comparison uses **2DG**, not bromosporine. S6's preserved neutrophil survival/ROS/degranulation concerns conditioned medium from **2DG- or LDH-knockdown-treated cancer cells**, not direct bromosporine exposure. These cannot be relabelled as bromosporine normal-cell safety. SGC's HeLa cytotoxicity note concerns cancer cells; the L-Moses PBMC result concerns a different compound.

The reporting summary says specimens were randomly distributed, declares blinding “not applicable,” and gives no explicit power calculation. The animal figure legends identify one experiment, despite generic statements about independent repeats elsewhere. Pairwise testing across many endpoints/time points is reported without a clearly specified multiplicity correction. Fig5E's nine source values versus three mice and S8R's five available source values versus a six-mouse caption remain unresolved sampling/aggregation issues. These limit precision without erasing the observed tumor findings.

**BUB1B connection:** Zhang's Discussion cites Li2024 for H3K18la-dependent transcription of TTK/BUB1B in PDAC. Zhang's qualified bromosporine experiments report CXCL1/H3K18la/immune endpoints, not BUB1B abundance, BUBR1 function, checkpoint restoration or MVA rescue. The traced Li primary text used glycolysis inhibitors/genetic manipulations and linked reduction of H3K18la with lower TTK/BUB1B in cancer cells; it contains no bromosporine treatment. That citation provides a possible **suppression-direction cancer connection**, not a demonstrated restoration route. It also does not prove bromosporine lowers BUBR1 in normal tissue. [Li2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11071201/).

Constitutional BUB1B deficiency, residual-protein preservation, developmental disease and muscle function are untested. The data therefore do not change the provisional ordering of niacin restoration/apitegromab compensation or everolimus as a conditional probe; none becomes proven MVA therapy through this finding.

## Pending gaps and bounded follow-up disposition

1. Qualify the exact PCAF assay underlying the later 2.10 μM claim from a lawful full text/supplement; keep bromodomain binding separate from acetyltransferase and lactyltransferase activity.
2. Resolve PDAC cellular target engagement, catalytic inhibition and on-target rescue; neither PCAF biology nor a shared phenotype selects PCAF from bromosporine's many targets.
3. Clarify Fig5E/S8R sample aggregation, cell exposure times, vehicle composition and monotherapy schedule alignment. Existing statistical reporting supports benefit but leaves interaction strength and some experimental units unresolved.
4. Keep normal-cell/immune selectivity, systemic exposure and constitutional BUB1B-deficient tissue effects explicitly untested. No MVA promotion is warranted without functional benefit and preservation evidence in the relevant genotype/tissue.

**Compound-specific census action:** retain separate entries for (a) bromosporine monotherapy cancer activity, (b) bromosporine + RMP1-14 preclinical combination benefit and (c) unproven PCAF-mechanism attribution. Reject the stronger labels “selective catalytic PCAF inhibitor,” “formally demonstrated synergy,” “BUB1B restoration” and “MVA rescue” for this qualified evidence set. This is a source-based disposition, not an exhaustive literature-negative claim.


Integration note: this delegated source qualification contributes zero additional outer-corpus completion credit. It is not a clean-state review of the whole submission or biological replication.

[Structured qualification](bromosporine-qualification.json).
