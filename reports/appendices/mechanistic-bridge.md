> Package copy: local links/navigation are adapted for this finite review package. Labels identify underlying files that are not included. Scientific text is retained; approved status/source metadata changes are recorded separately; original and packaged hashes are recorded in PROVENANCE.

> **Retained dated method detail.** This appendix preserves the original analysis or proposal at its stated date. The current [Track 1 report](../jeweled_track1_report_submission.md), [Track 2 report](../jeweled_track2_report_submission.md) and [methods form](../jeweled_methods_form_answers.md) govern interpretation; see the [revision and evidence scope](../PROVENANCE.md#evidence-scope). Historical delivery or ranking statements are not current status.

# Track 2: mechanism, exposure and prior-work review

19 September 2026. This review supports the
[team report — archival document not included](../PROVENANCE.md#archival-material). It changes candidate prioritisation and
proposed experiments; it supplies no patient-cell or clinical drug result.

## What the BUBR1 manuscript establishes

In Suijkerbuijk et al., Figure 3C tests HSP90 inhibition in engineered U2OS cells expressing
I909T or L1012P. Geldanamycin reduces mutant abundance, and proteasome inhibition prevents
that clearance. A separate forced-expression experiment in Figure 4D restores the nocodazole
response. Both relevant panels are described as representative experiments; their legends
do not supply an independent-replicate count. Counts from neighbouring panels must not be
transferred to them. Neither panel demonstrates arimoclomol treatment, HSP70-mediated rescue,
N1002K behaviour or correction of spontaneous chromosome errors.
[Author manuscript, Figures 3C and 4D](https://www.hubrecht.eu/app/uploads/2015/08/Kops_Research_2010_Suijkerbuijk_Molecular-Causes-for-BUBR1-Dysfunction-in-the-Human-Cancer-Predisposition-Syndrome-Mosaic-Variegated-Aneuploidy.pdf).

The inference for this project is conditional: establish the endogenous allele defect before
testing a drug, and separate protein quantity, checkpoint response and segregation. A chaperone
family name is insufficient evidence of a shared rescue mechanism.

## Exposure check

A Gaucher fibroblast study reports HSPA1A induction at 100-400 micromolar arimoclomol,
GCase activity experiments at 50-800 micromolar, and emerging cytostasis at 400-800 micromolar.
It also examines longer treatment at lower concentrations. These experiments concern a
lysosomal enzyme, not BUBR1. The authors could not conclusively assign the effect to an
individual HSP70-family member through their knockdown experiment.
[Primary study, results and Figure 2](https://doi.org/10.1016/j.ebiom.2018.11.037).

The FDA's **August 2024 preapproval advisory briefing**, pages 42-43, questions the relevance
of several NPC-cell findings at 200-400 micromolar against a reported clinical maximum of
about 6.7 micromolar. Figure 10 contains eight NPC cell lines, with panel annotations of
n=3-5; these are not eight patients in a treatment trial. This is an advisory assessment,
not the final approval decision or current label.
[FDA briefing](https://www.fda.gov/media/180469/download).

The current MIPLYFFA label's Section 12.3 gives the following distinct populations and
statistics. Version 5, published 20 March 2026, was verified through DailyMed's history API.

| Label quantity | Reported concentration | Conversion to micromolar |
|---|---:|---:|
| Pediatric NPC, estimated steady-state serum trough, mean +/- SD | 206 +/- 60 ng/mL | 0.657 +/- 0.191 |
| Pediatric NPC, estimated steady-state serum peak, mean +/- SD | 523 +/- 194 ng/mL | 1.667 +/- 0.618 |
| Healthy subjects, steady-state plasma peak at twice the maximum recommended dosage, geometric mean | 2,090 ng/mL | 6.661 |

[Current label, Section 12.3](https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid=5feffc0e-453d-47fa-91dd-38d4952309bc&version=5).
Conversions use arimoclomol's active-moiety molecular weight **313.78 g/mol**, not the citrate
salt's molecular weight: `micromolar = ng/mL / 313.78`.
[PubChem CID 208924](https://pubchem.ncbi.nlm.nih.gov/compound/208924).

Thus 200-400 micromolar is approximately 120-240 times the label's estimated pediatric mean
serum peak. This is our arithmetic comparison of nominal concentrations, **not** a measured
media-to-tissue or unbound exposure ratio. Serum variability, exposure duration, binding in
culture medium, intracellular accumulation and tissue relevance remain unresolved. The adult
table cannot set a child's target, and neither statistic is a dose recommendation. The
calculations and source hashes are retained in `bridge-source-audit-2026-09-19.json`.

**Decision:** arimoclomol remains a conditional probe, without preferred-treatment status.
Before confirmatory assays, a pharmacologist must justify the exposure range and measurement
plan. A signal only in deliberately high mechanistic-control conditions cannot advance as a
practical repurposing result. Cytostasis must be measured alongside segregation. The separate
[NAC exposure review](nac-exposure.md) now examines its published human
exposures and analytes. Neither drug has an established case-specific exposure window.

## NAC precedents and model boundaries

| Primary experiment | Result relevant to the proposal | Boundary |
|---|---|---|
| MRC5 human fibroblasts, 2013; 1 mM NAC before bromate-induced oxidative stress | Prevented stress-associated reduction of BubR1/Mad2 in the reported experiment | Prevention of an artificial insult, not a BUB1B-mutant rescue; p53-dependent context matters |
| Sod1-deficient mouse oocytes, 2025; 500 micromolar NAC during maturation at 21% oxygen | Improved maturation and selected chromosome-segregation outcomes | Meiosis and Sod1 deficiency differ from somatic BUB1B disease; NAC-mediated BUBR1 restoration was not demonstrated |

[2013 primary study, Figure 1](https://doi.org/10.1002/cam4.101),
[2025 primary study, Figure 7](https://doi.org/10.1002/rmb2.12622).
The latter figure reports 50 control, 39 ascorbate and 52 NAC **oocytes**, with Fisher's exact
test and Bonferroni correction. These are not counts of independent mice or patients; the
donor count for that comparison was not established in this review. Not every chromosome-count
category changed significantly. The project should control culture oxygen and distinguish
baseline-defect reversal from protection against experimental stress.

## Expanded search and screening scope

Queries, timestamps, returned record identifiers and response hashes are recorded in
`expanded-search-audit-2026-09-19.json`. Full responses remain in the public-source cache
outside Git. Search retrieval and scientific review are separate activities.

| Europe PMC query family | Indexed matches / retrieved | Review status |
|---|---:|---|
| BUB1B/BUBR1/MVA terms and arimoclomol, brand or selected code names | 0 / 0 | No match under these terms; not worldwide novelty proof |
| BUB1B/BUBR1 and N1002K, HGVS or coordinate terms | 0 / 0 | Same limitation; alternate nomenclature remains possible |
| Gene/disease terms and acetylcysteine variants or NAC | 147 / 147 | Selected primary papers reviewed; complete full-text screening not performed; NAC is ambiguous |
| Gene and broad chaperone/proteostasis terms | 871 / 100 | Retrieval incomplete |
| Gene and chaperone/proteostasis terms restricted to title/abstract | 7 / 7 | Narrowing loses relevant full-text evidence; does not replace citation tracing |
| Arimoclomol and aneuploidy/chromosome-instability terms | 4 / 4 | Three review titles plus a CDC45 primary paper; its arimoclomol passage cites other work, not its own BUBR1 rescue |
| Gene/disease terms and phenylbutyrate terms | 15 / 15 | No comprehensive full-text screening; cannot support a drug-class exclusion |

The CDC45 source was checked in full-text XML:
[primary paper](https://doi.org/10.1111/cpr.13257). Review titles served as discovery leads,
not efficacy evidence.

An expanded ClinicalTrials.gov all-fields query returned two observational studies:
[NCT07713797](https://clinicaltrials.gov/study/NCT07713797), an endometriosis biomarker study,
and [NCT00539422](https://clinicaltrials.gov/study/NCT00539422), a breast-disease polymorphism
study. Their conditions and designs do not describe an MVA drug trial. The earlier exact MVA
condition query returned zero records. These results cover one registry and specified queries;
patents, other registries, unindexed and unpublished work remain outside the completed review.

The defensible proposed contribution is a case-specific, falsifiable comparison of allele
function and downstream protection with exposure and selection controls. Neither drug is
claimed to be a new therapy, a first use, or effective for this child.
