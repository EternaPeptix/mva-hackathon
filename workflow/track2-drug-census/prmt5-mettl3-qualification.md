# GSK591 and STM2457 primary qualification for the BUB1B MVA census

Assessment date: 2026-09-23. Scope: a bounded review beginning with Zhao et al., *Clinical Epigenetics* 18, 151 (2026), PMID 42121254, PMCID PMC13455430, DOI [10.1186/s13148-026-02144-7](https://doi.org/10.1186/s13148-026-02144-7). This is an evidence classification for compound research, not treatment guidance.

**Disposition: retain as a CRC mechanism and combination-effect record; do not promote GSK591, STM2457, or their combination to a BUBR1-restoration or MVA-rescue nominee on this evidence.** Full primary text, all three supplements, and the key experimental figures were accessible. The paper contains direct compound experiments, so its record is stronger than an abstract-only lead. Its only explicit BUB1B connection is a transcript co-expression entry. No BUBR1 protein or functional MVA endpoint is reported.

## What connects this paper to BUB1B

Main Table 1 reports a positive Spearman correlation between METTL3 and BUB1B expression: r = 0.57408837, p = 2.40 × 10⁻⁷. The methods describe a 34-patient SRA cohort for the co-expression analysis. This is an observational transcript association in CRC; the underlying expression data were not independently reanalysed here. BUB1B is not one of the proteins examined in the displayed intervention blots. The article and supplements do not report compound-induced BUB1B transcription, BUBR1 protein abundance or half-life, kinetochore localization, spindle-checkpoint competence, chromosome missegregation, aneuploidy rescue, or any BUB1B-mutant/MVA model.

METTL3 protein stabilization by PRMT5 is a distinct result. It cannot be relabeled as BUBR1 stabilization. The direction of either drug's effect on BUBR1 remains **UNKNOWN**. A correlation does not establish that inhibiting METTL3 or PRMT5 will raise BUBR1; CRC growth suppression is not a surrogate for MVA correction.

## Compound identity and probe qualification

| Item | Verified evidence | Limit |
|---|---|---|
| GSK591 identity | The [SGC probe record](https://www.thesgc.org/chemical-probes/gsk591) identifies GSK591 as EPZ015866 / GSK3203591. Formula C22H28N4O2; InChIKey TWKYXZSXXXKKJU-FQEVSTJZSA-N. The CRC methods give CAS 1616391-87-7 and TargetMol as supplier. | EPZ015666 is a distinct related probe and must not be merged with EPZ015866. Salt products must retain their own identity. No independent analysis of the study's actual compound lot was available. |
| GSK591 activity/selectivity | SGC reports PRMT5/MEP50 H4 methylation IC50 11 nM; cellular SmD3 symmetric dimethylation EC50 56 nM in Z-138; selectivity relative to its methyltransferase panel up to 50 µM. SGC2096 is the matched inactive control, reported inactive up to 10 µM. | This is assay/panel qualification, not proteome-wide specificity at 40 µM in CRC. The CRC study does not use SGC2096 or an orthogonal PRMT5 inhibitor. |
| STM2457 identity | CRC methods give CAS 2499663-01-1 and TargetMol. The [original Nature paper](https://www.nature.com/articles/s41586-021-03536-w) identifies STM2457 as a catalytic METTL3 inhibitor. Its accessible Extended Data Fig. 1 caption describes SAM competition, selectivity against other methyltransferases, and a 468-kinase panel at 10 µM. | The Nature publisher page exposed a subscription preview and extended-data captions, not the full main paper in this retrieval. Its full main/supplement package was not qualified here. The CRC study uses 40 µM, above that stated 10 µM screening condition. |

The SGC record and primary CRC methods are the main identity authorities. Vendor claims do not substitute for the compound-specific experiments or target engagement in the CRC study. The 40 µM CRC exposure is approximately 700 times the SGC Z-138 cellular EC50 for GSK591; this cross-assay comparison flags an attribution question, not a universal toxicity threshold.

## Evidence from the main experimental figures

| Evidence | What was actually tested and observed | Interpretation |
|---|---|---|
| Fig. 1G–H | METTL3 overexpression/depletion changes PRMT5 protein; PRMT5 overexpression raises METTL3 protein in HCT116 and HT29. | Genetic reciprocal regulation; no BUBR1 readout. |
| Fig. 2A–B | PRMT5 overexpression, with or without GSK591, does not significantly change METTL3 mRNA in the displayed comparisons. | A reported null transcriptional result. It is not evidence that GSK591 restores METTL3 or BUBR1. GSK591 concentration/time for this panel are not explicitly specified in the main caption. |
| Fig. 2C–G | Endogenous co-IP; PRMT5 WT versus catalytic E444Q mutant and PRMT5 knockdown; METTL3 sDMA and cycloheximide chase. WT PRMT5 supports greater METTL3 sDMA and slower METTL3 protein loss than the catalytic mutant. | Supports a PRMT5-dependent METTL3 protein-stability mechanism. The chase is a genetic catalytic comparison, not a GSK591 protein-rescue experiment. |
| Fig. 2H–L | METTL3 overexpression raises PRMT5 mRNA; adding GSK591 does not significantly reverse that increase. METTL3 knockdown lowers PRMT5 mRNA/m⁶A; WT METTL3 rescue performs better than D395A/W398A catalytic-mutant rescue; ActD decay supports PRMT5 mRNA stabilization. | Genetic METTL3→PRMT5 regulation through RNA methylation. Preserve the GSK591 null in panel H. |
| Fig. 6 | Genetic overexpression and METTL3 knockout affect mitochondrial morphology, ATP5A1, MFN1/2, FIS1 and a COX8 mitophagy reporter. CKAP4, GM130 and TMEM192 are reported unchanged. Reporter work includes HEK293; CRC work uses HCT116/HT29. | Mitochondrial/mitophagy evidence. These are not drug-combination experiments or normal-cell drug-safety controls. Reversal of morphology is partial in displayed Fig. 6G. |
| Fig. 7A–C | METTL3-overexpressing HCT116/HT29, with GSK591 40 µM; CCK-8 at 24 h and colony formation. GSK591 reduces the increased viability/colonies relative to METTL3 overexpression alone. | Direct GSK591 functional evidence in a cancer context. Displayed groups are NC, METTL3, METTL3+GSK591; there is no parental+GSK591 arm in these panels. Thus a differential-sensitivity interaction is incompletely controlled. |
| Fig. 7D–G | HCT116 xenografts: NC, METTL3 overexpression, METTL3+GSK591; n = 6 female BALB/c nude mice/group in methods. Tumor growth/weight and Ki67 decrease with GSK591 relative to METTL3-overexpressing tumors. | Direct **single-compound** animal evidence. No STM2457 or GSK591+STM2457 animal arm is shown. No MVA model. |
| Fig. 7J–L | HCT116/HT29: control, GSK591, STM2457, both. Caption specifies 40 µM of each drug; CCK-8 at 24 h and colonies after 15 days. Combination is lower than either single agent in the indicated comparisons. | Supports stronger growth/viability suppression under one combined-dose condition. Does not by itself establish pharmacologic synergy or apoptosis. |
| Fig. 7M–O | Two CRC patient-derived organoid cases with PRMT5/METTL3 blots, images, and CCK-8; n = 4 stated for viability. Combination lowers viability versus control in both. In PDO #1 it also lowers viability significantly versus each single agent; in PDO #2 both combination-versus-single-agent comparisons are marked **ns**. | Mixed pairwise evidence across only two patient models. Do not report significant superiority over monotherapy in both PDOs. No normal-colon PDO drug arm. |

The general statistical methods state at least three experiments, mean ± SD, t tests or ANOVA in GraphPad Prism 10.1.1. Fig. 7 specifies one-way ANOVA and three independent experiments, while also specifying n = 4 for PDO CCK-8; the relationship between replicate experiments and organoid wells is not adequately described. Exact raw numerical values and uncertainty for the drug effects were not supplied in a machine-readable source-data table in the retrieved package; plotted heights have not been promoted to precise measurements here.

## Dose, synergy and control audit

Supplementary Fig. S4 contains four **single-agent** CCK-8 curves:

| Model | STM2457 concentrations (µM) | GSK591 concentrations (µM) |
|---|---|---|
| HCT116 | 0, 5, 10, 20, 40, 80 | 0, 5, 10, 20, 40, 80 |
| HT29 | 0, 0.5, 1, 5, 10, 20, 40 | 0, 5, 10, 20, 40, 80, 160 |

S4 supports concentration-dependent loss of cell activity, but its caption does not state exposure duration. No fitted IC50s or two-drug matrix are supplied there. Forty µM STM2457 gives roughly 0.4 relative HT29 viability in S4 by visual reading, versus roughly 0.6 in Fig. 7J; these are estimates from plots, and S4 timing is unspecified. Consequently the blanket main-text description of 40 µM as only slightly inhibitory needs qualification.

No Bliss, Loewe, ZIP, HSA, Chou–Talalay/combination-index calculation, isobologram, or formal drug-interaction term was found in the full article or supplements. A significant one-way comparison against monotherapy tests a difference, not an excess over an explicit additivity model. **Formal synergy remains unestablished in this paper.**

Controls present include genetic negative controls, knockdown and WT/catalytic-mutant rescues for the molecular mechanism; NC/single-agent/two-agent arms for the in-vitro combination; and vehicle controls described for mice. Not found: matched inactive chemical probes, orthogonal inhibitors, CRC target-engagement measurements tying the 40 µM viability effect to sDMA/m⁶A inhibition, or resistance/target-rescue controls for the compound phenotypes. The cell drug vehicle and its final concentration are not clearly reported. The main methods do not give an exact PDO exposure duration, dosing-refresh schedule for 15-day colony experiments, or animal vehicle composition.

## Animal and normal-model boundaries

Methods report 1 × 10⁶ HCT116 cells implanted subcutaneously; GSK591 50 mg/kg by intragastric gavage every other day, beginning 14 days after inoculation, with euthanasia two weeks after treatment starts. This is recorded solely as an experimental exposure, not a dose for humans. Fig. 7E instead shows days 10–20 on its growth-curve axis; the time origin and its consistency with the stated 28-day endpoint are unresolved. Animal allocation randomization, blinding, systemic exposure, body-weight toxicity, blood counts and organ histology are not reported in the retrieved study package.

Normal/adjacent tissues and normal single-cell populations serve observational expression comparisons. HEK293 serves a genetic mitochondrial experiment. None supplies a matched normal-cell therapeutic window for GSK591+STM2457. The paper does not establish tolerability in BUB1B-deficient normal cells or developing tissues.

## Reporting discrepancies that affect confidence

These observations are reasons to keep claims narrow, not evidence of misconduct:

- Fig. 7's overall caption uses in-vivo combination language, but the methods and displayed animal groups contain GSK591 alone with METTL3 overexpression.
- Supplementary Fig. S2 shows PRMT5-alone survival log-rank p = 0.16, despite the main text describing a significant association for either marker. METTL3 alone has p = 0.006 and the four combined strata p = 0.05. Neither establishes that the joint marker statistically outperforms METTL3 alone.
- Fig. 2's caption calls the mRNA-decay drug doxorubicin; the methods, result narrative and panel-L axis identify actinomycin D. The latter is the supported reading, with the inconsistency preserved.
- Fig. 1H's caption says METTL5, while the displayed blot and narrative say METTL3. Fig. 7H includes DLD1, which is absent from the caption's cell-line list. Fig. 7I shows six plotted points with reported r = 0.6476, p = 0.0037; the calculation/replicate basis is unclear and was not independently validated.
- The PDO #2 null pairwise results and the S4 dose/timing limitation should survive any downstream summary.

## Current access classification

On 2026-09-23, the [SGC](https://www.thesgc.org/chemical-probes/gsk591) linked GSK591 suppliers and an inactive-control request route. Live [TargetMol GSK591](https://www.targetmol.com/compound/gsk591) and [STM2457](https://www.targetmol.com/compound/STM2457) product pages listed research products and explicitly limited them to research, excluding human use and individual sales. This supports current research-reagent availability, not patient medicine access, regulatory approval, clinical dosing or MVA benefit. No purchase was made. A clinical-trial/regulatory census was outside this bounded paper qualification, so no global absence-of-trials claim is made.

## Reading completion and provenance

Primary main text: Europe PMC fullTextXML successfully downloaded; methods, results, discussion, Table 1 and all figure captions read. Main Figs. 1, 2, 6 and 7 visually inspected from publisher originals. Figs. 3–5 and 8 were assessed through captions/text only, not visual validation. The primary publisher HTML was also preserved. PMC browser access itself presented a reCAPTCHA; no challenge was bypassed.

All three original supplementary DOCX files were downloaded from publisher MediaObjects. MOESM1 contains Tables S1–S2 and references; MOESM2 contains Tables S3–S4; MOESM3 contains Figs. S1–S4. All were text-extracted and rendered with the bundled document runtime. All rendered pages were inspected (3, 2, and 3 pages respectively; MOESM2's second rendered page is blank). S4 was additionally inspected as its embedded original image. No extra experiment methods or synergy supplement was discovered. Raw experimental data, uncropped blots and reanalysis of public cohorts were not obtained/performed.

SGC's full probe page was read. The original STM2457 Nature publisher preview and exposed extended-data captions were read, but not its full main text/supplements. Europe PMC fullTextXML attempts for PMC7613134 and the GSK discovery paper PMC4753547 returned HTTP 500; PMC pages presented browser checks. These auxiliary papers are **not** counted as completed full-primary qualifications. The 2026 CRC paper itself receives full primary-text/supplement qualification, bounded by the specific visual-reading limits above.

the structured qualification records URL, retrieval time, HTTP outcome, source filename, byte count and SHA-256 for retrieved source assets. the external asset manifest hashes every retained original and derived asset; it distinguishes downloaded source bytes from extraction/render outputs. The primary article declares CC BY-NC-ND 4.0; publisher originals are preserved without modification, with local text/renders used for analysis. This assessment is original analytical prose.

## Census decision

Record GSK591 and STM2457 as research probes with compound-specific CRC evidence, and the pair as a preliminary in-vitro/PDO growth-inhibition combination. Assign **no demonstrated BUBR1 restoration**, **no demonstrated MVA functional rescue**, **no established normal-cell window**, and **no established formal synergy** from this paper. Do not infer benefit from the BUB1B co-expression row or from genetic METTL3 stabilization. This record does not upgrade provisional niacin-restoration or apitegromab-compensation nominees, nor change everolimus from a conditional probe. Any MVA relevance would require direct molecular and chromosome-segregation measurements in an appropriate BUB1B-deficient model; those data are presently absent here.


Integration note: this delegated source qualification contributes zero additional outer-corpus completion credit. It is not a clean-state review of the whole submission or biological replication.

[Structured qualification](prmt5-mettl3-qualification.json).
