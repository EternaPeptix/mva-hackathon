> Package copy: local links/navigation are adapted for this finite review package. Labels identify underlying files that are not included. Scientific text is retained; approved status/source metadata changes are recorded separately; original and packaged hashes are recorded in PROVENANCE.

> **Retained dated method detail.** This appendix preserves the original analysis or proposal at its stated date. The current [Track 1 report](../jeweled_track1_report_submission.md), [Track 2 report](../jeweled_track2_report_submission.md) and [methods form](../jeweled_methods_form_answers.md) govern interpretation; see the [revision and evidence scope](../PROVENANCE.md#evidence-scope). Historical delivery or ranking statements are not current status.

# Acetylcysteine: exposure evidence and assay implications

19 September 2026. Public-source review for the [Track 2 proposal — archival document not included](../PROVENANCE.md#archival-material).
These are published research exposures, not a regimen for the child. No drug experiment has
been performed by this project.

## What concentration was measured?

Reduced NAC, low-molecular-weight disulfides, protein-associated forms, cysteine and glutathione
must remain separate measurements. A six-volunteer crossover study reported different oral
bioavailability for reduced NAC (4.0%) and total NAC (9.1%). Only its indexed abstract was
reviewed here; the full analytical method was not inspected.
[Olsson et al., 1988](https://pubmed.ncbi.nlm.nih.gov/3360052/).

The current ACETADOTE label describes total-acetylcysteine pharmacokinetics after intravenous
administration and reports protein binding of 66-87%. It also describes formation of cysteine
and disulfides, with cysteine contributing to glutathione metabolism. Its established use is
acetaminophen poisoning; this does not establish a chronic MVA treatment window. DailyMed's
history API identifies version 17, published 30 July 2025, as the latest returned version.
[Label, Sections 1 and 12.3](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=472f158a-5ab9-4308-8e49-1116e6ea3d39).

**Inference for the assay:** multiplying a published total concentration by a generic
"unbound fraction" does not establish the concentration of reduced NAC inside cells. Nor
does low parent-NAC exposure exclude a cysteine/glutathione-mediated effect. Both the chemical
species and the functional response need measurement.

## Human studies used for comparison

A 30-person adult study used oral 600-mg tablets, initially once and then twice daily for
three days before the final sampling dose. Its analytical method reduced intermolecular
disulfides with dithiothreitol before measurement. Table 2 gives baseline-corrected means of
individual peak concentrations; these differ from the maximum of a cohort's mean time curve.
The two reported groups each contained 15 adults. The following rows use its repeated-dose
phase, without pooling groups or inferring a child's exposure.
[Papi et al., Table 2 and Assessments](https://doi.org/10.1007/s12325-020-01542-4).

| Research population and formulation | Reported peak, mean +/- SD | Calculated micromolar |
|---|---:|---:|
| 600-mg study, Chinese cohort, n=15 | 3,513.1 +/- 1,142.6 ng/mL | 21.53 +/- 7.00 |
| 600-mg study, Caucasian cohort as reported, n=15 | 3,419.3 +/- 920.0 ng/mL | 20.95 +/- 5.64 |

The signal after disulfide reduction cannot be equated with free reduced NAC at collection.
Nominal cell-culture concentrations of 500-1,000 micromolar in the previously reviewed BUBR1-related
studies are about **23-48 times** these group mean peaks. This is a nominal comparison across
different matrices and analytes, not an intracellular exposure ratio or proof of inefficacy.

A separate adult crossover study deliberately examined an **11-g single dose** chosen to
approximate an acetaminophen-antidote loading dose. Its 29-person pharmacokinetic population
had mean total-NAC peaks of 26.5 micrograms/mL for effervescent tablets and 28.4 micrograms/mL
for oral solution, approximately 162 and 174 micromolar. The analysis measured total NAC
after reduction/derivatisation. This acute high-dose context cannot justify chronic pediatric
exposure, and its bioequivalence result is not MVA efficacy evidence.
[Greene et al., Table I and bioanalytical methods](https://doi.org/10.1016/j.curtheres.2016.06.001).

All conversions use **163.20 g/mol** for acetylcysteine:
`micromolar = ng/mL / 163.20`, or `micrograms/mL * 1000 / 163.20`.
[PubChem CID 12035](https://pubchem.ncbi.nlm.nih.gov/compound/12035).
Source values, uncertainty types, calculations and hashes are recorded in
`nac-exposure-source-audit-2026-09-19.json`. No concentration-to-dose extrapolation is made.
The DailyMed CETYLEV name query returned no current records; historical research data are
retained as such, without inferring current market availability or withdrawal from that result.

## Exposure does not establish biochemical benefit

In a separate study, 24 adults receiving hemodialysis and seven healthy controls received
sustained-release oral NAC. Renal-disease participants had markedly altered exposure.
After 14 days, homocysteine decreased, but several other endpoints, including measured total
glutathione and protein carbonyls, did not significantly change. These results do not describe
the child's renal function and cannot determine an MVA dose or response.
[Nolin et al., Methods and Pharmacodynamics](https://doi.org/10.2215/CJN.00210110).

The inference is to measure engagement and benefit separately. Drug exposure, a lower ROS
probe signal, a change in a circulating thiol, and rescued cellular function are different
outcomes. The completed literature review establishes these distinctions; it does not establish
a clinically usable exposure window for this case.

## Consequences for the proposed experiment

1. **Define the analyte and route before choosing conditions.** Record formulation, nominal and
   measured media exposure, sampling time, reduced versus total NAC, and analytical recovery.
   Do not transfer an acute intravenous or antidote-loading exposure into a chronic proposal.
2. **Test the metabolic route.** Measure relevant intracellular thiols and glutathione redox
   state with appropriate sample handling. Use a functional endpoint alongside engagement.
   A low parent-drug concentration alone must not falsely reject an indirect mechanism.
3. **Exclude assay chemistry.** Include cell-free drug/probe controls, orthogonal endpoints,
   matched media pH and culture oxygen, and a washout comparison. Distinguish an effect on the
   measurement from an effect in viable cells.
4. **Keep exposure and selection gates.** A response seen only in deliberately high laboratory
   conditions stays mechanistic evidence. A fall in chromosome-error frequency caused by arrest,
   death or selective expansion does not count as repair. Persistent segregation errors with
   better function must be described as downstream protection.

Acetylcysteine remains an exploratory comparator. The literature portion of this exposure
review is complete within the stated sources. The case-specific exposure window, feasibility
and functional benefit remain to be established through the proposed work.
