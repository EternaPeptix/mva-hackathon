# Apitegromab: requested fracture surveillance is new regulatory context; statistical gaps remain

**Source-qualified regulatory-context update, 27 September 2026.** The official BLA 761463 approval letter was newly obtained and read. It adds **requested enhanced pharmacovigilance (EPV)** for fractures and a separate required pregnancy study. It supplies no new motor-function estimate, fracture event count, causal adjudication or incidence estimate. The existing provisional compensation priority is unchanged; this packet does not propose a new appraisal ID or a clinical regimen.

## Source identity and bounded access

The [official Drugs@FDA application](https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=overview.process&ApplNo=761463) identifies ISEMBYLD, apitegromab-mstn, Scholar Rock, BLA 761463, and original approval on 11 September 2026. Its currently exposed history links the label and an approval letter; **no efficacy, statistical or multidisciplinary review dossier link is present**. The previously failed label PDF was not requested again. No journal route was reopened.

The actually linked [nine-page approval letter](https://www.accessdata.fda.gov/drugsatfda_docs/appletter/2026/761463Orig1s000ltr.pdf) has FDA Reference ID **5868054** and an electronic signature dated **11 September 2026, 04:42:12 PM**. It addresses BLA 761463 for Isembyld/apitegromab-mstn and the same SMA-on-SMN2-treatment indication as the retained label. The downloaded file has a valid PDF signature, nine pages, FDA author metadata, and September 14 PDF creation/modification metadata; these file dates are not substituted for the letter's signed approval date. Its SHA-256 is `846b3931180fc2fcf7370ccb316ae37b5b7710c264543c1a7df8f3b0a788a4b1`.

## What is new, and what it does not mean

| Letter source | New source-qualified information | Interpretation boundary |
|---|---|---|
| Pages 7–8, **REQUESTED ENHANCED PHARMACOVIGILANCE (EPV)** | FDA requests that the sponsor submit all serious and non-serious domestic and/or foreign fracture cases as **15-day Alert reports**. The first request has no separately stated end date. | Sponsor-to-FDA reporting request, not new fracture events, a clinician dosing mandate, an MVA-specific contraindication or a fracture 505(o) postmarketing requirement. It is not described as a REMS. |
| Pages 7–8, EPV | FDA requests fracture narrative summaries/analyses within required periodic safety reports: quarterly during the first three post-approval years, then annually through the fifth year after initial US approval. The requested analysis includes interval/cumulative cases and causality context, including timing, treatment duration, confounders, underlying risk, outcome and dechallenge/rechallenge. | This requests future assessment; it does not supply an adjudicated cause, effect size or risk estimate. |
| Page 8, EPV continuation | Annual cumulative utilization and exposure reporting, including doses/units administered, is requested through the fifth year after initial US approval. | This is a monitoring plan. No unique-patient/person-time denominator or new incidence is supplied. The five-year language applies here and to the narrative reports, not an invented termination of the separate 15-day request. |
| Pages 4–6, **POSTMARKETING REQUIREMENTS UNDER 505(o)** | Separate **required study 4891-1**: worldwide single-arm pregnancy safety follow-up for at least ten years; draft/final protocol dates January/September 2027, study completion March 2038 and final report March 2039. | Distinct pregnancy PMR, not a fracture PMR or new observed pregnancy-outcome evidence. |
| Page 4, **ADVISORY COMMITTEE** | The letter explains the application was not referred because it did not raise significant safety or efficacy issues in the intended population. | That regulatory statement concerns the intended SMA population. It does not erase the labeled fracture signal or establish universal, MVA, bone-development or oncology safety. |

The retained label's randomized fracture counts remain **5/53 at 10 mg/kg, 3/53 at 20 mg/kg, and 1/50 placebo**. Its open-label 20-mg/kg extension reports **22 patients with fractures**, including five serious femur-fracture events. These are **inherited label findings**, not newly discovered events. The randomized comparison and higher-dose extension remain separate, and the latter is not a labeled-dose incidence estimate. The new letter strengthens the documentation that fracture follow-up is an active regulatory concern; it does not quantify an additional benefit–risk penalty or change the nominee ordering by itself.

## Trial documents, history and unresolved statistics

The fresh [ClinicalTrials.gov API record for NCT05156320](https://clinicaltrials.gov/api/v2/studies/NCT05156320) authenticates **SAPPHIRE / SRK-015-003**, Scholar Rock, 188 participants, and the primary HFMSE outcome in the main efficacy population. It distinguishes the 2–12-year main population (10 mg/kg, 20 mg/kg, placebo) from the 13–21-year exploratory population (20 mg/kg, placebo). The record's **last update posted is 22 January 2026**, with submission 16 January; derived `versionHolder: 2026-09-25` is an API-data version marker, not the study's update date.

The current API response has **`hasResults: false`**, no `resultsSection` and no `documentSection`. Therefore it supplied no posted outcome tables, full protocol or SAP. This is a statement about the retrieved current public record, not proof that results or an SAP never existed, were never submitted, or were never prospectively specified.

The public history page returned only an application shell to the web reader. No browser surface was available; native access was unavailable while the host was locked. To identify the site's actual history-data request, the directly linked public scripts were read as text, without execution. One ordinary unauthenticated request to the exact observed public website data path returned **403**. It was not retried or bypassed. Two exact-source official-domain search queries returned mismatched trial records; none was followed or used as SAPPHIRE evidence. **No historical version was successfully read**, and no trial-amendment chronology is claimed.

| Requested question | Outcome of this follow-up |
|---|---|
| Prespecified pooled/20-mg/10-mg test ordering and familywise error control | **UNRESOLVED.** Current endpoint listing and approval letter do not supply the SAP hierarchy. Retain pooled positive, 20-mg null and nominal 10-mg result exactly as before. |
| Protocol/SAP amendments and their dates | **UNRESOLVED.** No protocol/SAP recovered; history inaccessible in this bounded pass. |
| Missingness, analysis sets and sensitivity/tipping-point results | **UNRESOLVED.** No new results tables or reviewer analysis obtained. |
| New fracture event timing, adjudication and incidence | **UNRESOLVED.** The letter requests future analyses; it supplies none. |
| Whether new source changes comparative MVA efficacy or nominee ordering | **No justified change.** Regulatory-context delta only. No MVA rescue, BUBR1 restoration or cross-trial superiority is established. |

## Reading and disposition

The entire nine-page letter text was read. Native pages **1 and 4–9** were inspected, including the complete relevant safety/EPV/PMR sections, timetable and signature. Pages 2–3 were read as extracted administrative text only. There are no native trial efficacy or fracture-result tables in the letter. Current registry identification/status/design/arms/outcomes/eligibility/reference/IPD fields were read; no historical record, full SAPPHIRE scientific main, SAP, reviewer dossier or trial-native outcome table was read in this pass. Existing abstract/label findings retain their earlier provenance.

Independent review accepted this narrow regulatory-context refinement. The [existing apitegromab appraisal](apitegromab-review.md) retains its efficacy figures, dose/age boundaries and provisional research priority. The [structured qualification](apitegromab-regulatory-qualification.json) records the corrected source metadata and separate reading scopes. This adds no candidate, observed outcome, original-trial reading or original-queue credit. Raw sources and earlier frozen packets remain external and unchanged.
