# FDA approval refresh, 4 October 2026

The [current official CDER novel-approval list](https://www.fda.gov/drugs/novel-drug-approvals-fda/novel-drug-approvals-2026) contains **45 numbered products**, five beyond the saved 40-product layer. All prior table fields remain unchanged. The [dated refresh](approval-refresh-20261004.json) supplements the immutable [historical annotations](approval-overrides.json); it does not overwrite them or enlarge the **15,357-entity** catalog denominator. [Structured source review and acceptance audit](fda-approval-refresh-20261004-audit.json).

| FDA row | Product / active ingredient | Approval date | Existing catalog cross-reference |
|---|---|---|---|
| 41 | Onswik / insulin efsitora alfa-gobe | 23 September 2026 | CHEMBL5095357; suffix requires product-identity review |
| 42 | Lyrfigtu / lirafugratinib | 23 September 2026 | CHEMBL5314555; administered hydrochloride remains distinct |
| 43 | Juvmo / tavapadon | 25 September 2026 | CHEMBL3697617 |
| 44 | Atebrioz / zilurgisertib | 25 September 2026 | CHEMBL5314579; administered fumarate dihydrate remains distinct |
| 45 | Emcitate / tiratricol | 28 September 2026 | CHEMBL41632; separate EMA product overlap |

The inventory has 41 whole-product/name cross-references, two component-only combination cross-references and two unmatched products. Its 47 ingredient occurrences are not unique molecules. All five additions already have catalog records and no exact appraisal joins in the current ledger. This refresh adds **zero molecules, curated appraisals or accepted abstract screens**. Numbers/dates above derive from the official table; identity and matching checks derive from the saved catalog/alias layers bound in the structured audit.

The new products received bounded label-section review of identity, indication, mechanism, warnings and clinical studies. Whole-label or complete-primary-paper review is not claimed. FDA-hosted labels were recovered for [Lyrfigtu](https://www.accessdata.fda.gov/drugsatfda_docs/label/2026/220425Orig1s000lbl.pdf) and [Emcitate](https://www.accessdata.fda.gov/drugsatfda_docs/label/2026/220963s000lbl.pdf). The inspected Drugs@FDA application pages explicitly lacked labels for Onswik, Juvmo and Atebrioz; [Onswik](https://pi.lilly.com/us/onswik-uspi.pdf?s=pi), [Juvmo](https://www.rxabbvie.com/pdf/juvmo_pi.pdf) and [Atebrioz](https://files.mirumpharma.com/atebrioz/Atebrioz-USPI-MedGuide.pdf) sponsor-hosted PIs have explicit separate provenance.

| Product | What the inspected evidence measures | Pending scientific qualification |
|---|---|---|
| Onswik | Label glycemic endpoints. A **partial** QWINT-2 PRO report measures SF-36 physical-function domains and reports comparative HRQoL nulls. | Domain estimates, denominators, formulation and cohort linkage remain unqualified. Treatment-burden positives do not establish physical-function gains. [Primary report](https://doi.org/10.1007/s13300-025-01833-5). |
| Lyrfigtu | Tumor response/duration in a single-arm study; baseline ECOG is not treatment function. | Congress-level compound/gene cooccurrence has no verified same-experiment link; its body remains unread. [Label section 14](https://www.accessdata.fda.gov/drugsatfda_docs/label/2026/220425Orig1s000lbl.pdf). |
| Juvmo | Parkinson's motor-related daily living and diary ON/OFF durations. | Primary TEMPO papers, hierarchy, denominators and harms require appraisal. Label narrative/table P-value and interval-label discrepancies are retained. [Sponsor label section 14](https://www.rxabbvie.com/pdf/juvmo_pi.pdf). |
| Atebrioz | CT heterotopic-ossification volume; CAJIS supplies baseline eligibility/stratification. | Function change cannot be inferred from structural imaging. Original PROGRESS outcomes and unresolved development-code aliases require review. [Sponsor label section 14](https://files.mirumpharma.com/atebrioz/Atebrioz-USPI-MedGuide.pdf). |
| Emcitate | MCT8-deficiency thyroid endpoints and peripheral physiological signs. | Brain or motor rescue does not follow from T3 normalization. Original Triac studies and age-dependent neurological outcomes need qualification. [FDA label section 14](https://www.accessdata.fda.gov/drugsatfda_docs/label/2026/220963s000lbl.pdf). |

Searches retain exact queries, source IDs, dates and limits. The 79 retrieved bibliographic pointers have zero exact source/ID matches against the 1,640 accepted outcome records; they are **pending retrieval leads**, without abstract-reading credit. Congress aggregates are not same-experiment evidence. No-hit snapshots do not establish absent efficacy. The inspected sources establish no BUB1B/MVA rescue or ranking against niacin, NMN or apitegromab.

Both raw transport-byte hashes and the historical newline-normalized text hashes are preserved. The source has two CRLF sequences; silently comparing the two hash representations would produce a false mismatch. Root replayed the parser and checked complete identities, counts, raw/text hashes, unchanged historical bytes and empty new-product appraisal joins. `fda_freshness.py --output <dated-file>` now permits a refresh without replacing the historical annotations. The archive retains exact reading scopes and source-specific safety/formulation details.
