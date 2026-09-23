# Batch09: 128 complete saved-string revalidations

All **128 selected source records** were read and annotated through **32 untruncated outputs of four records**, with every end marker observed: **181,225 abstract characters**. The selection is the next 128 in the immutable original 3,598 order after excluding all 1,274 accepted records, including the explicit food 64 admission. No biological, date, publication-type or outcome filter was applied.

The completed batch gives **1,402 accepted revalidations, 3,563 outer-ledger records and 2,196 pending saved abstracts**. The original 3,598 selection objects and all 1,274 prior acceptances remain unchanged. The 143 initially empty abstract fields retain separate recovery statuses; this batch does not re-audit or erase those recoveries. This batch does not exhaust the corpus, establish 128 independent studies, or create 128 drug appraisals.

Actual benefits, nulls and adverse effects are preserved under their tested interventions:

| Source | Retained result | Qualification boundary |
|---|---|---|
| MED:41603595 | A local PEGSN-PGA/rhBMP-2/rapamycin hydrogel reportedly improves aged bone integration/regeneration; rhBMP-2 alone aggravated senescent-cell adipogenic imbalance. | Species, doses, component controls and quantified bone-function endpoints are absent from the abstract. This is not rapamycin monotherapy or everolimus equivalence. |
| MED:31310234 | Electroporated recombinant Ndc80/Mis12 complexes functionally complemented depleted proteins to support checkpoint signaling and segregation; modified Spindly restored localization. | Living cultured-cell complementation, despite “in vivo” in the title; no animal treatment or BUB1B delivery. |
| MED:31122175 | An unnamed Aurora-B inhibitor rescued misalignment induced by an unnamed CENP-E inhibitor. | Preserve actual positional rescue without inventing compound names or postmitotic fidelity. |
| MED:40661163 | NP-110 increased Tert mRNA in the stated deficient mouse model; PD98059 changed shelterin expression in fibroblasts. | Marker effects differ from genetic RHAMM isoform/telomerase results. Organ function and PD98059 target attribution need qualification. |
| MED:40180916 | Ginkgolide B reduced ovarian metastatic tumor burden in mice. | Actual antitumor benefit; dose, target specificity and normal-cell fidelity are not given. |
| MED:19554570 | Ispinesib produced tumor delays and maintained complete responses in several xenografts. | High toxicity, maximum-tolerated-dose exceedances and evaluability exclusions remain visible; no normal-tissue rescue. |
| MED:39204193 | Isolated 10-HDA suppressed lymphoma-cell survival at a reported 24 h IC50 of 496.8 µg/mL; comparator IC50s were about 1,000 µg/mL in source-labelled LO2 cells and above 1,000 in fibroblasts. | Retain the comparator difference without inferring a clinically useful window, authenticated normal tissue, consumer royal-jelly efficacy or MVA benefit. |
| MED:27558154 | Reduced-dose BI6727+MK1775 matched lethal-dose BI6727 activity with no regrowth after washout; lapatinib was cytostatic. | In vitro result; equal potency is not automatically quantified synergy. Candidate-gene alterations lacked a significant overall-survival association. |
| MED:21216274 | Mad2 haploinsufficiency protected progenitors against Ara-C and reduced E3330 responsiveness. | Genetic/slow-cycling context; no protection after cyclophosphamide/irradiation and no altered short/long-term HSC repopulation. |
| MED:27082996 | OTSSP167 abrogated the checkpoint through Aurora B/BUB1/Haspin activity not recapitulated by MELK RNAi. | Do not treat it as a selective MELK probe, substitute BUB1B for BUB1, or infer restoration. |

Other genuine cancer-context positives include bortezomib-mediated MCL apoptosis, ONC201 radiosensitization, an unspecified EZH2-inhibitor/YM155 combination, and paclitaxel/navitoclax responses under genetic p31comet depletion. Their actual agents and limits remain in each record. Reviews, biomarkers, endogenous hormones, gene/protein manipulations, assay dyes, diagnostic antibodies and plant growth regulators are separately identified. Megator “Mtor” is not the mTOR drug target; kinase/target names do not identify absent compounds. Senescence-marker reductions and survival of selected cells do not by themselves establish faithful functional recovery.

There are **14 bounded source questions**, not 14 drug nominations: MED:38589728, MED:36612203, MED:21340479, MED:40180916, MED:19554570, MED:36213002, MED:28299189, MED:41603595, MED:27558154, MED:31122175, MED:40661163, MED:40924783, MED:31310234 and MED:39204193. The survivor-risk abstract MED:40924783 contains an internally inconsistent meningioma confidence interval; its bound is not silently corrected. No full-paper chase was performed in this batch.

Saved publication status is retained exactly:

- **Retracted:** MED:31410029, linked to MED:32273722 as `Retraction in`. The reported findings are labelled as retracted, and no retraction reason or notice-reading credit is inferred.
- **Selected preprint:** PPR:PPR666725, linked as `Preprint of` MED:39789219. No explicit version list/number is supplied in the snapshot. The journal counterpart was not newly read here.
- **Errata:** MED:33475988→MED:33950420 and MED:33482819→MED:33536016. Notice contents remain unassessed.
- **Other types:** two conference abstracts, one editorial/comment, and 38 records bearing review-related publication labels. The exact `Review` tag occurs 36 times; multi-valued publication-type counts are not disjoint categories. Eleven `Preprint in`, one `Preprint of`, and two `Update of` links remain source lineage, not independent-study credit.

The root declared complete reading and annotation comparison for eight predetermined records at zero-based indices 0, 15, 32, 63, 64, 79, 96, 127. Literal output chunks were `f923fb` and `82c06c`; the valid annotation comparison was `82d415`. The earlier wrong-key lookup `521207` printed no annotations and receives no comparison credit. No substantive correction was required. This is an eight-record spot-check, not an independent review of all 128 or an error-rate estimate.

The full verifier replay passed source identity/hash equality, every literal string and end marker, original annotation reconstruction, complete metadata preservation, both root sample logs and the existing food-admission/batch08 gates. Machine checks validate artifacts and declarations; they cannot prove viewing or scientific correctness. No new primary qualification or demonstrated MVA efficacy is claimed. Niacin and apitegromab remain provisional and unproved for MVA; everolimus remains conditional.

External files: `selection.json`, `review.json`, `batch-evidence.json`, `validation.json`, original `notes-batch-01..32.json`, literal `batch-01..32.jsonl`, and the root sample files. `verify_revalidation_batch09.py` exports `check(root=...)` for compact repository integration. Full replay additionally accepts `--corpus` and `--evidence-root`, with `--selection` and `--review` for external paths. `build-batch09.py` rebuilds the external package from the unchanged source, selection and notes.
