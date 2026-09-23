# Recovered decoder aliases: exact-source role review

22 September 2026, America/Phoenix. Reader: `/root/final_coverage_audit`.

**The 32 recovered catalog–record matches are not 32 new administered-drug findings.** All 11 complete saved abstracts were read using `read_abstracts.py`, without HTML stripping, in batches of 4, 4 and 3; all three end markers were checked. The 24,929 displayed abstract characters belong to corpus SHA-256 `e5b8b017bb68af38818f1cc8c1bee47adb16afe0367981104a6b17886d2954d9`. The accompanying JSON binds each annotation and recovered alias role to its source hash and catalog ID.

There are **33 alias–record roles across 32 catalog–record matches** because mitomycin and mitomycin C identify the same matched entity in one record:

- 23 aliases: computational candidates, sensitivity predictions or secondary response/regimen associations.
- 5 aliases: actual clinical/preclinical treatments across four records.
- 2 aliases: an assay stain and an administered cell-cycle assay perturbation.
- 3 aliases: endogenous metabolites or a secretory phenotype marker.

No recovered alias is itself a genetic intervention. Some papers contain genetic experiments, whose effects must remain assigned to those experiments.

| PMID | Recovered aliases | Role and source-supported result |
|---|---|---|
| 40847407 | Eosin; L-pyroglutamic acid; pyridoxine | Eosin is a histological stain. The other two are measured liver metabolites in a genetically induced HCC model. No supplementation experiment. |
| 39512465 | Refametinib; selumetinib; trametinib | Secondary GSCAlite drug-response analysis associated high model-gene expression with resistance. No new drug administration or causal resistance experiment established. |
| 37959410 | Pembrolizumab | Administered in one adult SLE/PML case. Neurological/imaging deterioration continued and the patient died 15 months after PML diagnosis. This is nonresponse in a complex case, not proof that the drug caused harm or that all PD-1 therapies fail. |
| 32228488 | Alvocidib; nelarabine | DrugBank candidate nominations from HNSCC expression/network work; neither drug was tested. |
| 36685838 | Erlotinib; lapatinib | Secondary HCC prognostic-signature/response association. The abstract does not resolve historical measured response versus predicted sensitivity or establish a new administration experiment. |
| 19002265 | Cisplatin | Actual cisplatin-based second-line therapy in 31 NSCLC patients. Low BRCA1 expression associated with longer time to progression; this does not isolate a cisplatin effect or show restorative benefit. |
| 38143740 | Durvalumab, olaparib, paclitaxel, cyclophosphamide, doxorubicin, ixabepilone; plus fluorouracil, bosutinib, gemcitabine, gefitinib, methotrexate, mitomycin/mitomycin C, temozolomide, bortezomib, elesclomol | First group: secondary dataset associations used to propose regimens. Second group: GDSC sensitivity predictions. The source also reports functional assays of signature genes; those results do not belong to the drugs. No newly administered drug regimen is established. |
| 34771604 | Epirubicin; docetaxel | Actual six-cycle breast-cancer combination. Early serum TK1 response stratified later survival; the comparator was a biomarker-defined group, not untreated patients or either drug alone. |
| 26676187 | Nocodazole | Administered cell-cycle assay perturbation producing a substantial G2/M population. Reduced proliferation was shown with miR-210 overexpression and YES1 siRNA, not credited to nocodazole. |
| 25023465 | Substance P | Endogenous secretory tumor-subtype descriptor. No peptide administration or drug-response test. |
| 25564371 | Eplerenone | Administered mineralocorticoid-receptor blockade inhibited endothelial proliferation in a transverse-aortic-constriction heart-failure model. Dose, timing, magnitude and cardiac-function results are absent; the abstract's seven-day exposure applies to aldosterone, not eplerenone. |

## Functional-lead assessment

**No new demonstrated BUB1B, chromosome-fidelity or muscle-function rescue emerges from these 11 abstracts.** This conclusion concerns the reported endpoints, without assuming patient phase, a specific missense-protein defect, or a requirement to outperform existing nominees.

Eplerenone is a real newly recovered named-drug cellular result, and the checked decision ledger contains no eplerenone entry. It merits a narrow primary-source follow-up to determine dose/timing, model species, vascular remodeling and cardiac performance. The direction is **suppressed endothelial proliferation**; its benefit or harm cannot be determined from the abstract. Do not promote it to a beneficial repurposing candidate from this result.

The pembrolizumab case is negative clinical evidence worth preserving if the existing PD-1 route is discussed. Any primary follow-up should establish treatment timing, dose, other interventions and immune trajectory. This case must not be generalized to checkpoint inhibition in another disease or to all antibodies. The remaining recovered treatments describe cancer-treatment biomarker associations, and the rest do not contain a new treatment experiment.

## Scope and provenance limits

- These 11 readings are a separate decoder-role audit. During integration, source/ID and hash comparison added six new records to the corrected queue ledger. The other five lie outside that queue and add no cumulative reading credit here.
- No repository, original retrieval, alias policy, or historical annotation was edited. Only external audit outputs were written.
- Saved publication metadata for these 11 records does not contain a retraction label. No live retraction or publication-version check was performed.
- Current saved strings were read directly. This does not establish what the original 3,598-record annotators read or why their notes differ from the current source.
