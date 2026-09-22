# BUB1B assay records: source interpretation and corrections

**21 September 2026.** All 24 retrieved ChEMBL BUB1B activity rows, spanning
15 molecule identifiers and seven documents, are retained in
[`bub1b-all-target-activities.tsv`](bub1b-all-target-activities.tsv). This extends
the earlier two-clinical-drug table to experimental compounds. The table preserves
the downloaded annotations; it is not a validated list of direct binders or
therapeutic effects. ChEMBL data retain their source license stated in the census.

## A source table reverses four apparent inhibitory results

Forster et al.'s JAK3 inhibitor study supplies the original **Table S2**, file
`mmc2.xls`, sheet `Tabelle1`. The assay is explicitly singlicate. The following
cells report **residual activity, percent of control**, not percent inhibition.
Compound 4 maps to CHEMBL4564280 and compound 5 to CHEMBL5188882 through the
document's compound records. [Primary paper and supplements](https://pmc.ncbi.nlm.nih.gov/articles/PMC5119931/),
[compound mapping](https://www.ebi.ac.uk/chembl/api/data/compound_record.json?document_chembl_id=CHEMBL5653628&limit=1000).

| Original cell | Compound | Concentration | Residual activity (%) | ChEMBL activity containing the rounded value | Discrepancy |
|---|---|---:|---:|---|---|
| F53 | 4 | 100 nM | 92.3963963963964 | 26230300 | Labeled inhibition |
| G53 | 4 | 500 nM | 83.31531531531532 | 26230299 | Labeled inhibition; assigned to compound 5 at 100 nM |
| H53 | 5 | 100 nM | 109.98494731560461 | 26231114 | Labeled inhibition; assigned to compound 4 at 500 nM |
| I53 | 5 | 500 nM | 115.85549422980432 | 26231115 | Labeled inhibition |

The source heading is A5; compound headings are F10 and H10; concentrations
are F11:I11. Thus all four annotations reverse the endpoint meaning and two
also transpose compound/concentration assignment. Values above 100% remain raw
singlicate measurements: they do not establish activation. The correction is
recorded separately in [`target-assay-source-checks.json`](target-assay-source-checks.json);
the downloaded API rows are unchanged. These data do not support describing
either JAK3 compound as a potent BUBR1 inhibitor or restorative lead.

## TG101209 merits a bounded stability experiment

The 2025 renal-cancer study reports opposite thermal responses: its derivative
8h destabilizes BUB1B, whereas TG101209 stabilizes it. This is a real reported
signal worth retaining. The abstract's experimental objective and demonstrated
functional outcome concern cancer-cell killing, not healthy progenitor rescue.
[Primary abstract](https://pubmed.ncbi.nlm.nih.gov/39818011/).

ChEMBL annotates TG101209's CETSA result as a **2.19-fold signal**, interpreted
there as increased heat resistance, following
20 µM for one hour and 45°C for three minutes in Caki-1 cells, and its biochemical
IC50 as **67.3 µM**. BAY1816032 has a 1.33-fold CETSA annotation and IC50 >100 µM.
These quantitative values were checked against the API, **not independently
against the original figure or its unheated controls**. The complete article
was not obtained. Thermal stabilization does not by itself demonstrate increased
steady-state abundance, slower cellular turnover, selective direct binding,
correct chromosome segregation or usable mutant protein.

**Decision:** retain TG101209 as an experimental protein-stability probe. Compare
unheated abundance, turnover and orthogonal engagement with chromosome fidelity
and productive output before promotion. A validated exposure/eligible-product
bridge remains absent. It does not yet replace niacin's restoration nomination;
neither agent has established rescue of the actual alleles. TG101209 is not
TG101348/fedratinib. BAY1816032's established BUB1 pharmacology remains distinct
from a later BUB1B thermal observation; neither name alone proves selectivity.

The [deeper qualification](lead-qualification.md#tg101209-competing-pharmacology-matters)
adds the original compound's JAK2/FLT3/RET pharmacology and failed full-manuscript
retrieval. It strengthens the need for competing-target controls without
inventing a clinical toxicity threshold or calling an inaccessible paper negative.

## Interpretation of all seven source groups

| ChEMBL document / source | Rows | What is recorded | Current interpretation and verification depth |
|---|---:|---|---|
| CHEMBL4251633 / [PMID 29306206](https://pubmed.ncbi.nlm.nih.gov/29306206/) | 2 | BubR1 hyperphosphorylation after combretastatin-related compounds in HeLa | Cellular mitotic perturbation; no restorative result in these API annotations. Original full assay figures not verified. |
| CHEMBL4806967 / [WO2020083926A1](https://patents.google.com/patent/WO2020083926A1/en) | 3 | Residual activity >75% at 1 µM, including a dasatinib source form | Limited inhibition, not rescue. Relation and units preserved from curation; experimental source validation incomplete. |
| CHEMBL5465560 / [EUbOPEN selectivity dataset](https://doi.org/10.6019/CHEMBL5465560) | 2 | CHEMBL4082278: 95% of control at 1 µM and 80% at 5 µM | Concentrations are in `activity_properties` under `Selectivity Remarks`, not the short assay description. Curation links [PMID 27977190](https://pubmed.ncbi.nlm.nih.gov/27977190/); the original assay table was not independently verified. No restorative outcome. |
| CHEMBL5653628 / [PMID 27840070](https://pmc.ncbi.nlm.nih.gov/articles/PMC5119931/) | 4 | Apparent inhibition values for JAK3 compounds 4 and 5 | Original supplemental workbook verified: all are residual activity; two assignments also wrong, as above. |
| CHEMBL5661966 / [PMID 12719470](https://pubmed.ncbi.nlm.nih.gov/12719470/) | 3 | BubR1 phosphorylation decreases with ZM447439 | Primary abstract reports Aurora-B inhibition with failed alignment/segregation/cytokinesis and compromised checkpoint. This is not proof of direct BUBR1 inhibition. Other BUBR1-perturbed attachment-rescue settings are appraised separately. |
| CHEMBL6121040 / [PMID 39818011](https://pubmed.ncbi.nlm.nih.gov/39818011/) | 9 | Five biochemical IC50s and four CETSA values | Primary abstract supports TG101209 versus 8h thermal directions; exact figure-level validation remains open. Admit the bounded probe above, not therapeutic rescue. |
| CHEMBL6134771 / [PMID 39493221](https://pubmed.ncbi.nlm.nih.gov/39493221/) | 1 | CHEMBL6163966: 4.38% inhibition at 1 µM | Weak annotated inhibition, no restorative result. Original full assay table not verified. |

These are document groups, not seven independent demonstrations of useful BUB1B
modulation. The source errors justify checking direction, concentration and
compound identity before ranking any target-database hit.
