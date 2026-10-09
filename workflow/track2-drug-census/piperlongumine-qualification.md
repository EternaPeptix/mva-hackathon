# Piperlongumine: cellular anticancer activity, unproven selective CDK1 mechanism

**Piperlongumine has measured cellular anticancer activity in this paper: reduced MTT signal, fewer colonies after a drug exposure followed by withdrawal, and increased Annexin-V/PI-positive cells.** The evidence exceeds a docking prediction, but does not establish direct or selective CDK1 inhibition. There is no treated-animal experiment, matched normal-cell comparison, BUBR1 response, chromosome-fidelity rescue, or useful-normal-function outcome. Retain an experimental cancer comparator, not an MVA lead or a pepper-product recommendation.

Source: Zou et al., *Frontiers in Oncology* (2021), [PMID34745968](https://pubmed.ncbi.nlm.nih.gov/34745968/), [PMC8568458](https://pmc.ncbi.nlm.nih.gov/articles/PMC8568458/), [DOI10.3389/fonc.2021.745590](https://doi.org/10.3389/fonc.2021.745590). This qualification concerns this paper only. Earlier abstract annotations for PMID23764397 and PMID28190297 are separate, incompletely qualified sources; they are not evidence of animals, normal-cell selectivity or BUB1B response in the present experiment.

## Compound, models and exposure

The administered material was piperlongumine (PL), purchased from NeOnc Technologies. No catalog, purity, batch analysis or detailed formulation is given. A 10 mM stock in DMSO was stored at −20°C; the authors state final DMSO never exceeded 0.5%. Piperlongumine is not piperine. No black-pepper, long-pepper extract, dietary preparation or oral exposure was tested, and the paper does not establish product or exposure equivalence. The description of a CMAP collection as FDA drugs does not establish approval of PL.

The three named models are SK-OV-3, CA-OV-3 and HO-8910. All are cancer-cell models, with no normal ovarian/fallopian epithelial or other healthy-cell drug comparator. The current official [Cellosaurus HO-8910 record](https://www.cellosaurus.org/CVCL_6868) identifies it as a contaminated HeLa derivative, registered ICLAC-00606. Preserve the measured response of the named cells, but do not count HO-8910 as an authenticated independent ovarian model. The paper reports a cell-bank source but no STR authentication or mycoplasma test.

| Experiment | Actual exposure and outcome |
|---|---|
| MTT | 24 or 48 h PL. Native Figure6 spans 0–100 µM for SK-OV-3/HO-8910 and 0–20 µM for CA-OV-3. Reported IC50 values at 24/48 h: SK-OV-3 49.32/16.28 µM; CA-OV-3 18.76/11.58 µM; HO-8910 12.70/6.80 µM. These are author-reported fitted values, not independently refitted. |
| Colony formation | Figure6D/E compares 0, 10 and 20 µM in each line, with significantly fewer colonies at both doses. Methods describe 200 cells/well, 24 h treatment, withdrawal, then 14 days normal culture, counting colonies >50 cells. Caption instead describes plating after treatment and 7 days growth. The durable colony suppression is positive, but exact sequence/duration is unresolved. |
| Protein abundance | Figure6F: 10/20 µM for 24 h; representative CDK1, CCNB1 and survivin bands decrease, with β-actin loading control. No plotted densitometry, dose-specific protein P-values or raw full blots are supplied in the retrieved package. |
| Apoptosis | Annexin-V/FITC and PI flow cytometry is a distinct death-associated assay. Figure6G shows increasing apoptotic-cell percentages in SK-OV-3 and CA-OV-3 at 10/20 µM, and HO-8910 at 5/10 µM; significance marks accompany both active doses. Its axes do not distinguish 24 from 48 h. Prose separately reports SKOV3 37.6% and 53.4% at 20 µM after 24/48 h; those two time-specific numbers are not separately verified from native plotted panels. |

All Figure6 data are labeled mean±SD, n=3; methods specify at least three independent experiments and explicitly three for colonies/apoptosis. Technical wells, culture passages and cytometry event counts are not stated. The statistical methods describe Student's t-test for two groups and one-way ANOVA for multiple groups; Figure6 names a multiple-comparison test without specifying which procedure. Exact P-values and individual data are not supplied. These are cell-culture experiments, not human biological replicates.

MTT alone is a metabolic readout, but the colony and apoptosis assays support cellular anticancer effects independently of MTT. No cell-free optical/reduction-interference control is described; this omission does not establish an artifact. The unexposed/0-µM comparator is present, but matched DMSO concentrations are not explicit. A direct dilution of the printed 10 mM DMSO stock to 100 µM would yield 1% DMSO, exceeding the stated 0.5% maximum; no alternate stock or solvent-exchange step is described. This is a conditional arithmetic/reporting conflict for high-dose MTT exposures, not proof that vehicle caused the 10–20 µM findings.

## What the target and BUB1B findings establish

CDK1 was selected from tumor-versus-normal expression/network analysis. PL docking scored −8.121 kcal/mol versus −8.24 for AZD5438. Those are calculated docking scores, not measured binding energies or inhibitory potencies. AZD5438 is a docking reference, not an administered positive-control arm. Figure5 visualizes a proposed pocket; there is no biochemical CDK1 kinase assay, binding measurement, cellular occupancy/thermal-shift assay, selectivity panel, CDK1 genetic phenocopy, resistant-target rescue or overexpression rescue. Lower CDK1/CCNB1 abundance can accompany altered proliferation or dying-cell composition and cannot alone assign the cause of death to direct CDK1 inhibition. No DNA-content cell-cycle analysis or direct G2/M-arrest endpoint is reported here.

BUB1B is one of ten computational hub genes and appears among the **115 tumor-upregulated genes** in SupplementaryTable1; the other list contains **181 downregulated genes**. Both lists were fully read and counted. They provide no BUB1B-specific fold change, P-value or drug-response measurement. There is no PL-treated BUB1B RNA assay, BUBR1 protein/localization assay, checkpoint test, aneuploidy measurement or MVA model. Enrichment for chromosome segregation is not a completed segregation experiment. Clinical CDK1 prognosis analysis likewise is not a treated-patient benefit.

## Source inconsistencies and limits

The native evidence should not be silently repaired:

- Figure6E's SK-OV-3 control mean is above 400 colonies and CA-OV-3 control is above 200, whereas Methods/caption state 200 plated cells. The graph is labeled colony number; no aggregation or alternative denominator explains this mismatch. Retain the reported direction while questioning the absolute colony denominator.
- Methods give 14 days after drug withdrawal; Figure6 caption gives 7 days and a different plating sequence. Figure6G has dose-only bars, leaving its mapping to the two apoptosis time points incompletely reported; HO-8910 also has substantial control apoptosis.
- Figure6F control lanes are labeled “Vector,” but no transfection/vector experiment is described. That label does not create a genetic rescue arm.
- Results/Figure1/supplement support 296 common genes (115 up/181 down); Discussion instead states 247 (103/144). Methods use |log2FC|>1, while Figure1 caption prints a ≥2 threshold. These are unresolved analysis-reporting discrepancies, not reviewer-reconstructed DEGs.
- The final Discussion says only one cell line was selected, despite Methods/Results/native Figure6 showing three named lines. The registry issue further limits the ovarian-model interpretation.

No administered-animal tumor efficacy, exposure/PK, tissue distribution, systemic toxicity, therapeutic window, resistance reversal, chemotherapy combination or normal-tissue function was measured in this paper. Claims about selective tumor killing, animal effects or synergy in its Discussion refer to other publications. They were not independently qualified here. Missing normal controls do not establish toxicity, and the cancer-cell effects do not establish MVA benefit or harm.

## Decision and reading boundary

**Piperlongumine** receives a source-qualified experimental anticancer assessment with an unproven CDK1 mechanism, using only the exact compound-name alias. No existing decision object matched the compound or this PMID at the inspected ledger; prior abstract annotations and the Batch15 follow-up remain traceable. This is one compound appraisal, with no new outer-abstract credit. Do not add “piperine,” generic pepper products or an unverified chemical alias.

Against provisional niacin restoration and apitegromab compensation, PL has real drug-treated cancer-cell outcomes but no test of the relevant BUB1B/MVA defect or of compensatory normal function. Neither provisional nominee is proven in MVA; this paper supplies no basis for promoting PL above them. A decisive next experiment would require a matched normal/BUB1B-deficient functional comparison and direct target-mechanism discrimination, rather than more docking. No clinical dose or treatment recommendation follows.

Actually read: complete scientific main text and all captions in three untruncated endmarked outputs; all six native main figures; all text/cells of SupplementaryTable1 DOCX; current PubMed/EuropePMC status metadata and the complete official HO-8910 registry record. Supplementary DOCX layout was not rendered, cited studies and raw transcriptomic data were not reread or reanalysed, and no biological replication was performed. Fresh PubMed/EuropePMC metadata contain no correction/retraction relationship; a bounded official-domain DOI-notice search found no separate notice. This is a current bounded status check, not an integrity guarantee. Raw XML, images and DOCX remain external; hashes and read-evidence references are in the accompanying manifests.

Companion: [structured qualification](piperlongumine-qualification.json). Raw-source and reading-manifest references remain in the external `piperlongumine-followup-20260923/` cache.

Integration review: Root read the complete scientific main and captions, extracted supplementary gene table, official HO-8910 registry and native Figure 6. A separate reviewer checked the core methods/results, native Figure 6, all native supplementary table cells and registry. The source reviewer inspected the remaining main figures. Raw source assets and reading logs remain external.
