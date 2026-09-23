# PDLSC conditioned medium: bounded source qualification

Integration state: retained source-only comparison; all 207 compound assessments remain unchanged. Independent and root reading records are tracked in the companion JSON.

**Source-level outcome: normal-cell biological-preparation comparator with a BUB1B RNA signal; no defined compound/medicine identity, BubR1 protein-supply result, or MVA rescue.** Do not give the preparation's measured benefit to individual secretome constituents or background culture reagents.

Suh HN, Ji JY, Heo JS. *Translating proteome and transcriptome dynamics of periodontal ligament stem cell-derived secretome/conditioned medium in an in vitro model of periodontitis.* BMC Oral Health 24, 390 (2024). [Publisher primary article](https://link.springer.com/article/10.1186/s12903-024-04167-z), [original PDF](https://link.springer.com/content/pdf/10.1186/s12903-024-04167-z.pdf), PMID **38539170**, PMCID **PMC10967149**. Batch-15 original review **n64**, unfiltered index **4818**, review-array index **63**. This report paraphrases/analyzes the paper and supplements under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

Original question: Which secretome components and controlled mineralization outcomes support the LPS osteogenic rescue, and does benefit extend beyond this culture model?

## Preparation, controls and functional result

Commercial human periodontal ligament stem cells (PDLSCs; CELPROGEN), passages 4–6, supplied both the preparation and recipient model. Cells at 70% confluence were washed, then cultured in serum-free α-MEM without antibiotics for 48 hours. Supernatant was centrifuged at 3,000 rpm for five minutes. The resulting **PDLSC-CM is an incompletely defined biological mixture**, not an isolated molecule, established medicine, or ordinary consumer product. No total-protein-normalized treatment dose, active-component concentration, isolated vesicle dose, or dose-response series is supplied.

Four recipient-cell conditions were compared: osteogenic medium control, 5 µg/mL *Porphyromonas gingivalis* LPS, CM, and LPS+CM. Media were replaced every two days. The osteogenic background contains 5% FBS, **50 µg/mL ascorbic acid, 1 µM dexamethasone and 3 mM β-glycerophosphate**; the CM arms received corresponding osteoinductive factors, but their final FBS content is not independently itemized. These named molecules are background culture ingredients, not individually tested rescue interventions. LPS is the experimental inflammatory challenge.

Figure 2 supplies a controlled, non-tumor culture benefit:

| Readout | Timing/statistics | Result and boundary |
|---|---|---|
| ALP activity | Figure 2A; ALP/calcium/mineralization caption n=5, mean±SD, *P<0.05 versus control and #P<0.05 versus LPS | LPS reduces the signal; CM and LPS+CM increase it. Reported ALP units are mM/100 µg protein. This is an osteogenic assay, not mitotic checkpoint function. |
| Intracellular calcium | Seven days; Figure 2B, n=5 | LPS lowers and CM restores/increases the signal. **Units conflict:** original methods p.3 says mg/100 mg protein, while native Figure 2B/original p.6 says mM/100 mg protein. Do not silently harmonize or use the chart as a calibrated absolute calcium effect. |
| OCN, OSX, RUNX2 protein | Seven days; Figure 2C, western blot and β-actin-normalized densitometry, n=3 | LPS-associated decreases are relieved by CM. These are osteogenic proteins, not BubR1. RUNX2 in LPS+CM remains visually below the control-normalized level, so “all markers fully normalized” would overstate the graph. |
| Extracellular mineralization | Alizarin Red S after 14 days; Figure 2D, n=5, stained-area analysis | LPS lowers the control-relative mineralization readout; LPS+CM increases it. Rough visual values are ~0.4 for LPS and ~1.3 for LPS+CM with control=1; CM alone ~3.2. These are graph approximations, not digitized or tabulated raw means. |

The methods specify four treatment groups but no clear separate LPS-injury interval followed by delayed CM treatment or LPS withdrawal. Discussion uses preconditioning/recovery language; a rigorously staged post-injury treatment schedule is not established. The mock medium collected from cell-free dishes for 48 hours is explicitly the **proteomic comparator**; its use as an identically aged functional-treatment control is not explicitly stated. Donor independence is not documented by the figure n values. No animal periodontal repair, clinical cohort, BUB1B-deficient genotype, chromosome-segregation assay, or general safety window is shown in this source.

## What is actually known about BUB1B/BubR1

The two omics experiments must remain separate:

1. **Secretome proteomics** identifies preparation components. Main Table 1 and Supplementary Data 1/2 include ECM/cytoskeletal proteins such as fibronectin, collagen I and keratins. No component was isolated and administered, depleted and rescued, or shown necessary/sufficient for the Figure 2 benefit. Some top-peptide-count keratins have abundance ratios below one; being listed in Table 1 does not imply each was increased. A targeted search of all supplied protein-table cells found no BUB1B/BUBR1 or corresponding O60566 accession, MAD3L/SSK1 aliases or long protein-name synonyms, verified against the [official human UniProt record](https://www.uniprot.org/uniprotkb/O60566/entry). This finite-list absence does not prove biological absence from CM.
2. **Recipient-cell RNA sequencing** supplies the BUB1B signal. Main Figure 6/Table 2 infer protein-interaction networks from differentially expressed genes. The nine BUB1B connections in Table 2 are database-network topology, not observed binding or measured BubR1 protein.

The native [Supplementary Data 3 workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs12903-024-04167-z/MediaObjects/12903_2024_4167_MOESM3_ESM.xlsx), `total (25737)` **row 7664**, identifies BUB1B (file row identifier 7662; Entrez gene ID 701). Its relevant cells are:

| Condition | Normalized log2 signal (I:L) | Raw count (M:P) |
|---|---:|---:|
| Control | 3.2658386234 | 46 |
| LPS | 2.9444137572 | 21 |
| CM | 2.7628354243 | 19 |
| LPS+CM | 4.3460645789 | 83 |

Cell **G7664** gives **2.6420372765** for LPS+CM relative to LPS. Independent calculation `2^(L7664−J7664)` reproduces this ratio. The same BUB1B value appears at **B239:C239** of `LPS VS. LPS+PDLSC-CM (776)`. The ratio is **2.1143672082** versus untreated control, whereas CM alone versus control is **0.7056363570**. Thus this is context-dependent RNA association, not a general CM-induced increase or simple return to normal baseline.

The workbook supplies one summary/count value per condition without replicate/sample identifiers, replicate columns or per-gene P/FDR columns. The primary methods describe Bowtie2 alignment, Bedtools counts, quantile normalization using EdgeR, and ExDEGA visualization. They do not map the RNA summaries to individual biological/technical samples or clearly specify the RNA harvesting time. No GSE/GSM/SRR/PRJNA accession token was found in the supplied article XML; this does not prove that no external raw dataset exists. The paper states a DEG threshold of fold change≥2 and P<0.05, but a gene-specific significance estimate cannot be independently checked from these supplied columns. BUB1B is not among the six RT-PCR-validated genes. There is no BubR1 western/quantitative protein assay, protein delivery or half-life measurement, BUB1B knockdown/rescue, or mitotic function validation. Changes in cell-cycle distribution/proliferation or bulk cell composition remain possible explanations, not mechanisms experimentally separated here.

## Eligibility and remaining limits

This paper materially answers the normal-cell culture component of n64: the whole CM preparation improves controlled osteogenic readouts in an LPS setting, and native data establish a BUB1B RNA association. It does **not** identify a defined commonly acquired compound or medicine with inherited efficacy, determine the causal constituent, establish BubR1 protein supply, or extend that benefit to whole-organism periodontal repair or MVA. Retained as a biological comparator/source annotation; no new small-molecule efficacy identity or ingredient-level benefit is assigned. The ledger and unique-reading counts remain unchanged.

## Source reviewer reading and provenance scope

- Read all scientific main text (Introduction, Methods, Results, Discussion, Conclusion) in three complete outputs with beginning/end markers, plus abstract, all six figure captions and both tables. Read all native publisher Figures 1–6 visually. Directly rendered and inspected original PDF pp.3, 6 and 11 for units, Figure 2 and BUB1B table layout.
- Retrieved all six supplements. For workbooks 1–4, read all sheet headers and structure plus targeted BUB1B rows; scanned every cell for BUB1B aliases/mitotic-checkpoint terms. This **is not a full row-by-row reading** of the 25,737-gene dataset or all protein annotations, nor raw sequencing reanalysis. The exact relevant cells and ratio were independently checked.
- Supplement 5: read the complete native single-slide text and visually inspected its four embedded raster blot images (two TIFFs converted solely for viewing). This does not claim full-slide rendering or a forensic blot-identity audit. Supplement 6: read and viewed the complete one-page workflow.
- Official Europe PMC metadata and Crossref relations were inspected; no article-level retraction/correction relation was returned in these checked records. An XML retraction link is nested within **cited reference 55**, not the present article. That cited paper/notice was not newly read and is not used as proof of a mechanism here. This is a bounded notice check, not an exhaustive publication-integrity assurance.
- Raw source URLs, hashes, failures and successful replacement figure URLs are in `retrieval-manifest.json`. The external source package contains the replay check for identity, successful asset integrity, complete main-text extraction and the targeted workbook result. Portable source URLs/hashes and review state are retained in the companion JSON; exact RNA cells and limitations are in [pdls-secretome-expression.json](pdls-secretome-expression.json). Provenance/automated scanning is not a substitute for the actual reading recorded above. The original source review made no repository edits, emails, contact requests or new reading/count claims.

## Portable descriptive replay

[replay_pdls_bub1b.py](replay_pdls_bub1b.py) uses only the Python standard library, pins the native Supplement 3 hash, and checks the exact target identity, headers, filtered row and all six descriptive ratios. Run it with `--workbook PATH_TO_NATIVE_SUPPLEMENT_3`. It neither changes the workbook nor estimates protein, causality, sample independence or differential-expression significance.

The independent check read all Methods, selected Results/Table 2, native Figure 2/original PDF page 3 and targeted supplement cells. Root separately read all Methods, the Figure 2 caption/Table 2 and selected Results/Discussion, viewed native Figure 2, inspected source headers/target cells and replayed all six ratios. Neither inherits the source reviewer’s wider reading scope. Root’s full-workbook import was intentionally stopped; the completed native-cell extraction and replay provide the stated evidence. The source workbook remains unchanged.
