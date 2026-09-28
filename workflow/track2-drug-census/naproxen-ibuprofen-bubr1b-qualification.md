# Naproxen and ibuprofen: adverse-direction BUB1B transcript comparators

**Source qualification, 28 September 2026.** A 2026 primary study adds a direct, but tissue- and exposure-limited, BUB1B RNA signal for two commonly acquired NSAIDs. It provides no positive treatment evidence for MVA.

## What the source measured

Crow et al. treated primary rhesus-macaque Sertoli cells with naproxen (4, 40, 400 µM) or ibuprofen (1, 10, 100 µM), using matched vehicles. RNA-seq used three biological replicates per condition after a five-day exposure protocol. The authors say the highest concentrations exceed typical human plasma concentrations; the lower-dose expression profiles were broadly similar to their vehicle controls. The study concerns Sertoli-cell and blood-testis-barrier biology, not muscle or MVA.

The deposited supplementary differential-expression tables report:

| Compound and condition | BUB1B result | Interpretation |
|---|---|---|
| Naproxen, 400 µM | Effect column −0.681759; P=5.26×10⁻⁸; adjusted P=6.08×10⁻⁶ | BUB1B RNA is lower in this high-concentration Sertoli-cell comparison. |
| Ibuprofen, 100 µM | Effect column −1.0040105; P=2.21×10⁻¹⁰; adjusted P=2.89×10⁻⁸ | BUB1B RNA is lower in this high-concentration Sertoli-cell comparison. |

The spreadsheet does not label the effect column in its visible header, so the values are transcribed without relabeling them as fold changes. The source’s DEG definition requires an adjusted P value below 0.01 and a 1.5-fold change threshold. The main text says 400 µM naproxen produced 142 differentially expressed genes, whereas its figure legend reports 158 (53 up and 105 down); that count discrepancy does not change the BUB1B row in Data File S1. The 100 µM ibuprofen comparison reports 1,567 DEGs.

The article also reports changes in transepithelial electrical resistance, but increased resistance is not itself proof of a healthier or more functional barrier. No BUBR1 protein abundance, checkpoint function, chromosome-segregation fidelity, muscle output, or BUB1B-MVA model was measured. No experiment tested Leu737Ter, Asn1002Lys, or a rescue hypothesis.

## Candidate decision

**Neither NSAID is competitive as a BUB1B-restoration or MVA muscle-benefit candidate on these data.** The direction is lower BUB1B transcript at the highest tested concentrations, opposite the proposed restoration goal. This is a mechanistic/adverse-direction comparator, not evidence that ordinary prescribed or over-the-counter use harms a person with MVA. The high-dose, nonhuman-primate Sertoli-cell result does not establish an effect at usual human exposure, in muscle, on BUBR1 protein or on chromosome fidelity. These findings do not justify changing medical use.

EPA CPDat v4.1's formulation table contains ibuprofen (DTXSID5020732; two composition rows) and naproxen sodium (DTXSID7045576; six rows). Both are labeled “Food and drug.” For v4.1, the bulk export lacks a stable product ID; those row counts are not product prevalence or exposure measures. Naproxen sodium is retained as a salt-form catalog match and is not silently treated as identical to the free naproxen used in the cell study.

## Sources and provenance

- Crow et al., [full primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC13072845/) ([PubMed record](https://pubmed.ncbi.nlm.nih.gov/41977223/)); open-access license stated in the article.
- [Europe PMC supplementary files](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13072845/supplementaryFiles), including the workbook containing Data Files S1 and S4; RNA-seq data are deposited as [SRA BioProject PRJNA1241271](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1241271).
- [EPA ChemExpo CPDat release information and bulk downloads](https://comptox.epa.gov/chemexpo/).

The local review cache records SHA-256 `715529790cff3c2e2bcfc83e16febf3142e099bdeb3ab8e77251718e92bc4a31` for the Europe PMC article XML, `41a21aea89a7e8dbe0bf26a32c4a571574e2b2c5b8c1a4547d60f76d14d261d6` for the supplementary archive, and `4dfb29bdd95e7e2c3ada041b7d931b3ea2bb1485f799d7b23c6e9710044fcdee` for its nested workbook. These hashes identify the retrieved source bytes; they do not independently validate the authors’ analysis.
