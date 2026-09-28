# Ozone: adverse-direction BUB1B transcript evidence in neonatal lung

**Source qualification, 28 September 2026.** A consumer-article name hit resolved to a direct ozone exposure experiment with a BUB1B transcript result. It is evidence about a specific high-concentration environmental exposure in neonatal mice, not evidence that therapeutic ozone rescues BUB1B deficiency.

## What the source measured

Gabehart et al. exposed three-day-old BALB/c mice with their dams to **1,000 ppb ozone for three hours**. Age-matched littermates received filtered air. Whole lungs were collected six or 24 hours later; the GEO design reports four selected lungs per condition and time point. The study used the Agilent GPL7202 whole-mouse microarray. The published supplementary table lists the single annotated Bub1b probe, A_51_P490509, as **2.16-fold lower at 24 hours** (signed fold-change −2.16, P=6.90×10⁻⁵, q=0.043). Bub1b was not a significant hit at six hours.

An independent reanalysis of the GEO series matrix and GPL7202 annotation identified the same sole Bub1b probe. At 24 hours, ozone/filtered-air signal ratio was **0.462** (log2 difference −1.113; two-sided Student t-test P=6.94×10⁻⁵; Benjamini–Hochberg q=0.04347 across 41,173 probes), matching the published Bub1b row to rounding. At six hours the ratio was 0.883 (P=0.319; q=0.726). The [analysis output](ozone-bubr1b-reanalysis.json) and [reanalysis script](reanalyze_ozone_bubr1b.py) bind the downloaded GEO matrix and annotation by SHA-256.

The per-gene result is clear, but global reproduction is not exact: applying the article's P<0.05 and 1.5-fold rule gives 620 versus the paper's 621 differentially expressed features at six hours, and 866 at 24 hours in both analyses. A simple BH-plus-1.5-fold screen gives 60 features at 24 hours, whereas the article reports 54 after its FDR procedure. The BUB1B P and q values nevertheless match the published supplemental row; this note does not claim an exact recreation of the authors' complete FDR pipeline.

## Candidate decision

**Ozone does not compete with the current Track 2 nominees.** The finding is an adverse-direction RNA observation in an acute neonatal whole-lung exposure model. The source also describes broad suppression of cell-cycle programs and reduced lung-cell proliferation. It did not measure BUBR1 protein, checkpoint function, chromosome segregation, muscle output, MVA, or either nominated variant. Whole-lung composition and the developmental context prevent a cell-intrinsic interpretation.

The main catalogue contains OZONE (CHEMBL2447938) as a clinical-development record. The live ChEMBL record lists maximum phase 3 but `therapeutic_flag=false` and no first-approval date; this metadata does not establish a particular effective regimen or an MVA indication. Inhaled environmental ozone in this mouse study must not be treated as equivalent to topical, injected, or device-generated medical ozone. CPDat lists this DTXSID under an Article-kind “Construction and building materials” category; that record does not establish common purchase or personal exposure.

Retain this as a source-specific adverse-exposure comparator. It does not establish ordinary ambient-exposure effects, human harm, or the efficacy or safety of a distinct medical-ozone protocol.

## Sources and provenance

- Gabehart et al., [full primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC3930365/) and [PubMed record](https://pubmed.ncbi.nlm.nih.gov/24336422/).
- [GEO series GSE45166](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45166), with log2-transformed, quantile-normalized values, and [Agilent platform GPL7202](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL7202).
- [ChEMBL OZONE record](https://www.ebi.ac.uk/chembl/explore/compound/CHEMBL2447938) and [EPA CPDat chemical record](https://comptox.epa.gov/chemexpo/chemical/DTXSID0021098).

Retrieved series-matrix SHA-256: `678c9c54ee28aa3dd8d42600b8eafdbf649238599fd9ebb8920aad8bf1e3d109`. Retrieved GPL7202 annotation SHA-256: `bd2f3199bea338b87d12b41f451ea899ac3e52bc5fb62d8e935514504e45c62d`. The source files are retained in the external project cache, not bundled in the repository.
