# Quercetin aglycone: bounded BUB1B transcript screen

**Research appraisal, 28 September 2026.** This tests whether publicly deposited quercetin-expression data support a BUB1B restoration hypothesis. It is a source-specific transcript check, not an efficacy or treatment claim.

## Finding

The primary study treated HepG2 hepatocellular-carcinoma cells with 3 μM non-conjugated quercetin after serum starvation and profiled four timepoints. The authors reported broad transcriptional changes, G1 arrest and reduced cell proliferation; they did not report a BUB1B-specific rescue. The deposited GEO series has eight processed arrays: one control and one quercetin array at 2, 4, 12 and 24 hours. GEO says the raw data were RMA-normalized with GeneSpring GX v9.0.

The official GPL570 annotation maps one probe, **203755_at**, to **BUB1B** (Entrez 701). Re-reading that row in GEO's processed series matrix gives small, mixed-direction changes:

| Time | Control log2 signal | Quercetin log2 signal | Difference | Quercetin/control |
|---:|---:|---:|---:|---:|
| 2 h | 11.0038 | 10.8706 | −0.1332 | 0.912 |
| 4 h | 10.6877 | 10.8944 | +0.2067 | 1.154 |
| 12 h | 10.9813 | 11.1916 | +0.2103 | 1.157 |
| 24 h | 11.2810 | 11.2262 | −0.0548 | 0.963 |

These are descriptive differences from one array per condition and timepoint. The deposited matrix has no replicate-level variance, so this analysis performs no statistical test and cannot establish either an increase or a null effect. It does not show a consistent BUB1B-transcript increase across time.

## Competitive interpretation

This evidence does not raise quercetin aglycone above the current Track 2 nominees. It supplies no replicated BUB1B-transcript increase, BUBR1 protein or checkpoint measurement, chromosome-fidelity result, muscle output, MVA model or variant-specific test. The cancer-cell G1-arrest result cannot be treated as beneficial checkpoint rescue in affected muscle progenitors.

Keep this aglycone appraisal separate from the exact enzymatically manufactured quercetin-glycoside training preparation and from dasatinib-plus-quercetin studies. None of those interventions inherits this transcript profile, and this transcript profile does not inherit their strength or muscle outcomes. Quercetin remains an unproven, lower-priority research comparator, not a BUBR1-restoration lead. The already identified [PMID 15309432 and PMID 21632981 expression pointers](quercetin-glycosides-qualification.md#fidelity-search-and-fair-competitive-judgment) remain unresolved; GSE15162 is a separate study and does not resolve their directions.

## Reproduction and limits

The [reanalysis script](reanalyze_quercetin_gse15162.py), [machine-readable output](quercetin-bubr1b-gse15162-reanalysis.json), and [source manifest](quercetin-bubr1b-source-manifest.json) bind the processed-data matrix and the BUB1B probe mapping. Reproduce with the official GEO matrix:

```sh
curl -L 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE15nnn/GSE15162/matrix/GSE15162_series_matrix.txt.gz' \
  -o /tmp/GSE15162_series_matrix.txt.gz
python3 workflow/track2-drug-census/reanalyze_quercetin_gse15162.py \
  --matrix /tmp/GSE15162_series_matrix.txt.gz \
  --out /tmp/quercetin-bubr1b-reanalysis.json
```

The script checks the matrix SHA-256 and refuses changed bytes, sample order or probe mapping assumptions. The full GPL570 annotation response is not bundled; its source URL, digest and exact BUB1B annotation row are preserved in the manifest. The series-matrix values are already processed; no raw CEL files or original differential-expression model were reanalyzed. This is a single cell line and one array per condition/timepoint, not a replicated BUB1B assay. It measures neither protein nor function and cannot predict a response in normal muscle or MVA.

**Primary sources:** [Notas et al., PMID 21782406](https://pubmed.ncbi.nlm.nih.gov/21782406/); [GEO GSE15162](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15162); [Affymetrix GPL570 annotation](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570&targ=self&form=text&view=full).
