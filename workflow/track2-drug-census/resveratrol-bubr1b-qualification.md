# Resveratrol: public BUB1B transcript data from A549 cells

**Research appraisal, 28 September 2026.** This is a source-specific analysis of the GEO data for the primary resveratrol/A549 paper. It tests a BUB1B-transcript question only; it does not assess treatment efficacy in MVA.

## Finding

Whyte et al. treated A549 lung-carcinoma cells with 25 μM resveratrol for 48 hours and deposited four Affymetrix arrays in [GEO GSE9008](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9008): one ethanol-control array and three resveratrol arrays. GEO reports MAS5.0 processing. The official GPL570 annotation maps probe **203755_at** to **BUB1B** (Entrez 701).

The deposited MAS5.0 signal for the control was 2471.85. The three treated arrays were 2064.50, 1975.53 and 1967.72: each is lower than that single control, by **16.5%, 20.1% and 20.4%**, respectively. Their mean treated/control ratio is **0.810**. These are descriptive array ratios. Because the dataset has only one control array, it cannot estimate control variance or support a valid inferential test; we report no P value or FDR. The GEO supplement contains selected GO-based gene lists and no BUB1B/probe row, but it is filtered rather than a complete expression table; that omission is not a genome-wide negative.

## Competitive interpretation

This result does **not** support a resveratrol-to-BUBR1-restoration claim. It is consistent with a directionally lower BUB1B transcript signal in this one tumor-cell experiment, but the one-control design, 25 μM exposure and lack of BUBR1 protein, checkpoint, fidelity, muscle or MVA measurements prevent a stronger conclusion. It does not establish harm from wine, foods, supplements or ordinary use.

Keep this molecular result separate from resveratrol's randomized exercise-associated power/torque findings and from the unresolved aged-rat BUBR1-staining paper. The transcript result weakens the case for resveratrol as a restoration lead; it does not erase its separate role as a functional muscle comparator. It remains below the current Track 2 nominees for an MVA-directed proposal because neither this dataset nor the muscle studies demonstrate benefit for the nominated genotype.

## Reproduction and limits

The [reanalysis script](reanalyze_resveratrol_gse9008.py), [machine-readable output](resveratrol-bubr1b-gse9008-reanalysis.json), and [source manifest](resveratrol-bubr1b-source-manifest.json) bind the GEO series matrix, filtered gene-list supplement and probe mapping. Reproduce from the official matrix:

```sh
curl -L 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE9nnn/GSE9008/matrix/GSE9008_series_matrix.txt.gz' \
  -o /tmp/GSE9008_series_matrix.txt.gz
python3 workflow/track2-drug-census/reanalyze_resveratrol_gse9008.py \
  --matrix /tmp/GSE9008_series_matrix.txt.gz \
  --out /tmp/resveratrol-bubr1b-reanalysis.json
```

The script checks the matrix SHA-256, series, publication link, platform and sample order. The full GPL570 annotation response is not bundled; its source URL, digest and exact probe mapping are preserved in the manifest. This reanalysis reads GEO's processed MAS5.0 matrix, not raw CEL files, and does not recreate the primary paper's full differential-expression analysis. Its three treatment arrays share a single control; differences are not statistical treatment estimates. No result transfers automatically to normal muscle, retail products, or the other resveratrol experiments.

**Primary sources:** [Whyte et al., PMID 18089832](https://pubmed.ncbi.nlm.nih.gov/18089832/); [GEO GSE9008](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9008); [GPL570 annotation](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570&targ=self&form=text&view=full).
