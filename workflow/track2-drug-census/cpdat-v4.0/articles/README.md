# CPDat v4.0 article-kind discovery screen

This folder holds a dated screen of composition records whose Product Use Category kind is **Article**, separate from the [v4.0 formulation screen](../README.md). “Article” includes consumer goods and industrial/article categories; it does not mean a pharmaceutical article. The source rows are public-document composition records, not evidence of current sales, use, dose or human exposure.

After requiring a DTXSID, the v4.0 table contains **63,457 composition rows**, **2,468 DTXSIDs**, **9,506 normalized names** and **14 product categories**. The fixed 7,362-record BUB1B/BUBR1 title-and-abstract screen retained **190 review-class rows** across **125 source/ID records** and **79 DTXSIDs**; low-specificity and multi-DTXSID matches are counted separately in `screen-summary.json`.

The companion [v4.1 Article-kind screen](../../cpdat-v4.1/articles/README.md) is based on the current ChemExpo bulk download. Categories and product-key schemas differ by release, so row counts are not treated as comparable purchase or prevalence estimates. All exact-name matches need context review; they are not efficacy results.

The ozone exact-name hit in both releases has a source-specific BUB1B RNA decrease in neonatal mouse lung after acute high-concentration exposure. The result and independent GEO reanalysis are documented in the [ozone qualification](../../ozone-bubr1b-qualification.md); this is an adverse-exposure comparator, not a BUB1B-restoration candidate.

Reproduce with the shared screen builder and `--puc-kind Article`. The [source manifest](source-manifest.json) binds the archive snapshot, corpus, decision ledger, code and outputs.
