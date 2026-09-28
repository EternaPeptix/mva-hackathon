# EPA consumer-product chemical overlay

This is a first, dated extension of the drug census to chemicals reported in consumer-product formulations. EPA's CPDat is built from publicly available product documents and maps named chemicals to product-use categories; a row is not proof that a product is currently sold, that anyone is exposed, or at what concentration. The EPA page now identifies **CPDat 4.1 (June 2025)** as current. This replay uses the publicly downloadable **CPDat v4.0, Figshare version 5 (posted 9 April 2025)** archive because that version's bulk file was retrievable and hash-bound in this pass. The newer current release remains to be acquired and compared. [EPA CPDat](https://www.epa.gov/chemical-research/chemical-and-products-database-cpdat), [ChemExpo](https://comptox.epa.gov/chemexpo), [v4.0 archive](https://doi.org/10.23645/epacomptox.5352997).

ChemExpo reports **40,459 unique DTXSIDs** in CPDat 4.1 overall. That total spans broader records than the 7,009-DTXSID v4.0 formulation subset below and is not a directly comparable product count. The v4.1 bulk composition file was not retrieved here.

The CPDat v4.0 product-composition table contains **1,513,875 formulation rows with a DTXSID**, covering **7,009 distinct DTXSIDs** and **36,933 normalized chemical-name strings**. These counts describe one table and release, not unique purchasable molecules or current exposure. The supplied [formulation index](consumer-formulations.tsv) retains the source DTXSID, reported names, category labels, source-row count and distinct product-ID count. [Category counts](product-category-summary.tsv) keep the source's product categories separate.

An exact-name screen against the frozen 7,362-record BUB1B/BUBR1 literature corpus returned **1,150 review-class match rows** across **725 source/ID records** and **233 DTXSIDs**. A source record can contribute separate title and abstract matches. It separately records **9,684 low-specificity hits** and **172 multi-DTXSID ambiguous hits**, which are excluded from the review list. Twenty-two current decision IDs share an exact alias with a review-class mention. These are co-mentions only: the source may discuss an assay reagent, an exposure, a therapeutic agent, a control, or background literature. The screen is not an efficacy ranking. See [the output table](bub1b-exact-name-co-mentions.tsv) and [hash-bound summary](screen-summary.json).

The only new source-level appraisal in this pass is BPA. CPDat v4.0 links BPA (DTXSID7020182) to **42 formulation product IDs** in “Food and drug” and “Home maintenance” categories. Separate cell studies report reduced BUB1B under exposed conditions, but with strong model and toxicity limits. BPA therefore joins the adverse-exposure comparator set, not the treatment shortlist; [full qualification](bisphenol-a-consumer-exposure-qualification.md).

Reproduce the lexical pass from the repository root after downloading the public v4.0 archive and the saved gene-query corpus:

```sh
python3 workflow/track2-drug-census/build_cpdat_bub1b_screen.py \
  --cpdat-zip /path/to/cpdat_v4.0.zip \
  --corpus /path/to/literature-gene-unrestricted.json \
  --decision-ledger workflow/track2-drug-census/decisions.json \
  --out /path/to/cpdat-screen-output
```

The [source manifest](source-manifest.json) binds the archive, corpus, script and published outputs. A co-mention, product-category link or chemical-name hit does not establish BUB1B target engagement, improved chromosome fidelity, MVA efficacy, or risk from ordinary use.
