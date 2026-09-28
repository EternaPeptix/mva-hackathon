# EPA CPDat v4.1 consumer-formulation screen

**Dated overlay retrieved 28 September 2026.** This processes the current CPDat 4.1 bulk-composition download surfaced by [EPA ChemExpo](https://comptox.epa.gov/chemexpo/), whose page labels the source release CPDat 4.1 (May 2025). The downloaded file's name, exact bytes and SHA-256 are recorded in the [source manifest](source-manifest.json). EPA's page reports 40,459 unique DTXSIDs across ChemExpo's broader data; that is not the size of the formulation-only subset below.

## Formulation subset and lexical screen

The bulk CSV has 1,083,810,855 uncompressed bytes. Filtering to `PUC Kind=Formulation` and nonblank DTXSID yields **1,600,483 rows**, **7,572 DTXSIDs**, **40,395 normalized names** and **12 product-use categories**. The v4.1 export has no stable product-ID column; the descriptive distinct-record key used here is `Data Document ID + Product Name`. It is not equivalent to v4.0's `cpdat_product_id`.

An exact-name screen against the same frozen **7,362-record BUB1B/BUBR1 Europe PMC corpus** yields **1,079 review-class field rows**, across **685 source/ID records** and **233 DTXSIDs**. It separately counts **10,432 low-specificity** and **193 ambiguous multi-DTXSID** matches, which are not included in the review-class table. Twenty-one current decision aliases match at least one review-class name. A field row is a chemical name in a paper title or abstract, while that source record entered the corpus through a BUB1B/BUBR1/MVA full-text query; it does not prove the chemical and gene were discussed together, that the substance was administered, or that it changed BUB1B.

Compared with v4.0, this v4.1 formulation table contains **657 additional DTXSIDs** and omits **94** previously represented DTXSIDs. Among the new IDs, **11 have review-class exact-name co-mentions** after common-word and name-fragment exclusions. The [initial triage table](new-exact-hit-triage.tsv) keeps assay/reagent, disease-model, material and drug-treatment contexts separate. These are lead-discovery hints, not efficacy findings. A separate [Article-kind screen](articles/README.md) adds consumer-article and industrial article records without merging them into the Formulation count.

## Drug-direction finding and candidate decision

The most decision-relevant new consumer-formulation match is **ibuprofen**. A primary 2026 rhesus-macaque Sertoli-cell study reports lower BUB1B RNA at 100 µM ibuprofen; its supplementary Data File S4 gives an effect-column value of −1.0040105, P=2.21×10⁻¹⁰ and adjusted P=2.89×10⁻⁸. The same paper reports lower BUB1B RNA at 400 µM naproxen (effect column −0.681759, adjusted P=6.08×10⁻⁶). These high-concentration results point away from BUB1B restoration. The study did not measure BUBR1 protein or checkpoint function, chromosome fidelity, skeletal muscle or MVA; it does not show harm from ordinary use. See the [source qualification](../naproxen-ibuprofen-bubr1b-qualification.md).

CPDat v4.1 records ibuprofen (DTXSID5020732; two formulation rows) and naproxen sodium (DTXSID7045576; six rows), both under “Food and drug.” The source's salt-form record is not silently equated with the naproxen preparation used in the cell study. These source categories do not establish current sale, prevalence or human exposure.

**Candidate ranking is unchanged:** niacin remains the proximal BUBR1-restoration hypothesis, and apitegromab remains the distinct downstream muscle-compensation hypothesis. Neither has demonstrated MVA rescue. Ibuprofen and naproxen are adverse-direction mechanistic comparators in this specific high-dose Sertoli-cell experiment, not treatment candidates. Everolimus remains a conditional probe.

## Reproduction and limits

The archive itself is not bundled here. Download the ChemExpo bulk composition file from its official [Get Data page](https://comptox.epa.gov/chemexpo/get_data/) and use the recorded SHA-256 in `source-manifest.json` to confirm the dated snapshot. Then run from the repository root:

```sh
python3 workflow/track2-drug-census/build_cpdat_bub1b_screen.py \
  --cpdat-zip /path/to/ChemExpo_bulk_composition_chemicals_20260928.zip \
  --corpus /path/to/literature-gene-unrestricted.json \
  --decision-ledger workflow/track2-drug-census/decisions.json \
  --out /path/to/cpdat-v4.1-screen-output
```

The [formulation index](consumer-formulations.tsv), [category counts](product-category-summary.tsv), [exact-name table](bub1b-exact-name-co-mentions.tsv) and [screen summary](screen-summary.json) preserve the screen outputs. The source-specific MnBP appraisal is tracked separately because CPDat's generic Phthalate article alias is not an identity match to mono-n-butyl phthalate. Hash checks establish byte and count consistency, not chemical identity, efficacy, current market presence, prevalence, purchase or exposure. CPDat does not enumerate every marketed medicine or every commonly acquired molecule; source review remains incomplete.
