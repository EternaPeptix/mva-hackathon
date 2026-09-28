# Consumer supplement label prevalence: a separate prioritization layer

**NIH DSLD export downloaded 28 September 2026; export version 9.4.0, dated 10 January 2026.** This layer counts how often exact ingredient strings appear in the structured Supplement Facts of U.S. dietary-supplement labels. It adds evidence about product-label breadth for commonly acquired nonprescription compounds; it does not measure sales, purchases, actual consumer exposure, therapeutic efficacy or BUB1B rescue.

The export readme reports **214,780 labels**: **121,959 marked “On Market”** and **92,821 marked “Off Market.”** The 2,020,128 structured ingredient rows contain 76,422 distinct exact spellings; 59,911 spellings occur in at least one on-market label. Ten on-market labels have no Supplement Facts rows in the export. Counts are deduplicated by DSLD label ID and literal ingredient string.

| Current route or consumer comparator | Exact ingredient spelling | On-market labels | Separate spelling(s) retained |
|---|---|---:|---|
| Niacin restoration hypothesis | Niacin | 10,617 | Nicotinic Acid 23; Niacinamide 627. The broad “Niacin” label does not identify the chemical form in this export. |
| NMN molecular challenger | Beta-Nicotinamide Mononucleotide | 45 | Nicotinamide Mononucleotide 23; NMN 9. These are not added together because label IDs may overlap and identity is not inferred from a string. |
| Vitamin K2 molecular follow-up | Vitamin K2 | 1,306 | Menaquinone-7 49; Menaquinone-4 19. The tested experimental subtype is not resolved by the source review. |
| Muscle-function comparator | Creatine Monohydrate | 1,180 | Other creatine salts, blends and proprietary forms remain separate. |
| Muscle-physiology reserve | Taurine | 2,644 | L-Taurine 809; retained separately by source spelling. |
| Consumer compound with bounded BUB1B RNA evidence | Quercetin | 1,788 | Quercetin Dihydrate 403; formulations are separate. |
| Chromosome-fidelity caution | Caffeine Anhydrous | 2,347 | Caffeine 1,816; caffeine-containing mixtures and salts are not pooled. |
| Consumer-compound fidelity caution | Fisetin | 52 | Source-specific product spellings are kept separate. |

These are label counts, not patient counts or therapeutic-market shares. The FDA-approved protein biologic apitegromab is outside the DSLD dietary-supplement scope, so it receives no denominator or “zero prevalence” interpretation here.

## What this changes in candidate prioritization

It supports the case that **commonness and mechanism fit are separate axes**. “Niacin” is very common as a supplement-label term, but the export identifies only 23 labels with the literal “Nicotinic Acid” string, versus 10,617 with the broader “Niacin” string. That does not mean only 23 products contain nicotinic acid; the structured export may omit the form. It does mean we cannot use the 10,617 figure as evidence that those labels match the higher-dose nicotinic-acid intervention studied in adults or its exposure.

NMN has a more direct experimental BUBR1-abundance link than niacin but substantially fewer exact on-market label entries in this snapshot. Vitamin K2 has a bounded direct Bub1b transcript/tissue-signal precedent and a visible supplement-label footprint, making it a **high-priority source-form and assay follow-up**, not a replacement nominee: the existing review still lacks verified subtype/exposure, per-cell BUBR1 function and muscle benefit. Creatine, taurine, caffeine, quercetin and other frequent entries have their own qualified functional or fidelity findings; prevalence does not transfer those results to MVA.

The elected Track 2 pair remains **niacin for a BUBR1-restoration test and apitegromab for a distinct downstream muscle-compensation test**, performed separately. The [current route comparison](../decision-report.md) remains the authority for efficacy and nomination; this label layer does not change its evidence grades.

## Reproduction and data handling

The [pinned source manifest](source-manifest.json) records the DSLD CSV ZIP URL, embedded readme version/date, HTTP metadata, byte size and SHA-256. The public asset URL is served from the `staging-s3` path; the embedded readme identifies it as NIH ODS DSLD version 9.4.0, January 2026. The separately observed production API was version 9.5.0 (May 2026), so we do not claim these snapshots are identical. The ZIP itself is not stored in this repository. Download it and replay with:

```sh
curl -fL -o /tmp/dsld-full-database-2026-01-10-CSV.zip \
  https://api.ods.od.nih.gov/staging-s3/data/DSLD-full-database-CSV.zip
python3 workflow/track2-drug-census/dsld-label-prevalence/build_label_prevalence.py \
  --zip /tmp/dsld-full-database-2026-01-10-CSV.zip \
  --out-dir workflow/track2-drug-census/dsld-label-prevalence
```

The builder verifies the input hash and the source's label, market-status and ingredient-row totals. It emits [all exact ingredient spellings](ingredient-label-prevalence.tsv), the selected literal-name comparison [candidate table](candidate-prevalence.tsv), and [machine-readable totals](summary.json). Ingredient amounts are not normalized or summed. Ingredient synonyms, salts, extracts, brand blends, and active molecules are not silently merged. Proprietary blends may hide constituents, and text outside Supplement Facts is not counted.

The NIH describes DSLD as label information from supplements sold in the United States and includes on- and off-market products; its inclusion is not an endorsement or guarantee. The source's on-market flag is therefore a database field, not verification of present retail stock. A full export is a snapshot, not a live product census. See the [NIH ODS database overview](https://ods.od.nih.gov/Research/Dietary_Supplement_Label_Database/) and [DSLD API guide](https://dsld.od.nih.gov/api-guide).
