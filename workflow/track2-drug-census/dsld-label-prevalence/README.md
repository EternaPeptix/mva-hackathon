# Consumer supplement label prevalence: a separate prioritization layer

**NIH DSLD export downloaded 28 September 2026; export version 9.4.0, dated 10 January 2026.** This layer counts how often exact ingredient strings appear in the structured Supplement Facts of U.S. dietary-supplement labels. It adds evidence about product-label breadth for commonly acquired nonprescription compounds; it does not measure sales, purchases, actual consumer exposure, therapeutic efficacy or BUB1B rescue.

The export readme reports **214,780 labels**: **121,959 marked “On Market”** and **92,821 marked “Off Market.”** The 2,020,128 structured ingredient rows contain 76,422 distinct exact spellings; 59,911 spellings occur in at least one on-market label. Ten on-market labels have no Supplement Facts rows in the export. Counts are deduplicated by DSLD label ID and literal ingredient string.

| Current route or consumer comparator | Exact ingredient spelling | On-market labels | Separate spelling(s) retained |
|---|---|---:|---|
| Niacin restoration hypothesis | Niacin | 10,617 | Nicotinic Acid 23; Niacinamide 627. The broad “Niacin” label does not identify the chemical form in this export. |
| NMN molecular challenger | Beta-Nicotinamide Mononucleotide | 45 | Nicotinamide Mononucleotide 23; NMN 9. These are not added together because label IDs may overlap and identity is not inferred from a string. |
| Vitamin K2 molecular follow-up | Vitamin K2 | 1,306 | Menaquinone-7 49; Menaquinone-4 19. The tested experimental subtype is not resolved by the source review. |
| Amla / *Phyllanthus emblica* comparator | Amla | 340 | Separate forms: Amla Fruit Extract 92; Amla extract 77; Amla fruit extract 57; Phyllanthus emblica fruit extract 7. Counts are not pooled; none identifies the study extract. |
| Muscle-function comparator | Creatine Monohydrate | 1,180 | Other creatine salts, blends and proprietary forms remain separate. |
| Muscle-physiology reserve | Taurine | 2,644 | L-Taurine 809; retained separately by source spelling. |
| Consumer compound with bounded BUB1B RNA evidence | Quercetin | 1,788 | Quercetin Dihydrate 403; formulations are separate. |
| Chromosome-fidelity caution | Caffeine Anhydrous | 2,347 | Caffeine 1,816; caffeine-containing mixtures and salts are not pooled. |
| Consumer-compound fidelity caution | Fisetin | 52 | Source-specific product spellings are kept separate. |
| NAC redox/fidelity comparator | N-Acetyl-L-Cysteine | 557 | Other exact NAC spellings: N-Acetyl Cysteine 529; N-Acetyl L-Cysteine 210; N-Acetyl-Cysteine 99; N-Acetylcysteine 60. Counts are kept separate and not summed. |
| L-ascorbic-acid fidelity comparator | L-Ascorbic Acid | 8 | Ascorbic Acid 161; broad Vitamin C 19,525. These label strings are not pooled or treated as the same tested preparation. |
| Linoleic-acid checkpoint-context comparator | Linoleic Acid | 875 | The BUBR1/SAC result belongs to the B7E2-plus-vinorelbine mixture; isolated linoleic acid was tested for cancer-cell apoptosis synergy, not BUBR1 rescue. |

These are label counts, not patient counts or therapeutic-market shares. The FDA-approved protein biologic apitegromab is outside the DSLD dietary-supplement scope, so it receives no denominator or “zero prevalence” interpretation here.

## What this changes in candidate prioritization

It supports the case that **commonness and mechanism fit are separate axes**. “Niacin” is very common as a supplement-label term, but the export identifies only 23 labels with the literal “Nicotinic Acid” string, versus 10,617 with the broader “Niacin” string. That does not mean only 23 products contain nicotinic acid; the structured export may omit the form. It does mean we cannot use the 10,617 figure as evidence that those labels match the higher-dose nicotinic-acid intervention studied in adults or its exposure.

NMN has a more direct experimental BUBR1-abundance link than niacin but substantially fewer exact on-market label entries in this snapshot. Vitamin K2 has a bounded direct Bub1b transcript/tissue-signal precedent and a visible supplement-label footprint, making it a **high-priority source-form and assay follow-up**, not a replacement nominee: the existing review still lacks verified subtype/exposure, per-cell BUBR1 function and muscle benefit. Creatine, taurine, caffeine, quercetin and other frequent entries have their own qualified functional or fidelity findings; prevalence does not transfer those results to MVA.

Amla is a notable common-consumer challenger: one study of a prepared aqueous fruit extract in nontransformed colon epithelial cells reported higher BubR1 mRNA alongside improved spindle-checkpoint and chromosome-related outcomes. This supports a focused preclinical comparison, but not a BUBR1-protein, allele-rescue, muscle or clinical claim. The supplement-label spellings above cannot be equated to the study material.

The source-reviewed NAC and L-ascorbic-acid oocyte arms are stronger chromosome-fidelity observations than a generic “antioxidant” hypothesis: in Sod1-deficient oocytes, 500 µM of either compound during 18-hour in-vitro maturation increased the normal 20-chromosome category from 60.0% to 84.6% and reduced premature sister-chromatid separation. However, this was an imposed oxidative-stress/meiosis model at 21% oxygen; neither treatment's effect on BUBR1 was measured, and there was no skeletal-muscle or MVA model. This makes NAC a **high-priority common-compound fidelity comparator**, not a BUBR1-restoration nominee. The project already reviewed this study; joining it to product-label counts adds prevalence context, not a new experiment or independent study. The [source qualification](../ascorbate-oocyte-qualification.md) preserves denominators and assay limits.

Linoleic acid is also visibly represented on supplement labels, but the high-prevalence name does not create an MVA treatment case. The [B7E2 paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12570107/) isolated linoleic acid among extract components and tested purified-component/vinorelbine combinations for cancer-cell apoptosis. Its BUBR1 phosphorylation and mitotic-checkpoint-complex measurements were for **B7E2 plus vinorelbine**, not isolated linoleic acid; the effect belongs in a tumor-chemosensitization context, not a muscle-rescue ranking.

The elected Track 2 pair remains **niacin for a BUBR1-restoration test and apitegromab for a distinct downstream muscle-compensation test**, performed separately. Amla is added to the high-priority preclinical comparator set, not as a replacement nomination. The [current route comparison](../decision-report.md) remains the authority for efficacy and nomination; this label layer does not change its evidence grades.

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
