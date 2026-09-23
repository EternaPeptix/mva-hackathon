# Common-name food matches: source roles and remaining questions

Both common-name parts are consolidated without changing their record annotations. **154 source IDs and all 161 FDA common-name query rows are represented: 153 saved abstracts are declared read, and MED:18082453 remains explicitly unread because the saved corpus contains no abstract.** Its title and Comment-on metadata receive no experimental inference or reading credit.

Structural checks passed against the frozen corpus, selection, FDA mention rows, both reviews, both manifests and all 39 literal logs. The logs contain 260,518 abstract characters. Exact strings, hashes, source order, end markers, per-record references and declared chunk IDs match. These checks verify artifacts; they do not prove actual viewing or scientific correctness. The original reading declarations belong to the two part reviewers. The original consolidation preceded the root sample audit. The separate [admission check](../food-admission-findings.md) records twelve predetermined source/annotation comparisons, with no substantive correction; this is targeted checking, not an error-rate estimate.

**Coverage remains separate from admission.** Eight IDs overlap the previous literal-food overlay. The common layer therefore represents 146 IDs absent from that overlay, of which 145 have saved abstracts. Across both food layers, 67 + 153 − 8 = **212 unique saved-abstract IDs**, while 67 + 154 − 8 = **213 total hit IDs**. All **234 distinct query rows** are annotated, including one title/metadata-only row. Source-ID counts are not independent-study counts. This package awards **zero automatic outer-corpus credit and zero new appraisal credit**; the separate [64-ID admission](../food-admission-findings.md) reconciles already-read records against the remaining queue without counting an ID twice.

Saved common-layer status: **0 retraction flags, 2 preprint records, 22 review publications, 2 case-report publications, 3 erratum links on 2 records, and 6 journal-to-preprint links**. The errata are attached to MED:22457777 (saved reference only, no linked ID) and MED:40744332 (MED:40883118 and MED:41807160). The preprints are PPR178023 and PPR624548; only PPR178023 has explicit saved version metadata (version 1 and its version list). Linked notices and versions were not newly read, so their contents/reasons are not inferred.

The original readers classified 23 abstracts as reviews/syntheses, while publication metadata labels 22 as reviews. MED:25844087 accounts for that distinction: its reader classification is preserved alongside its saved research-article/Journal Article types. The combined 213-ID food union has 1 saved retraction flag, 5 selected preprints, 7 erratum links on 6 records, and 4 records with explicit version metadata. The full exact metadata inventory is retained in both JSON outputs; retraction reasons and linked-version equivalence are not invented.

The original role vocabulary is retained: **137 distinct query-role labels across 161 rows**. Exact counts for every role are in `summary.exact_query_role_counts`; counts are not collapsed into a new efficacy taxonomy. Common-name query totals are:

| Query | Rows |
|---|---:|
| ALBUMIN | 6 |
| ASCORBIC ACID | 3 |
| BIOTIN | 7 |
| CARBON DIOXIDE | 1 |
| CASEIN | 6 |
| COLLAGEN | 43 |
| DEXTROSE | 1 |
| DIMETHYL SULFOXIDE | 1 |
| FATTY ACIDS | 23 |
| FOLIC ACID | 2 |
| FORMALDEHYDE | 3 |
| GLYCINE | 6 |
| HYDROGEN PEROXIDE | 8 |
| L-ARGININE | 5 |
| L-CYSTEINE | 4 |
| L-GLUTAMINE | 1 |
| L-LYSINE | 1 |
| LACTIC ACID | 8 |
| NITROGEN | 8 |
| PULP | 3 |
| RIBOFLAVIN | 1 |
| SUCROSE | 6 |
| TAURINE | 2 |
| THIAMINE | 1 |
| UREA | 1 |
| VITAMIN A | 2 |
| VITAMIN D | 2 |
| YEASTS | 6 |

The **29 bounded follow-up IDs** below are questions, not new drug appraisals or nominations. The missing comment is included as a source-recovery question. Original model, compound/derivative/mixture/reagent distinctions, functional positives, nulls and BUB1B limitations remain in the unchanged source records.

MED:16040755, MED:18082453, MED:20951947, MED:22457777, MED:22886373, MED:23936839, MED:25229346, MED:26347199, MED:28223486, MED:29769406, MED:30693443, MED:31067459, MED:31540287, MED:31880094, MED:33068521, MED:34504786, MED:34685725, MED:36638165, MED:39506684, MED:39845481, MED:39990212, MED:40197001, MED:40401507, MED:40610403, MED:40744332, MED:41353334, MED:41803973, MED:42019507, PPR:PPR178023.

The most consequential preserved distinctions include actual unsaturated-fatty-acid hepatocyte rescue with unnamed species at the abstract level; Sod1-oocyte segregation-error reduction by ascorbic acid or NAC without explicit drug-restored BUBR1 localization; folic-acid cell protection versus cow expression markers; collagen biomaterial/islet function versus dietary collagen; and Wy-14643-associated short-term proliferation/protein changes, where fatty acids are downstream measurements and functional recovery was not assigned to the drug. Several source-specific BUB1B effects instead involve genetic perturbation or cancer progression. These observations do not establish MVA rescue or justify promotion over provisional niacin/apitegromab; everolimus remains conditional.

Replay the structural checks with `python3 workflow/track2-drug-census/verify_food_common.py --workdir /path/food-common-abstracts-20260923 --check-only`. Optional path arguments support relocation. The verifier rejects missing or changed identity/metadata/log/end-marker evidence; its sole allowed absent abstract is the explicit zero-reading MED:18082453 record. The saved package includes no raw abstracts; the original literal strings remain only in the part logs.

The [lipid-remodeling primary](../lipid-rescue-qualification.md), [bromate adverse study](../bromate-qualification.md), and [missing fatty-acid comment](../fatty-acid-binding-source-limit.md) are separately qualified. Their full-paper scope does not convert the remaining questions into validated efficacy.

The Wy-14643 flag (PMID 42019507) is already answered by the [same-paper primary appraisal](../protein-and-tissue-qualification.md): drug-associated short-term marker/protein changes, with no drug-assigned mass recovery, function or survival endpoint. The separate PHx recovery series cannot be assigned to Wy. This is a resolved duplicate follow-up, not a new appraisal.

Further primary follow-up retains [taurine](../taurine-qualification.md) as an access-limited muscle-physiology reserve and qualifies [oocyte metabolic timing](../oocyte-metabolism-qualification.md), including adverse ARA exposure and genetic rescue. These source qualifications add no abstract-reading credit. The folic-acid flag (PMID39506684) was already addressed at targeted Methods/Results/caption depth in the [existing appraisal](../catalog-independent-review.md); rediscovery adds no new qualification credit or cow clinical-benefit claim.
