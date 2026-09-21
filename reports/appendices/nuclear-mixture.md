> Package copy: local links/navigation are adapted for this finite review package. Labels identify underlying files that are not included. Scientific text is retained; approved status/source metadata changes are recorded separately; original and packaged hashes are recorded in PROVENANCE.

> **Retained dated method detail.** This appendix preserves the original analysis or proposal at its stated date. The current [Track 1 report](../jeweled_track1_report_submission.md), [Track 2 report](../jeweled_track2_report_submission.md) and [methods form](../jeweled_methods_form_answers.md) govern interpretation; see the [revision and evidence scope](../PROVENANCE.md#evidence-scope). Historical delivery or ranking statements are not current status.

# Nuclear mixture model output — 20 September 2026

The single approved Mac Studio 512S2 assay completed, and the independent saved-output check passed at **21:30:26 UTC**. GATK returned the boundary contamination value **0.0**, with its native error field **4.76837158203125 × 10⁻⁷**. This is a limited model output. It does **not** establish that the sample is unmixed, pure, clinically negative or below a validated detection limit. The output was not the native insufficient-information fallback (0,1).

This adds nuclear mixture model evidence separately from the earlier mitochondrial NO/ND result. It does not validate individual variants, determine cell fractions, establish donor identity or resolve mosaic aneuploidy. The raw native error is retained without converting it into a clinical confidence interval or biological detection limit.

## Counts and scope

All 37,267 public-panel markers were selected before reading the sample. There was no selection using patient genotypes. The native table and independent reconciliation gave:

| Marker category | Sites |
|---|---:|
| Fixed public panel |37,267|
| Emitted native rows |37,265|
| No native row emitted |2|
| Emitted total count 0 |5|
| Emitted total count 1–10 |14|
| Preliminary total count >10 |37,246|
| Excluded by native relative-coverage gates |97|
| Retained by native relative-coverage gates |37,149|

The categories partition the panel as 37,149+97+14+5+2=37,267. Coverage statistics were computed from the 37,246 rows with total count >10: mean 46.6924501959944 and median 46.0. Retention required total count strictly above 0.5×median and strictly below 3×mean. The table contains 1,739,160 REF+ALT+other-ALT read observations; these are not unique fragments or cells. The two missing rows and five emitted-zero rows remain distinct from reference genotypes.

The 37,149 coverage-retained markers span all 22 autosomes. They are the coverage-filter population; the native model’s internal homozygous-site count and selected estimation strategy were not exported by this CLI and are not inferred here.

| Autosome | Panel | No row | Zero |1–10| Coverage excluded | Coverage retained |
|---|---:|---:|---:|---:|---:|---:|
|1|3,055|2|0|0|5|3,048|
|2|3,170|0|0|2|7|3,161|
|3|2,865|0|0|0|7|2,858|
|4|2,643|0|0|2|9|2,632|
|5|2,389|0|1|0|9|2,379|
|6|2,463|0|0|0|7|2,456|
|7|2,178|0|0|1|5|2,172|
|8|1,865|0|1|0|3|1,861|
|9|1,529|0|1|0|8|1,520|
|10|1,923|0|0|0|8|1,915|
|11|1,880|0|0|1|6|1,873|
|12|1,892|0|1|0|4|1,887|
|13|1,387|0|0|2|2|1,383|
|14|1,204|0|0|0|3|1,201|
|15|1,044|0|1|1|2|1,040|
|16|1,010|0|0|1|2|1,007|
|17|978|0|0|0|2|976|
|18|1,066|0|0|2|3|1,061|
|19|897|0|0|0|2|895|
|20|880|0|0|1|0|879|
|21|512|0|0|0|2|510|
|22|437|0|0|1|1|435|

## Evidence and validation

[Raw native table — underlying file not included](../PROVENANCE.md#evidence-scope), [raw model output — underlying file not included](../PROVENANCE.md#evidence-scope), [producer summary — underlying file not included](../PROVENANCE.md#evidence-scope), [independent verification — underlying file not included](../PROVENANCE.md#evidence-scope), [guard receipt — underlying file not included](../PROVENANCE.md#evidence-scope), and [compact result — underlying file not included](../PROVENANCE.md#evidence-scope) are retained. The 22-row [minor-allele-fraction segmentation — underlying file not included](../PROVENANCE.md#evidence-scope) is a model intermediate, not a copy-number or tumor conclusion.

The standalone checker freshly hashed 30 distinct input/source/output files totaling 30,473,290,114 bytes; checked exact configuration, source, qualification, runtime and command links; independently reconstructed every marker’s membership and coverage category; retained both source and native AF literals while requiring identical finite binary 64 values; and verified sample/header/dictionary and raw result formats. It did not import producer helpers. Its scope is independent saved-output reconciliation, not an independent read counter or estimator algorithm. The [complete 20-file output mirror — underlying file not included](../PROVENANCE.md#evidence-scope) matches fresh remote hashes before and after copying, with no mismatches.

The source and wrapper passed [34 native fictional controls — underlying file not included](../PROVENANCE.md#evidence-scope) and [peer source review — underlying file not included](../PROVENANCE.md#evidence-scope). The earlier failed fictional fixture-copy qualification and its correction remain preserved. The previous [fictional table pilot — underlying file not included](../PROVENANCE.md#evidence-scope) and [public-panel transformation — underlying file not included](../PROVENANCE.md#evidence-scope) were not rerun. None of these controls establishes end-to-end aligned donor-read mixture accuracy.

## Resource use and limitations

The [estimator log — underlying file not included](../PROVENANCE.md#evidence-scope) reports a 96-point segmentation input below its 100-point approximation/window requirement and no changepoint candidates for that input. It does not identify that input by chromosome in those warnings. The 22 reported segment rows therefore provide no evidence that CN/LOH is absent. Both native logs also record fallback from the optional incompatible x86 compression library to JdkDeflater/JdkInflater on ARM; both commands completed successfully.

The guard measured 737.806 seconds end to end, a sampled peak RSS of 1,605,451,776 bytes and a sampled peak result-plus-launch output of 14,375,812 bytes. The two native commands returned 0; the pileup took 708.310 seconds and estimation 2.787 seconds. The declared 90-minute/4 GiB RSS/2 GiB output/minimum 40 GiB-free limits were not reached. Only the new owned process group was controlled. The original caller’s completed metadata identity was unchanged before and after; this is separate from its prior whole-callset verification.

The assay uses the existing qualified GATK 4.6.1.0 runtime, minimum MQ 50, no added per-base BQ threshold and no overlap collapse. Native quality/filter behavior remains as qualified. This estimator is distinct from VerifyBamID2 FREEMIX. The public panel’s mask/spacing/frequency construction is not claimed to be a clinically validated optimal contamination panel. Same-donor DNA, related donors, ancestry mismatch, CN/LOH/allelic imbalance, mosaic aneuploidy, mapping effects and correlated read observations may confound interpretation.

The [approved configuration — underlying file not included](../PROVENANCE.md#evidence-scope) differs from the prepared proposal only in its two execution booleans, as recorded by [root approval — underlying file not included](../PROVENANCE.md#evidence-scope). This later result is separate from the dated checkpoint completed at 21:29:27 UTC, which intentionally excluded the then-active assay. No submission, publication, outreach or laboratory action was performed.
