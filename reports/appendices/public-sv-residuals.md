> Package copy: local links/navigation are adapted for this finite review package. Labels identify underlying files that are not included. Scientific text is retained; approved status/source metadata changes are recorded separately; original and packaged hashes are recorded in PROVENANCE.

> **Retained dated method detail.** This appendix preserves the original analysis or proposal at its stated date. The current [Track 1 report](../jeweled_track1_report_submission.md), [Track 2 report](../jeweled_track2_report_submission.md) and [methods form](../jeweled_methods_form_answers.md) govern interpretation; see the [revision and evidence scope](../PROVENANCE.md#evidence-scope). Historical delivery or ranking statements are not current status.

# Saved public structural-variant residual annotation review

2026-09-20. Complete; independently verified at 23:20:31 UTC.

The existing tandem-repeat annotation distinguishes shared misses within several deletion size bins, most clearly the 100–499-base bin with 103 truth records. Insertions remain difficult across annotation strata, particularly at 1,000–9,999 bases. This finite diagnostic supports retaining an explicit structural-variant limitation; it does not identify a single cause or establish that an untested caller would improve performance. Stop this diagnostic here. No third caller or patient structural-variant rerun is launched from this result.

This is a descriptive analysis of saved public HG002 chromosome 20 GIAB v5.0q benchmark files. It adds annotation strata to the previous fixed decisions. The prior totals remain 52 truth alleles matched by both callers, 17 by DELLY only, 119 by SvABA only and 614 missed by both. The union remains 188/802 truth alleles; this is not a new recall result or merged-call precision.

## Shared misses within the existing class and size bins

Each fraction is the number missed by both callers divided by all truth alleles in that cell. TRF is the literal existing tandem-repeat flag. Its absence is **not** evidence that an allele is outside repeats. In particular, 78 truth alleles without the TRF flag still carry a RepeatMasker annotation. The table marginalizes over RepeatMasker state; the complete joint TRF × RepeatMasker tables and all four caller outcomes, including zero cells, are in the saved summary.

| Class | Existing size bin (bases) | TRF flag present | TRF flag not encoded |
| --- | --- | --- | --- |
| DEL | 50-99 | 56/94 (59.6%) | 8/21 (38.1%) |
| DEL | 100-499 | 59/66 (89.4%) | 9/37 (24.3%) |
| DEL | 500-999 | 6/6 (100.0%) | 4/9 (44.4%) |
| DEL | 1000-9999 | 4/5 (80.0%) | 1/14 (7.1%) |
| DEL | 10000-50000 | 0/0 (no records) | 0/2 (0.0%) |
| INS | 50-99 | 76/101 (75.2%) | 10/23 (43.5%) |
| INS | 100-499 | 216/237 (91.1%) | 35/45 (77.8%) |
| INS | 500-999 | 57/67 (85.1%) | 4/5 (80.0%) |
| INS | 1000-9999 | 48/49 (98.0%) | 21/21 (100.0%) |
| INS | 10000-50000 | 0/0 (no records) | 0/0 (no records) |

For 100–499-base deletions, 59/66 with TRF versus 9/37 without TRF are shared misses. The joint annotation cells differ as well: TRF with no RepeatMasker fields has 50/52 shared misses, TRF with a RepeatMasker value has 9/14, and RepeatMasker without TRF has 4/27. Therefore “repeat annotated” is not an interchangeable binary label across these sources. Annotations describe these benchmark alleles; they do not prove why a caller missed them.

The insertion limitation extends across annotation groups. At 100–499 bases, shared misses are 216/237 with TRF and 35/45 without it. At 1,000–9,999 bases, these are 48/49 and 21/21. The literal REMAP values in the latter bin also have high shared-miss fractions: interspersed 17/17, partial 50/51, tandem 1/1 and novel 1/1. The last two cells are too small to support generalization. These counts do not support treating the gap as solely tandem-repeat associated.

All 802 alleles have finite LCR entropy values and literal REMAP annotations. LCR was kept continuous with no post hoc cutoff. The saved ranges overlap between at least one matched group and the shared-miss group in all eight class/size bins containing both matched and missed alleles; no diagnostic LCR threshold or effect estimate was inferred. TRF is present for 625 alleles and RM_clsfam has a value for 301. Missing annotations remain explicit in the complete tables.

## Scope and verification

The run used 12 saved files totaling 1,724,588 compressed bytes. It read the original 197,162-record chromosome context and joined the exact 802 selected alleles using chromosome, position, ID, full REF and ordered ALT. No allele was excluded or renormalized, and the existing DELLY-v2 and SvABA-v4 final decisions and class/size margins were preserved. All selected full REF spans fall within a single official benchmark interval; containment was descriptive, not an additional eligibility filter. Existing annotation presence and typed values were checked against both final truth files.

Genotype allele content agrees under the prior order-independent policy. Exact original and final GT encodings are retained: all 802 original phased separators became unphased separators in both final files (158 1|1→1/1, 333 0|1→0/1, 311 1|0→1/0). Unchanged phase is not claimed and no patient phase conclusion follows.

The producer's 29 focused fictional controls passed on Mac Studio 512S2 with native Python 3.12.14. One actual run completed in approximately 0.55 seconds under a 300-second, 2-GiB memory and 64-MiB output guard. Peak sampled process-group memory was 31,604,736 bytes. All pinned inputs, source, runtime, configuration and guard files were unchanged afterward. No alignments, BAMs, reference sequences, new callers, rematching, comparator reruns or network downloads were used.

A separate standard-library VCF reader reconstructed every full ledger row, the original context ordinals, retained raw genotype/annotation values, headers and all aggregate tables exactly. The [actual-output verification — underlying file not included](../PROVENANCE.md#evidence-scope) completed at 23:20:31 UTC. Its nine fictional helper checks are a separate scope from the producer’s 29 controls; they are not a claim of full assay or clinical qualification. The [report interpretation check — underlying file not included](../PROVENANCE.md#evidence-scope) confirms the table arithmetic, annotation denominators, genotype-encoding counts and range-overlap wording.

## Method decision and limits

Close this annotation diagnostic. Retain the established structural-variant sensitivity limitation and both caller-specific results. The class/size tables reveal descriptive annotation separation, but they do not measure improvement from another method, explain every shared miss, or justify confidence that a patient has no structural variant. Fixed broad size bins do not remove every confounder. This single draft chromosome benchmark contains correlated alleles; no confidence intervals, significance tests or population effect estimates were calculated. Small cells and absent annotations must remain visible.

The submission's causal ranking, EPCR, Track 2 proposal, pitch and frozen package are unchanged by this optional methods result. The exact analysis source is now archived in the [verified private v28 source snapshot](https://github.com/EternaPeptix/mva-hackathon/tree/45bb5c47cd1c2f92893ff6f8965051428965bfb3/workflow/source-snapshots/public-sv-residual-review-20260920). Actual result evidence remains in the local artifacts below. Portable package v1 is an earlier frozen snapshot.

## Reproducible local evidence

- [Result summary — underlying file not included](../PROVENANCE.md#evidence-scope), SHA256 `eaeefd901790f4a7d75039b75282d9b12a4a07c482b16e7c719808ebdcf4435a`.
- [Full 802-row ledger — underlying file not included](../PROVENANCE.md#evidence-scope), SHA256 `509d6016afafbe05b49d77d64e91781be7f545a9c619dded8b7fc4c5e36c220c`.
- [Original and final headers — underlying file not included](../PROVENANCE.md#evidence-scope).
- [Frozen producer source — underlying file not included](../PROVENANCE.md#evidence-scope) and [29 native controls — underlying file not included](../PROVENANCE.md#evidence-scope).
- [Approved configuration — underlying file not included](../PROVENANCE.md#evidence-scope), [exact launch plan — underlying file not included](../PROVENANCE.md#evidence-scope) and [guard completion — underlying file not included](../PROVENANCE.md#evidence-scope).
- [Independent source review — underlying file not included](../PROVENANCE.md#evidence-scope).

The launch plan records the exact native command and every input hash. The reviewed guard preserves existing outputs and requires a new owned output directory for a deliberately authorized reproduction.
