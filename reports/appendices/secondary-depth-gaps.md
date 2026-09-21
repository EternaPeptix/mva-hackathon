> Package copy: local links/navigation are adapted for this finite review package. Labels identify underlying files that are not included. Scientific text is retained; approved status/source metadata changes are recorded separately; original and packaged hashes are recorded in PROVENANCE.

> **Retained dated method detail.** This appendix preserves the original analysis or proposal at its stated date. The current [Track 1 report](../jeweled_track1_report_submission.md), [Track 2 report](../jeweled_track2_report_submission.md) and [methods form](../jeweled_methods_form_answers.md) govern interpretation; see the [revision and evidence scope](../PROVENANCE.md#evidence-scope). Historical delivery or ranking statements are not current status.

# Secondary-panel depth gaps: locations and transcript membership

The single approved saved-output map completed on Mac Studio 512S2 at **22:11:15 UTC, 20 September 2026**. It located low-depth bases within the unchanged 84-gene, all-transcript CDS design and reconciled the original window/gene/panel depth statistics. No alignment was opened, and no depth calculation or variant caller was rerun.

There are **2,749 bases below depth 20 at MQ20/BQ20**, compared with **5,572 at MQ60/BQ20**, within the same 393,880-base panel union. The strict view has **2,834 zero-depth bases** in five genes. This map makes those locations and their transcript memberships explicit; it does not establish clinical callability, a negative screen, correct locus assignment, phase or variant function. The separate 17-gene analysis, including BUB1B and CRYAA, is outside this map.

[Complete location appendix — underlying file not included](../PROVENANCE.md#evidence-scope) · [All gene/transcript intervals — underlying file not included](../PROVENANCE.md#evidence-scope) · [Original mapped summary — underlying file not included](../PROVENANCE.md#evidence-scope) · [455 membership/category intervals — underlying file not included](../PROVENANCE.md#evidence-scope)

## Same-base comparison

Every base belongs to exactly one cell below. Both modes retain BQ20; only the minimum MQ changes. These are paired views of the same saved read evidence.

| MQ20/BQ20 category | MQ60/BQ20 zero | MQ60/BQ20 depth 1–19 | MQ60/BQ20 depth ≥20 | Row total |
|---|---:|---:|---:|---:|
| Zero | 0 | 0 | 0 | 0 |
| Depth 1–19 | 1,076 | 1,673 | 0 | 2,749 |
| Depth ≥20 | 1,758 | 1,065 | 388,308 | 391,131 |
| Column total | 2,834 | 2,738 | 388,308 | 393,880 |

Thus **2,823 bases move from ≥20 to below 20** under the stricter view, including 1,758 that have zero strict retained depth. No base crosses upward into the ≥20 category. Nevertheless, **63 positions have a numerically higher strict-view depth** within their category. The mapper explicitly retains that saved-output observation; no cause is assigned here, and depths were not forced to be monotonic.

The original 455 intervals split at both depth-category and annotated CDS-membership boundaries. For readable gene-union locations, adjacent positions under one criterion collapse to **342 low-depth runs in the main view** and **357 in the strict view**. The strict zero subset has **32 runs**. These are descriptive segments dependent on the target geometry and threshold, not separate variants or biological events.

## Gene and transcript location guide

All 84 genes and all **2,578 coding transcripts** remain in the machine-readable result, including those without gaps. Twenty-five genes have at least one main-view low-depth base; 29 do in the strict view. The corresponding transcript counts are **317** and **442**. Transcript bases overlap, so these are annotation-membership counts, not independent observations or clinical transcript counts.

The table lists every gene with a low-depth base in either view. `Transcripts` means counts with any low-depth base at MQ20 / MQ60 / all annotated coding transcripts for that gene. The final column is the longest contiguous strict-view low-depth gene-union run, with coordinate-first tie-breaking; it is an entry point to the full appendix, not the complete gap list. **All locations use GRCh38, 0-based half-open coordinates: start included, end excluded.** No clinical transcript was selected.

| Gene | Gene CDS-union bases | <20 MQ20 | <20 MQ60 | Zero MQ60 | Transcripts MQ20/MQ60/all | Longest MQ60 <20 run |
|---|---:|---:|---:|---:|---:|---|
| ABCD1 | 2789 | 684 | 686 | 0 | 7/7/7 | X:153736738–153736849 |
| ACVRL1 | 1865 | 2 | 7 | 0 | 20/20/29 | 12:51920888–51920893 |
| APOB | 13815 | 39 | 39 | 0 | 2/2/9 | 2:21006471–21006489 |
| BAG3 | 2019 | 16 | 16 | 0 | 9/9/10 | 10:119651813–119651820 |
| BMPR1A | 1871 | 0 | 182 | 12 | 0/48/56 | 10:86919260–86919400 |
| CALM2 | 1126 | 4 | 4 | 0 | 1/1/28 | 2:47162171–47162175 |
| DSG2 | 4339 | 1 | 1 | 0 | 22/22/30 | 18:31498252–31498253 |
| FLNC | 8263 | 1 | 271 | 195 | 9/9/9 | 7:128858007–128858221 |
| GLA | 1660 | 236 | 236 | 0 | 16/16/16 | X:101397895–101397931 |
| KCNH2 | 3864 | 92 | 92 | 0 | 10/10/10 | 7:150977837–150977865 |
| KCNQ1 | 2232 | 109 | 109 | 0 | 7/7/9 | 11:2445256–2445324 |
| MEN1 | 1978 | 3 | 3 | 0 | 42/42/42 | 11:64810091–64810094 |
| MSH6 | 4220 | 1 | 1 | 0 | 20/20/27 | 2:47783491–47783492 |
| MYBPC3 | 3828 | 1 | 1 | 0 | 3/3/3 | 11:47346640–47346641 |
| MYH7 | 6215 | 0 | 307 | 209 | 0/56/62 | 14:23419887–23420186 |
| OTC | 1370 | 342 | 342 | 0 | 22/22/22 | X:38411936–38411999 |
| PCSK9 | 2864 | 1 | 1 | 0 | 1/1/14 | 1:55056559–55056560 |
| PMS2 | 3094 | 73 | 673 | 461 | 24/30/70 | 7:5977587–5977790 |
| PTEN | 2173 | 0 | 99 | 0 | 0/6/13 | 10:87965373–87965472 |
| RBM20 | 3714 | 5 | 5 | 0 | 4/4/4 | 10:110644612–110644613 |
| RET | 3593 | 4 | 4 | 0 | 10/10/19 | 10:43114646–43114647 |
| RPE65 | 1729 | 8 | 8 | 0 | 5/5/5 | 1:68448683–68448686 |
| RYR2 | 14979 | 3 | 3 | 0 | 9/9/9 | 1:237503325–237503326 |
| SCN5A | 6471 | 5 | 5 | 0 | 11/11/16 | 3:38566415–38566417 |
| SDHC | 1027 | 4 | 136 | 0 | 19/23/23 | 1:161362409–161362540 |
| SDHD | 861 | 0 | 86 | 0 | 0/5/17 | 11:112094886–112094970 |
| SMAD3 | 1496 | 2 | 2 | 0 | 7/7/20 | 15:67066271–67066272 |
| TGFBR1 | 1951 | 3 | 3 | 0 | 22/22/30 | 9:99105260–99105261 |
| TTN | 114595 | 1110 | 2250 | 1957 | 15/15/15 | 2:178654444–178654534 |

No gene has a zero-depth base in the main view. The strict zero locations can be approached directly as follows:

| Gene | Strict zero bases | Contiguous zero runs | Exact starting location(s), 0-based half-open |
|---|---:|---:|---|
| BMPR1A | 12 | 1 | 10:86919319–86919331 |
| FLNC | 195 | 1 | 7:128858026–128858221 |
| MYH7 | 209 | 1 | 14:23419946–23420155 |
| PMS2 | 461 | 4 | 7:5973398–5973542; 5977587–5977720; 5978623–5978680; 5982864–5982991 |
| TTN | 1957 | 25 | Disjoint runs from 2:178653426–178653507 through 2:178663432–178663444; every run is listed individually in the appendix. The intervening span is not one gap. |

The next useful review, if pursued, can use these exact intervals to establish which frozen transcript CDS segments are involved and where further locus-specific evidence would be needed. Merged gene intervals preserve transcript IDs overlapping any part; they do not imply that every listed transcript spans the entire interval. The full atomic rows retain exact changing membership and GFF line ordinals. This map makes no assertion about canonical transcript choice or biological exon numbering. X-linked genes retain their observed depth without diploid rescaling.

## Evidence limits and execution

The unchanged Ensembl 116 all-transcript policy was reconstructed from the frozen public GFF and FAI and matched to the saved target design. The target union is entirely ACGT. Gene and transcript denominators deduplicate shared CDS positions within each unit; only the panel union counts each genomic position once across the panel. Transcript totals must not be summed to estimate panel coverage. The map does not assess noncoding splice flanks, promoters, UTRs, repeat expansions, structural variants or low-fraction sensitivity. A depth threshold of 20 is a descriptive bin, not a validated clinical threshold.

Every saved target position was required exactly once in each mode; missing output would fail, not become zero. Per-window and gene/panel histograms, counts, sums, means and bins matched the completed depth summaries. This is further descriptive use and reconciliation of the existing measurements, not a new independent read counter or biological experiment. Existing reports, original depths and the separate 17-gene result were not edited.

The method passed **39 fictional native controls** and independent source review before root approved the exact config. Fresh prelaunch checks bound 13 input/receipt files, six source/test files, runtime, qualification and approval. The parent guard ran for **31.65 seconds**, with sampled peak group RSS **2,637,283,328 bytes** and sampled output **16,923,012 bytes**, below its 15-minute/4-GiB/512-MiB budgets. These are measured values for this pass; sampled limits are not instantaneous kernel caps. Eight terminal output files were mirrored with exact hashes, unchanged in native readbacks before and after copying. Only the new owned process group was controlled.

**Dated closure, 20 September 2026:** the [independent saved-input reconstruction — underlying file not included](../PROVENANCE.md#evidence-scope) passed at **22:21:07 UTC**, checking all 84 genes, 2,578 transcripts, 455 mapped rows, paired categories and 63 numeric increases without importing producer helpers. The [readable-output review — underlying file not included](../PROVENANCE.md#evidence-scope) passed at **22:23:18 UTC**, checking all 29 main gene rows, 731 location-appendix rows and exact merged intervals. No corrections were requested. The earlier pending-review report is preserved byte-for-byte in [RESULT-before-independent-review.md — archival document not included](../PROVENANCE.md#archival-material); this closure changes no analysis or result. The original depth audit, method source review and new derived-result review retain their distinct scopes.

## Receipts

- [Approved config — underlying file not included](../PROVENANCE.md#evidence-scope), [root execution review — underlying file not included](../PROVENANCE.md#evidence-scope), [source review — underlying file not included](../PROVENANCE.md#evidence-scope), [39-control receipt — underlying file not included](../PROVENANCE.md#evidence-scope) and [fresh prelaunch checks — underlying file not included](../PROVENANCE.md#evidence-scope).
- [Mapper terminal status — underlying file not included](../PROVENANCE.md#evidence-scope), [guard status — underlying file not included](../PROVENANCE.md#evidence-scope), [eight-file mirror verification — underlying file not included](../PROVENANCE.md#evidence-scope) and [readable-location derivation — underlying file not included](../PROVENANCE.md#evidence-scope).
- Summary SHA256: `db6d50cc8e115cbdca697f581d7ce5a754a419b50ff282276b8f9a52ba828eb4`.
- Original gap-interval SHA256: `b2212a358145a3c70c223e4b2a5ea6c26349b69bfe0e0e0d76d091c10085549c`.

No classification, causal probability, treatment proposal, source publication or submission was changed by this map.
