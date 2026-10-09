> Package copy: local links/navigation are adapted for this finite review package. Labels identify underlying files that are not included. Scientific text is retained; approved status/source metadata changes are recorded separately; original and packaged hashes are recorded in PROVENANCE.

> **Retained dated method detail.** This appendix preserves the original analysis or proposal at its stated date. The current [Track 1 report](../jeweled_track1_report_submission.md), [Track 2 report](../jeweled_track2_report_submission.md) and [methods form](../jeweled_methods_form_answers.md) govern interpretation; see the [revision and evidence scope](../PROVENANCE.md#evidence-scope). Historical delivery or ranking statements are not current status.

# Native secondary allele inventory — 2026-09-20

The independently verified inventory retained one exact called-ALT identity per caller, solely
in the haplotype/genotype-only review stratum. No exact called-ALT identity entered
a clean aggregate P/LP stratum in the frozen resources. This is not a clinical negative, panel-callability result or a
reportability decision.

## Scope and accounting

The inventory uses 22 completed autosomal source files per caller, plus one
separately completed fixed-ploidy X/Y file per caller: 46 original files in total.
It preserves the actual X/Y model and unavailable GT encodings. Source status,
original ordinals, full native records and normalization lineage are retained.
Noncoding concrete exact matches are included. FILTER, GQ, DP and called-ALT
participation are descriptive fields and do not exclude a matching record.

| Caller | Original records | Normalized records | Nonconcrete normalized records | Exact matching occurrences | Unique called-ALT identities |
|---|---:|---:|---:|---:|---:|
| GATK | 4,929,998 | 5,029,069 | 18,381 | 1 | 1 |
| bcftools | 4,943,443 | 5,012,648 | 7 | 1 | 1 |

These are caller-specific record counts from the same specimen/alignment context,
not independent samples or unique biological variants pooled across callers.
Both retained identities were fully called; there were no duplicate normalized
matching identities and no matching identity without a called ALT.

The frozen public resource contains 4,471,671 records, of which 335,093 mention
panel genes. The retained public inventory has 63,316 source records and exact
normalized identities, all from biallelic originals. Another 120 signal records
outside the primary/concrete scope were counted explicitly. The panel has 84
genes. The resource filename is `clinvar_20260913.vcf.gz`, its snapshot label is
`20260913`, and its literal header is `fileDate=2026-09-13`; these identify the same
SHA-pinned resource, not different downloads.

## Retained record and its limits

The retained normalized key is `3:38603929:T:C`, in **SCN5A**, public record
`48289` (`rs1805124`). In this dated public record the individual-allele
`CLNSIG` is `Benign/Likely_benign`. Selection occurred because its separate
`CLNSIGINCL` contains `440849:Pathogenic`; the method therefore keeps it in
`haplotype_or_genotype_only_review`. That included assertion does not become a
pathogenic classification of the standalone allele or establish the included
haplotype in this case. This allele is not nominated as a secondary finding on
present evidence.

| Native field | GATK | bcftools |
|---|---|---|
| GT | `1/1` | `1/1` |
| DP | 46 | 46 |
| AD | `[0,46]` | `[0,46]` |
| GQ | 99 | unavailable |
| PL | `[1858,138,0]` | `[255,138,0]` |
| FILTER | unfiltered `.` | unfiltered `.` |

Unavailable GQ remains unavailable; unfiltered `.` is not recoded to PASS. The
caller records are preserved separately. Sequence identity and matching displayed
GTs are not independent biological validation, parental phase or clinical
reportability evidence.

Absence of another exact called-ALT match applies only to these explicit native
sources, this frozen panel/resource and this concrete exact-allele inventory.
Variant-only record absence is not a reference genotype. The analysis does not
establish complete 84-gene callability, sensitivity to novel/structural/repeat or
mosaic variants, current clinical classification, disease-specific restrictions
or inheritance/reportability. The separate coding-depth work does not establish
84-gene callability here.

## Verification and reproducibility

The immutable producer passed 35 native fictional controls with no skips,
failures or errors. Public-only preflight separately exercised the actual
frozen resource's header, assertion selection and normalization before the
corrected case run. The first case attempt failed before matching at a literal
date-format gate; the first public-only attempt failed on an undeclared optional
CLNALLE field. Both attempts, all source versions and their controls are retained.
The corrections preserve strict metadata identity and explicit unknown values.

The corrected producer completed at `2026-09-20T19:07:58Z`. Its guard completed
with exit 0 after about 331 seconds, peaking at 585,498,624 bytes of group RSS and
440,123,353 bytes of new output. Caller 79898 and its four current children were
still live at the end.

The separate independent implementation passed at `2026-09-20T19:14:36Z`,
reconstructing original records, full native/public contexts, normalization
lineages, descriptive quality fields and all occurrence/group counts. It checked
175 fresh file bindings across all 46 native sources and reproduced the two
caller-specific occurrences. Its own 36 native fictional controls and 10 native
guard controls passed. It uses the same pinned bcftools normalization engine;
independent code and invocation are not an independent normalization algorithm or
biological validation. Its guarded run completed with exit 0 after about 313
seconds, with peak sampled RSS of 560,693,248 bytes and no resource-limit violation.

Internal artifact appendices use absolute private paths so this document can be
integrated into `repo/report` without breaking its links:

- [Producer summary — underlying file not included](../PROVENANCE.md#evidence-scope), SHA256 `49e9c55ebb28c10cc36f4da2909b1414e0c6b24ff76d35a31f991c05506d139a`.
- [Identity groups — underlying file not included](../PROVENANCE.md#evidence-scope), SHA256 `32be7d6b01fe070d3e083b2ba50b35cde391e200c1b6f1460025b7e8da4d23c8`.
- [Full occurrence evidence — underlying file not included](../PROVENANCE.md#evidence-scope), SHA256 `47d956270967a50f0cbceebcd8f5c2654b72edd0285d40047ed161a161d35f8e`.
- [Source censuses — underlying file not included](../PROVENANCE.md#evidence-scope), SHA256 `8aaa0f97bd5f0ea408abea06a1046e5ff03ef9c3708a9818f429e9ded03ea7a6`.
- [Public-only preflight — underlying file not included](../PROVENANCE.md#evidence-scope), SHA256 `52e02994aed5a055a643836594500e957b78597be897e9db5de6425b5398e7fb`.
- [Producer controls — underlying file not included](../PROVENANCE.md#evidence-scope), SHA256 `9c88971c0f73a52287687e067e07dcc08d798f34c3c94cd983208c3562ed805a`.
- [Run guard — underlying file not included](../PROVENANCE.md#evidence-scope), SHA256 `3b91b55b164891cd06fee305da6626aac31c2d70f31933aba1f13f591441841e`.
- [Independent verification — underlying file not included](../PROVENANCE.md#evidence-scope), SHA256 `d6aaedde3adc302db5e410beb59ee1ab47081d9033c2f00119a61f64a7bd215a`.

Producer source SHA256: `5e190155a01a9dc6144c123f79e6630a4110e8d062862c5b247115189d7e3786`.
Independent verifier SHA256: `ff438cd0c5aa14b066b83645ea09b56bcff2dc198dd8198191f3c48556413185`.
The individual case artifacts remain private. No new database, genotyper,
alignment extraction or publication was performed by this inventory.
