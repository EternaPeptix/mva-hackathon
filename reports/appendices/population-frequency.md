> Package copy: local links/navigation are adapted for this finite review package. Labels identify underlying files that are not included. Scientific text is retained; approved status/source metadata changes are recorded separately; original and packaged hashes are recorded in PROVENANCE.

> **Retained dated method detail.** This appendix preserves the original analysis or proposal at its stated date. The current [Track 1 report](../jeweled_track1_report_submission.md), [Track 2 report](../jeweled_track2_report_submission.md) and [methods form](../jeweled_methods_form_answers.md) govern interpretation; see the [revision and evidence scope](../PROVENANCE.md#evidence-scope). Historical delivery or ranking statements are not current status.

# Exact-key population-frequency resource audit — 20 September 2026

The single approved 25-key public-database audit completed its full guarded run on Mac Studio 512S2 at **20:19:13 UTC**, with **16 keys containing ALFA count rows excluded by the earlier gnomAD-only export, seven found keys containing no stored frequency rows, and two exact keys absent**. All 25 original prefix-only responses were reproduced exactly. The database and reader inputs retained their pinned bytes before and after the query. No original variant inventory, model, genotype or annotation was changed.

This resolves the principal uncertainty raised by the [saved-evidence triage — underlying file not included](../PROVENANCE.md#evidence-scope): for most keys, the empty export came from source-selection policy. It was not evidence that the database lacked all population information.

## Retained population evidence

There are **134 stored frequency rows, all ALFA**, across the 16 keys. Every row has positive AC and AN getters. The table gives the stored ALFA_TOT ratio and the maximum individual ALFA row ratio separately. Population groups can overlap, so counts are never pooled; no ancestry is inferred. ALFA_TOT is one stored source label, not a total calculated here. All full row values, source enum numbers, schema and presence/default evidence remain in [all-frequency-sources.json — underlying file not included](../PROVENANCE.md#evidence-scope).

| Exact normalized key | Stored ALFA_TOT AC/AN | Maximum ALFA row AC/AN (source) |
|---|---|---|
| `20:63488237:CCG:C` | 6319/11862 (53.2709%) | 53/94 (56.3830%; ALFA_SAS) |
| `21:13615550:G:T` | 2493/14050 (17.7438%) | 2406/9690 (24.8297%; ALFA_EUR) |
| `21:44539674:A:T` | 22530/34192 (65.8926%) | 75/92 (81.5217%; ALFA_AFO) |
| `21:44539675:T:G` | 23351/35280 (66.1876%) | 94/114 (82.4561%; ALFA_AFO) |
| `21:44539676:C:A` | 21876/33060 (66.1706%) | 89/108 (82.4074%; ALFA_AFO) |
| `21:44591443:T:C` | 618/20356 (3.0360%) | 572/15036 (3.8042%; ALFA_EUR) |
| `21:44592276:A:G` | 12/20840 (0.0576%) | 1/116 (0.8621%; ALFA_AFO) |
| `21:46334535:G:A` | 94/31874 (0.2949%) | 4/100 (4.0000%; ALFA_AFO) |
| `21:46334557:G:A` | 81/14781 (0.5480%) | 4/102 (3.9216%; ALFA_AFO) |
| `22:15700206:G:A` | 289/17132 (1.6869%) | 263/12706 (2.0699%; ALFA_EUR) |
| `22:15708029:T:C` | 1838/6392 (28.7547%) | 1835/5350 (34.2991%; ALFA_EUR) |
| `22:18606984:C:T` | 6281/16288 (38.5621%) | 4914/12036 (40.8275%; ALFA_EUR) |
| `22:23982904:G:A` | 275/1806 (15.2270%) | 66/206 (32.0388%; ALFA_OTR) |
| `22:42140642:T:C` | 2779/28202 (9.8539%) | 13/38 (34.2105%; ALFA_OAS) |
| `22:42141242:G:A` | 1398/20148 (6.9387%) | 14/34 (41.1765%; ALFA_OAS) |
| `22:50064040:A:AGCACCCCCACCCCACAGGCCACTCACCTCCCCG` | 1187/25297 (4.6923%) | 288/2756 (10.4499%; ALFA_OTR) |

Every one of these 16 keys has a maximum stored ALFA AC/AN getter ratio above 0.001. **Fifteen of the original 22 caller-specific model entries, representing 12 distinct gene/allele-set identities, have at least one such component allele.** This is a separate descriptive resource view; it is not a rerun or reclassification of the frozen gnomAD-only model tables. It also does not measure the frequency of a combined sequence.

For PCNT, both original SNVs have retained population counts: ALFA_TOT 94/31,874 and 81/14,781 (about 0.295% and0.548%), with maximum subpopulation rows 4/100 and 4/102. Their previous empty gnomAD-only exports should no longer be described as absence of public population evidence. The GATK alternative insertion and local sequence ambiguity documented in the triage remain separate questions.

EEF1A2, the three KRTAP10-1 components, both KRTAP10-6 and POTEH components, POTED, both CYP2D7 components, GSTT2, MLC1 and the second RIMBP3 component likewise have stored ALFA counts. Caller overlap is not independent replication; counts for KRTAP10-1 components cannot be assigned to the literally equivalent constructed joint sequence. No benign/pathogenic or causal conclusion follows from these counts alone.

## Remaining unresolved resource information

| Exact normalized key | Result |
|---|---|
| `20:34091776:T:TTTGG` | Exact key absent |
| `21:26844273:GAGCCCCTTCGTCCTCGCCCTC:G` | Exact key absent |
| `21:26844297:T:G` | Exact key found; no stored frequency rows |
| `21:43418483:A:C` | Exact key found; no stored frequency rows |
| `22:18606898:G:A` | Exact key found; no stored frequency rows |
| `22:23973345:T:A` | Exact key found; no stored frequency rows |
| `22:39686094:A:C` | Exact key found; no stored frequency rows |
| `22:39686097:G:C` | Exact key found; no stored frequency rows |
| `22:39686098:A:C` | Exact key found; no stored frequency rows |

These nine keys remain unresolved. All three CACNA1I keys were found but have no stored frequency rows; the prior clustered-sequence and unknown-phase limitations remain. EIF2S2 and the previously reviewed ADAMTS1 deletion were absent as exact keys, which does not distinguish a representation gap from resource omission. The ADAMTS1 SNV, SIK1, DDT and first RIMBP3 allele have found-but-empty records. None is assigned AF 0 or proved rare by this audit.

At the model-entry level, the seven entries with no ALFA-positive component above 0.001 are EIF2S2, the three CACNA1I pairs, ADAMTS1, SIK1 and DDT. This list is an accounting remainder, not a candidate ranking. The already completed ADAMTS1 review is not repeated; EIF2S2 retains its multi-ALT caller disagreement, CACNA1I retains one five-base cluster, and SIK1/DDT retain their original DP 2/DP 1 limitations.

## Presence, provenance and limits

All five fields in the installed frequency protobuf have implicit scalar presence. Positive AC/AN getter values are nondefault decoded count evidence. By contrast, all 134 `frequency` scalar getters and homozygote-count getters are zero: those zeros **cannot distinguish explicitly stored zero from absent/default**. The frequency scalar was not interpreted as AF, and zero homozygote counts were not treated as confirmed absence. Even when AN is positive, AC getter zero in a hypothetical record would still lack explicit presence information. The 18 fictional controls qualify that distinction; no actual queried row had AC 0.

Reserialized decoded message bytes and unknown fields are preserved, but they are not claimed to be the original MVStore wire encoding. The assay queried only the fixed 25 exact keys in the existing database; it did not normalize aliases, retrieve external data, inspect reads, reannotate variants or estimate clinical probability. Resource coverage, population callability and locus assignment are not validated by this audit.

Native controls: 18 passed, 0 failed/skipped, database unopened. The result receipt was timestamped at 20:18:59 before final input rehashing; worker terminal status completed at 20:19:12 and the enclosing guard at 20:19:13 UTC. Approved actual execution: one attempt, 27.46 seconds, peak sampled owned-group RSS 301,957,120 bytes. The guard used nice+10, 600 seconds, 4 GiB group RSS, 64 MiB output and minimum 40 GiB free, and bound Java/LIRICAL plus all 62 class-path dependencies. Both database hashes matched `191585391252f5b9f86a6fefa2b891feda8058bb0f273a661de9e101c8104c1a`. The original caller had completed naturally and was not restarted or signaled.

Exact source/config/command and qualification are retained in [PREPARATION.md — underlying file not included](../PROVENANCE.md#evidence-scope), which remains a preserved pre-run document. The [execution receipt — underlying file not included](../PROVENANCE.md#evidence-scope), SHA `1ac10afaafc19baed9f68e82489dcef6b5e810986acb44c15c096eae9c3d0041`, binds every input and output. [Terminal status — underlying file not included](../PROVENANCE.md#evidence-scope), SHA `211584a09b18dceeae7047245b8cca806ac0367ef1235363aa69288f86ad667b`, links that receipt; [guard status — underlying file not included](../PROVENANCE.md#evidence-scope), SHA `48a1059695e2bcd5927f6e1f5a9466b6a6a0aeb3cac0d9ce87b777b9986c2893`, confirms successful bounded completion. Raw output SHA is `7f9e9111e973260687962d591ec8fb26b2135ff66232f92ca695addee5dc438c`; the exact per-allele/model arithmetic is in [bounded-summary.json — underlying file not included](../PROVENANCE.md#evidence-scope).

The next decision should use this separate all-source evidence together with the preserved representation and quality limits. No additional read analysis is required to establish why these 16 exports were empty. No existing evidence archive or frozen checkpoint was reopened.
