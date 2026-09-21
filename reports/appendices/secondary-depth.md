> Package copy: local links/navigation are adapted for this finite review package. Labels identify underlying files that are not included. Scientific text is retained; approved status/source metadata changes are recorded separately; original and packaged hashes are recorded in PROVENANCE.

> **Retained dated method detail.** This appendix preserves the original analysis or proposal at its stated date. The current [Track 1 report](../jeweled_track1_report_submission.md), [Track 2 report](../jeweled_track2_report_submission.md) and [methods form](../jeweled_methods_form_answers.md) govern interpretation; see the [revision and evidence scope](../PROVENANCE.md#evidence-scope). Historical delivery or ranking statements are not current status.

# Secondary-panel coding-depth review — 20 September 2026

**Producer complete; independent technical audit passed.** No clinical negative or complete panel-callability result is inferred.

The corrected bounded run completed on Mac Studio 512S2 at 2026-09-20T21:02:17Z; its guard completed at 2026-09-20T21:02:19.198797+00:00 after 88.98 seconds. It examined the frozen 84-gene panel’s all-transcript CDS union: **2,112 disjoint windows and 393,880 ACGT reference bases**. These are annotated CDS targets, not a clinically selected transcript/exon set or every reportable variant position. The particular gene unions do not overlap, so their base totals equal the panel union. All four original read groups contribute to **245,400 extracted alignment records**.

## Panel depth and quality sensitivity

| Quality view | Mean depth | Zero-depth bases | Bases ≥10 | Bases ≥20 | Bases ≥30 | Genes with bases <20 |
|---|---:|---:|---:|---:|---:|---:|
| MQ0 / BQ1 | 43.81 | 0 | 393,845 | 392,617 | 375,770 | 8 |
| MQ20 / BQ20 | 42.46 | 0 | 393,495 | 391,131 | 369,203 | 25 |
| MQ60 / BQ20 | 42.11 | 2,834 | 390,460 | 388,308 | 366,466 | 29 |

At MQ20/BQ20, 391,131/393,880 bases (99.30%) have depth at least 20, but 25 genes contain bases below 20. At MQ60/BQ20, 388,308 bases (98.59%) reach 20; 29 genes contain bases below 20 and **2,834 bases have zero retained depth**. Strong panel-wide percentages therefore conceal local gaps.

**All strict-view zero-depth bases** occur in BMPR1A (12), FLNC (195), MYH7 (209), PMS2 (461) and TTN (1,957). These are losses of retained aligned evidence under the stated mapping-quality filter, not proof of a deletion, absent allele, reference genotype or correct locus assignment. The MQ20/BQ20 and MQ60/BQ20 views hold the base-quality threshold fixed. Comparing MQ0/BQ1 with either also changes BQ, so that comparison alone cannot attribute a loss solely to MQ.

## Complete gene table

The table gives all 84 selected genes, without suppressing low-depth rows. The three X-linked targets are reported at their observed depth, without diploid rescaling. Depth 20 is a descriptive bin, not a validated sensitivity or clinical threshold. Complete three-mode values and transcript IDs are retained in [the compact summary — underlying file not included](../PROVENANCE.md#evidence-scope) and [producer summary — underlying file not included](../PROVENANCE.md#evidence-scope).

| Gene | Chr | CDS bases | Mean MQ20/BQ20 | Min MQ20/BQ20 | Zero MQ20/BQ20 | <20 MQ20/BQ20 | Min MQ60/BQ20 | Zero MQ60/BQ20 | <20 MQ60/BQ20 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ABCD1 | X | 2789 | 23.15 | 6 | 0 | 684 | 6 | 0 | 686 |
| ACTA2 | 10 | 1332 | 42.02 | 25 | 0 | 0 | 25 | 0 | 0 |
| ACTC1 | 15 | 1314 | 43.58 | 27 | 0 | 0 | 27 | 0 | 0 |
| ACVRL1 | 12 | 1865 | 44.21 | 19 | 0 | 2 | 18 | 0 | 7 |
| APC | 5 | 9124 | 45.25 | 22 | 0 | 0 | 22 | 0 | 0 |
| APOB | 2 | 13815 | 39.04 | 17 | 0 | 39 | 17 | 0 | 39 |
| ATP7B | 13 | 4627 | 44.07 | 23 | 0 | 0 | 23 | 0 | 0 |
| BAG3 | 10 | 2019 | 46.74 | 14 | 0 | 16 | 14 | 0 | 16 |
| BMPR1A | 10 | 1871 | 47.84 | 26 | 0 | 0 | 0 | 12 | 182 |
| BRCA1 | 17 | 5893 | 50.66 | 25 | 0 | 0 | 25 | 0 | 0 |
| BRCA2 | 13 | 10662 | 43.54 | 21 | 0 | 0 | 21 | 0 | 0 |
| BTD | 3 | 2337 | 42.46 | 23 | 0 | 0 | 23 | 0 | 0 |
| CACNA1S | 1 | 5905 | 41.62 | 24 | 0 | 0 | 24 | 0 | 0 |
| CALM1 | 14 | 780 | 42.22 | 27 | 0 | 0 | 27 | 0 | 0 |
| CALM2 | 2 | 1126 | 45.39 | 10 | 0 | 4 | 10 | 0 | 4 |
| CALM3 | 19 | 982 | 45.52 | 22 | 0 | 0 | 22 | 0 | 0 |
| CASQ2 | 1 | 1711 | 40.76 | 22 | 0 | 0 | 22 | 0 | 0 |
| COL3A1 | 2 | 4516 | 38.93 | 21 | 0 | 0 | 21 | 0 | 0 |
| CYP27A1 | 2 | 1840 | 43.77 | 21 | 0 | 0 | 21 | 0 | 0 |
| DES | 2 | 2545 | 44.61 | 22 | 0 | 0 | 22 | 0 | 0 |
| DSC2 | 18 | 2844 | 47.35 | 25 | 0 | 0 | 25 | 0 | 0 |
| DSG2 | 18 | 4339 | 43.27 | 19 | 0 | 1 | 19 | 0 | 1 |
| DSP | 6 | 8931 | 43.50 | 21 | 0 | 0 | 21 | 0 | 0 |
| ENG | 9 | 2711 | 45.99 | 27 | 0 | 0 | 27 | 0 | 0 |
| FBN1 | 15 | 8746 | 44.50 | 22 | 0 | 0 | 22 | 0 | 0 |
| FLNC | 7 | 8263 | 40.93 | 19 | 0 | 1 | 0 | 195 | 271 |
| GAA | 17 | 3165 | 43.80 | 24 | 0 | 0 | 24 | 0 | 0 |
| GLA | X | 1660 | 25.07 | 13 | 0 | 236 | 13 | 0 | 236 |
| HFE | 6 | 1311 | 36.83 | 25 | 0 | 0 | 24 | 0 | 0 |
| HNF1A | 12 | 2117 | 45.74 | 28 | 0 | 0 | 28 | 0 | 0 |
| KCNH2 | 7 | 3864 | 38.40 | 12 | 0 | 92 | 12 | 0 | 92 |
| KCNQ1 | 11 | 2232 | 37.91 | 6 | 0 | 109 | 6 | 0 | 109 |
| LDLR | 19 | 3532 | 47.41 | 24 | 0 | 0 | 24 | 0 | 0 |
| LMNA | 1 | 3169 | 48.50 | 26 | 0 | 0 | 26 | 0 | 0 |
| MAX | 14 | 1136 | 43.61 | 23 | 0 | 0 | 23 | 0 | 0 |
| MEN1 | 11 | 1978 | 44.29 | 18 | 0 | 3 | 18 | 0 | 3 |
| MLH1 | 3 | 2424 | 41.90 | 23 | 0 | 0 | 23 | 0 | 0 |
| MSH2 | 2 | 3768 | 47.82 | 20 | 0 | 0 | 20 | 0 | 0 |
| MSH6 | 2 | 4220 | 43.31 | 18 | 0 | 1 | 18 | 0 | 1 |
| MUTYH | 1 | 1763 | 45.04 | 27 | 0 | 0 | 24 | 0 | 0 |
| MYBPC3 | 11 | 3828 | 43.19 | 19 | 0 | 1 | 19 | 0 | 1 |
| MYH11 | 16 | 6468 | 45.62 | 27 | 0 | 0 | 27 | 0 | 0 |
| MYH7 | 14 | 6215 | 46.11 | 25 | 0 | 0 | 0 | 209 | 307 |
| MYL2 | 12 | 1203 | 44.12 | 20 | 0 | 0 | 20 | 0 | 0 |
| MYL3 | 3 | 1147 | 44.56 | 27 | 0 | 0 | 27 | 0 | 0 |
| NF2 | 22 | 2088 | 48.96 | 30 | 0 | 0 | 29 | 0 | 0 |
| OTC | X | 1370 | 23.83 | 13 | 0 | 342 | 13 | 0 | 342 |
| PALB2 | 16 | 3911 | 45.35 | 22 | 0 | 0 | 22 | 0 | 0 |
| PCSK9 | 1 | 2864 | 39.98 | 19 | 0 | 1 | 19 | 0 | 1 |
| PKP2 | 12 | 3075 | 42.70 | 23 | 0 | 0 | 23 | 0 | 0 |
| PLN | 6 | 159 | 51.04 | 32 | 0 | 0 | 32 | 0 | 0 |
| PMS2 | 7 | 3094 | 43.54 | 4 | 0 | 73 | 0 | 461 | 673 |
| PRKAG2 | 7 | 2340 | 44.85 | 21 | 0 | 0 | 21 | 0 | 0 |
| PTEN | 10 | 2173 | 42.82 | 20 | 0 | 0 | 3 | 0 | 99 |
| RB1 | 13 | 3155 | 44.66 | 23 | 0 | 0 | 23 | 0 | 0 |
| RBM20 | 10 | 3714 | 41.94 | 13 | 0 | 5 | 13 | 0 | 5 |
| RET | 10 | 3593 | 40.57 | 18 | 0 | 4 | 18 | 0 | 4 |
| RPE65 | 1 | 1729 | 42.03 | 17 | 0 | 8 | 17 | 0 | 8 |
| RYR1 | 19 | 15117 | 43.26 | 21 | 0 | 0 | 21 | 0 | 0 |
| RYR2 | 1 | 14979 | 44.12 | 18 | 0 | 3 | 18 | 0 | 3 |
| SCN5A | 3 | 6471 | 38.65 | 17 | 0 | 5 | 17 | 0 | 5 |
| SDHAF2 | 11 | 1144 | 46.99 | 24 | 0 | 0 | 24 | 0 | 0 |
| SDHB | 1 | 1362 | 48.80 | 25 | 0 | 0 | 25 | 0 | 0 |
| SDHC | 1 | 1027 | 44.53 | 17 | 0 | 4 | 3 | 0 | 136 |
| SDHD | 11 | 861 | 45.38 | 27 | 0 | 0 | 4 | 0 | 86 |
| SMAD3 | 15 | 1496 | 40.72 | 19 | 0 | 2 | 19 | 0 | 2 |
| SMAD4 | 18 | 1990 | 44.18 | 21 | 0 | 0 | 21 | 0 | 0 |
| STK11 | 19 | 2683 | 49.51 | 27 | 0 | 0 | 26 | 0 | 0 |
| TGFBR1 | 9 | 1951 | 40.69 | 18 | 0 | 3 | 18 | 0 | 3 |
| TGFBR2 | 3 | 1917 | 42.37 | 26 | 0 | 0 | 25 | 0 | 0 |
| TMEM127 | 2 | 1143 | 41.53 | 22 | 0 | 0 | 22 | 0 | 0 |
| TMEM43 | 3 | 1493 | 45.14 | 27 | 0 | 0 | 27 | 0 | 0 |
| TNNC1 | 3 | 728 | 45.77 | 34 | 0 | 0 | 34 | 0 | 0 |
| TNNI3 | 19 | 897 | 51.76 | 30 | 0 | 0 | 30 | 0 | 0 |
| TNNT2 | 1 | 1883 | 42.57 | 21 | 0 | 0 | 21 | 0 | 0 |
| TP53 | 17 | 1751 | 48.99 | 32 | 0 | 0 | 31 | 0 | 0 |
| TPM1 | 15 | 2637 | 40.12 | 20 | 0 | 0 | 20 | 0 | 0 |
| TRDN | 6 | 2956 | 43.09 | 22 | 0 | 0 | 22 | 0 | 0 |
| TSC1 | 9 | 4040 | 49.98 | 25 | 0 | 0 | 25 | 0 | 0 |
| TSC2 | 16 | 6071 | 48.04 | 25 | 0 | 0 | 25 | 0 | 0 |
| TTN | 2 | 114595 | 40.20 | 5 | 0 | 1110 | 0 | 1957 | 2250 |
| TTR | 18 | 679 | 47.73 | 24 | 0 | 0 | 24 | 0 | 0 |
| VHL | 3 | 1037 | 44.35 | 30 | 0 | 0 | 30 | 0 | 0 |
| WT1 | 11 | 1620 | 45.72 | 20 | 0 | 0 | 20 | 0 | 0 |

## Interpretation and preserved method

The unchanged qualified depth producer uses the pinned reference, Ensembl 116 GFF and samtools. It excludes unmapped, secondary, supplementary, duplicate and QC-failed records; refuses unavailable mapping/base qualities for otherwise eligible records; applies native overlapping-mate suppression; and counts no deleted or reference-skipped base as covered. It uses raw qualities without BAQ. These are aligned-read depth observations, not independently established molecule counts. GQ, genotype likelihood, allele balance, locus assignment and parental phase are different quantities.

All views reuse the same reads and cannot be added as independent evidence. The design does not evaluate noncoding sequence outside CDS, full splice flanks, repeat expansions, balanced rearrangements, contamination or low-fraction/mosaic detection limits. An all-transcript union can include CDS absent from a clinically selected transcript set; gene-specific reporting restrictions remain separate. No change to variant classification, EPCR or treatment proposal follows from these depth measurements. The original 17-gene review remains unchanged and separately documents BUB1B, which is not a target in this 84-gene secondary panel.

The first new guard attempt failed on the completed-caller metadata schema before opening patient alignments or creating a producer result. Its source, qualification and failed launch were preserved. The correction requires the actual typed completed shard-ID list 0–63 and passed 21 native fictional controls. The corrected run reused the original 17-gene producer bytes and ten retained native controls. New extraction used only the reviewed alignment and the frozen public target scope.

The sampled guard stayed below its declared limits. Peak group RSS was 2,640,330,752 bytes; sampled combined output was 39,415,507 bytes. These figures cover this guarded run, not the complete genomic project or a predicted future runtime. The controls and output checks qualify their specific software contract; they do not establish biological or clinical sensitivity.

## Provenance and verification

- [Public target design — underlying file not included](../PROVENANCE.md#evidence-scope) and [public terminal status — underlying file not included](../PROVENANCE.md#evidence-scope).
- [Corrected wrapper controls — underlying file not included](../PROVENANCE.md#evidence-scope) and [metadata-only preflight — underlying file not included](../PROVENANCE.md#evidence-scope).
- [Producer terminal status — underlying file not included](../PROVENANCE.md#evidence-scope), [guard status — underlying file not included](../PROVENANCE.md#evidence-scope), and [pre-extraction binding receipt — underlying file not included](../PROVENANCE.md#evidence-scope).
- [Byte-identical local mirror check — underlying file not included](../PROVENANCE.md#evidence-scope), including all 19 producer output pins.

The separately implemented audit passed at 2026-09-20T21:06:15.696931+00:00, with 55 fresh file bindings. It independently reconstructed all 84 genes, 2,112 windows and 393,880 reference positions, and matched all 245,400 original-versus-selected decoded record payloads, tags and multiplicities. Every per-position, window, gene and panel depth statistic reconciled. All 393,880 positions matched in each of the three technical recounts. The checker passed 40 native fictional controls; its dedicated guard passed eight. The actual audit guard completed in 167.82 seconds, with peak sampled RSS 2,689,925,120 bytes. Its recounts use the same pinned samtools engine and overlap policy, so they are technical reproduction, not an independent depth algorithm or biological validation.

[Completed independent verification — underlying file not included](../PROVENANCE.md#evidence-scope), [terminal audit status — underlying file not included](../PROVENANCE.md#evidence-scope), and [audit guard — underlying file not included](../PROVENANCE.md#evidence-scope).

Independent verification SHA-256: `43982f2f129d0bb6a85fb4448878f7ce9037e8f3a98c61cc90410e94381f72b9`.

Producer summary SHA-256: `614ebc4ad6aff0bda9b910b819843219c3b5f19f36e8093e8735a624cfc76aab`.
