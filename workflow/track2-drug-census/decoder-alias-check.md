# Saved decoder comparison — 22 September 2026

**Confirmed decoder bug; no historical files changed.** The arbitrary `<[^>]+>` expression in `literature_match.normalize` can begin at a mathematical less-than sign and delete intervening prose through the next greater-than sign, including the end of a later formatting tag. The restricted named-tag expression in `broad_hub_screen` preserves that prose.

The input contains **7,362 source records: 7,191 nonempty abstracts and 171 missing abstracts**. Both functions were extracted directly from the saved source by AST, without executing their module entrypoints. All abstract slots were compared. The existing 92,792-row alias policy was held fixed: 42,592 specific aliases were screened; 50,200 ambiguous aliases were excluded. No aliases were reclassified or rebuilt.

- Legacy normalization exactly reproduces the saved 1,779 catalog–record rows **and their matched alias sets**.
- New normalization gives **1,811 rows: 32 additions, zero losses**, comprising 33 alias–record additions across 11 source records. Mitomycin and mitomycin C are two aliases of the same entity in one record.
- All gains occur with titles held at legacy normalization. Titles normalize differently in 56 records, but those title differences change no specific matches here.
- A controlled ablation changes only generic tag removal to the restricted tag pattern, retaining the legacy unescape order and number of passes. It reproduces the new specific-alias match sets for **all 7,362 records**. Thus the observed match delta is attributable to generic tag stripping.
- Abstract normalization changes in 260 records. Legacy removal finds 338 non-format spans in 263 records; 334 spans in 260 records contain **15,696 normalized tokens** preserved by the new decoder. Four other spans have zero retained tokens. These are token counts, not unique words.

## Exact source-level alias gains

All records below are MED/PubMed. No specific alias losses occur. Catalog IDs, DOI and source identities appear in the compact JSON.

| PMID | Newly recovered specific aliases under the unchanged policy |
|---|---|
| 40847407 | eosin, l pyroglutamic acid, pyridoxine |
| 39512465 | refametinib, selumetinib, trametinib |
| 37959410 | pembrolizumab |
| 32228488 | alvocidib, nelarabine |
| 36685838 | erlotinib, lapatinib |
| 19002265 | cisplatin |
| 38143740 | bortezomib, bosutinib, cyclophosphamide, doxorubicin, durvalumab, elesclomol, fluorouracil, gefitinib, gemcitabine, ixabepilone, methotrexate, mitomycin, mitomycin c, olaparib, paclitaxel, temozolomide |
| 34771604 | docetaxel, epirubicin |
| 26676187 | nocodazole |
| 25023465 | substance p |
| 25564371 | eplerenone |

Direct-gene flags independently change from 1,455 to 1,460: MED 41983084, 39952875, 33417716, 19139432; PPR170445. No direct-gene flags are lost. This is a lexical flag count, not the size of the separate gene-query corpus.

## Interpretation and follow-up

The recovered matches remain **co-mentions**, not administered-drug experiments, rescue evidence, or newly qualified candidates. This audit does not review their intervention roles; some recovered terms can be endogenous molecules or assay reagents. The legacy decoder needs a separately versioned repair, regression cases for inequalities and formatting, and downstream alias/triage reconciliation before any historical counts are superseded. Preserve the current historical files and counts until that work is performed.

**The cause of the separate manual 3,598-record annotation errors remains UNKNOWN.** This reproducible lexical-decoder bug does not establish that the manual readers used the function, that it caused their errors, or that correcting it validates those annotations.

Outputs: `decoder-alias-delta.json` is the compact identity-level result with SHA-256 hashes of the raw corpus, both source files, saved alias policy and historical mentions. `decoder-alias-details.json` adds changed-source identities, normalized-text hashes and deleted-span offsets/hashes/counts, without raw abstracts. `audit_text_decoder.py` reproduces the comparison. All five input hashes were identical before and after the run. No repository files or historical counts were modified.

Replay with `python3 workflow/track2-drug-census/audit_text_decoder.py --corpus /path/to/literature-gene-unrestricted.json --out /path/to/external-audit`. Detailed metadata is regenerated beside the compact output, outside the historical tables.

The subsequent [source-role review](decoder-source-role-review.md) reads all 11 recovered records and distinguishes actual exposures from predictions and background. It establishes no new functional BUB1B or muscle rescue; eplerenone remains a primary-source question.
