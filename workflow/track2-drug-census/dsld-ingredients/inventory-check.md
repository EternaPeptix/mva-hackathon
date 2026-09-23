# NIH DSLD ingredient-group inventory and all-alias screen

Checked 23 September 2026. This is a source-defined ingredient-group inventory and automated discovery screen, with **zero new abstract reads, primary-paper reads or efficacy appraisals**. It uses the [official NIH DSLD API documentation](https://dsldapi.od.nih.gov/) and the hash-bound A–Z/Other responses. Ingredient groups are not distinct molecules, products, currently stocked supplements or verified medicines. Source synonym assignments remain identity-unverified.

## Scope and results

| Measure | Count |
|---|---:|
| Exact partitions / source groups | 27 / 6,467 |
| Preserved synonym cells / group-name cells | 78,788 / 6,467 |
| All indexed source cells | 85,255 |
| Distinct nonempty normalized aliases | 69,707 |
| Normalized aliases shared across source groups | 6,228 |
| Saved literature records / nonempty abstracts | 7,362 / 7,191 |
| Group–alias–record mention rows | 41,496 |
| Matched source groups / literature records | 1,143 / 6,578 |
| Source groups with no match | 5,324 |
| Unflagged mention rows / source groups / records | 1,077 / 371 / 774 |
| Unflagged mention rows with direct gene co-mention | 177 |

“Unflagged” means only that no listed lexical ambiguity flag applies. Every mention remains **context-unreviewed**, with efficacy UNKNOWN. Broad terms such as protein, complex, same, stem cells and muscle dominate many matches; short and multi-group aliases are also flagged. All hits are retained. A no-hit group is unknown, not evidence of inefficacy. There are 171 corpus records without saved abstracts; their titles can still match but provide no abstract-reading credit.

Exact normalized-name cross-references reach **1,121** groups in the main catalog, **43** in its status-unverified supplement, **545** in the FDA food table and **60** in the 168-entry appraisal ledger. These overlapping group counts do not add unique molecules to any catalog. FDA names, mixture components, products and source synonym relations are not verified chemical or exposure equivalents.

## Matching and preservation

The script uses the corrected `broad_hub_screen.clean` for titles and abstracts, preserving inequality text and drug names. It applies NFKC, casefolding and Unicode alphanumeric-token boundaries; punctuation and hyphens separate tokens. All overlapping phrases are matched, without substrings inside words. Thus a parent/component term can match inside a hyphenated derivative; the adjoining hyphen is explicitly flagged. Raw text, synonym index, duplicate source cells, category type/case and empty values remain in the inventory. No raw abstracts or snippets are emitted.

There are **6,080 exact group-name strings**, with **368 repeated names spanning 755 IDs**. IDs remain distinct. The source has 199 empty synonym arrays and one empty synonym cell (`dsld:2605:synonym:0`, Eleuthero). It contains no within-group repeated raw synonym string, but 124 additional identical cells arise when group names repeat as synonyms; all are retained. There are 7,971 within-group normalized collisions, also retained.

Category arrays contain **592 JSON null cells**, not the literal string “null”; `botanical` and `Botanical` remain distinct. All 6,467 nutrient-information lists are empty, so they cannot establish a nutrient classification or absence. The saved API docs identify version 9.4.0, while the matching before/after API responses identify **9.5.0, May 2026**, with processed index **2025-09-25n**. The existing guide-visible extraction is empty and was not treated as read prose. Separate partition requests are not a transactional snapshot.

## Checks and replay

All **15 fixtures passed**: corrected decoding/inequality preservation, NFKC/case/whitespace, overlapping phrases, no within-word match, hyphen-component flags, common/short/multi-group terms, inexact/partial totals, duplicate/empty/type-preserving cells, conflicting IDs, ID mismatch, tampered hash, late mutation and a missing partition. Failure tests preserve prior published outputs. Every source group and all 85,255 cells round-trip exactly from the generated tables. Every mention points to valid matching cells. A separate regex matching algorithm agrees on **5,295 alias–source–field hits** for ten declared aliases across all 7,362 records. This checks implementation, not scientific interpretation.

A second full build produced the identical immutable run and every output hash. All inputs are rehashed before publication. Complete files are staged and only then exposed through an atomically replaced relative `derived/current` symlink. Resolve that symlink once and use the immutable directory plus `summary.json` hashes; historical runs remain intact.

```sh
python3 build_dsld_inventory.py --source-dir /path/supplement-ingredients-20260923 \
  --catalog-dir /path/repo/workflow/track2-drug-census \
  --corpus /path/literature-gene-unrestricted.json --out /path/derived \
  --expected-appraisals 168
python3 build_dsld_inventory.py --catalog-dir /path/repo/workflow/track2-drug-census \
  --source-dir /path/supplement-ingredients-20260923 --out /path/derived --self-test
```

This check binds `derived/runs/7a86160f6d9a7faeef69` and a **168-appraisal ledger**. After the parent's separate ledger update, explicitly replay with `--expected-appraisals 170`; this historical check must not be silently relabeled as a 170-ledger check. The portable script uses only Python's standard library and the specified catalog directory's corrected cleaner. Sources, fetch script, retrieval manifest and repository were not modified.
