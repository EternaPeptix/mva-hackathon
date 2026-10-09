# Independent NIH DSLD source check

**The saved 27-partition enumeration reconciles.** All recorded response hashes/sizes match; each A–Z/Other response has an exact (`eq`) total equal to its observed row count. The union contains **6,467 unique group IDs**, with no repeated/conflicting IDs and no missing pages within the reported totals. The largest partition has 838 rows, below the requested 1,000. All `_id` values match `groupId`.

The manifest’s bounded wording—“all documented A–Z/Other by-letter responses retrieved to reported exact totals”—is justified. This is not independent proof of every possible backend record, all molecules, all labels or current availability. No global match-all total was compared. “Other” coverage beyond the documented browse behavior is not formally established; no specific missing record was found.

**One concrete correction:** the 592 category entries are **JSON null elements**, not the literal string `"null"`. For example, group 11 Ashwagandha has `["botanical", null]`. There are zero literal `"null"` category strings. Serializing a Counter keyed by null can misleadingly turn that key into the text `"null"`. Preserve the original arrays and element types.

| Source property | Independently counted |
|---|---:|
| Source group IDs | 6,467 |
| Distinct exact group names | 6,080 |
| Names shared across multiple IDs | 368 |
| Rows bearing those shared names | 755 |
| Synonym cells / distinct exact strings | 78,788 / 72,074 |
| Empty synonym arrays / empty string cells | 199 / 1 |
| Empty nutrientInfo arrays | 6,467 |
| Factsheet link/name entries | 5,392 |

Distinct IDs with the same name must survive: Glycine IDs 92/7364 have different categories, as do Glycerol IDs 1110/8895. This is not a conflicting duplicate-ID defect. The empty synonym is group 2605 Eleuthero, array position 0. Category case variants (`botanical` and `Botanical`), nulls, empty cells and ordering should remain intact in provenance-bearing output.

The [official API documentation](https://dsldapi.od.nih.gov/) describes ingredient **groups**, pagination and a group model with optional synonyms, factsheets and nutrientInfo fields. It is labeled version 9.4.0; the captured API reports 9.5.0, May 2026. The documented category type is `array[String]`, while actual source arrays contain nulls. All other observed source fields have stable types and the same six-field shape.

Both captured version responses are identical and name `dsldnxt_labels-2025-09-25n` and `dsldnxt_ing_groups-2025-09-25n`. Those are index-name strings, not proof of a transactional snapshot or every record’s update date. May 2026 is an API version timestamp, not a current-market validation.

Official ingredient pages call synonyms **“Related Terms on Labels”** and show multiple groups for the same term. These include forms, extracts, plant parts, mixtures, quantities, branded phrases and misspellings—not verified chemical equivalence. Factsheet links include topic searches, not automatically read efficacy studies. All nutrientInfo arrays are empty, so populated CAS/IUPAC/identity enrichment is absent in this snapshot. [Official example](https://dsld.od.nih.gov/ingredient/tin)

The separate [API guide](https://dsld.od.nih.gov/api-guide) response is an SPA shell with no rendered guide text; the saved visible-text file is empty and direct web retrieval likewise supplied none. A public linked script returned 403. No hidden-browser or access-bypass route was used. Thus the actual static API endpoint/model documentation was read; a complete rendered-guide reading is not claimed.

The exploratory `ingredient-example.json` is a correctly logged HTTP 500 HTML error using `method=letter`; it is excluded from the successful 27-response inventory. The successful calls use `method=by_letter`.

The initial source audit preceded derived output. The extension below now verifies a specific immutable run. Only these two check files were written; sources, builder and repository were left unchanged.

Reviewed retrieval-manifest SHA-256: `28409b93b597b4fc1633a47646c188fb530b1250b4a479f0324952cdaf911ad3`.



## Derived-table extension: pinned 168-appraisal run

**All source fields round-trip without loss** in `derived/runs/7a86160f6d9a7faeef69`. The summary SHA-256 is `7a86160f6d9a7faeef69423d63fa05b6d42c50caae6ee782f98b2967e5d34415`; its appraisal input has **168 entries**, digest `de8078cc684a7e0fd1ac3ccf01d834e0391fa1222250c5050965ac0c90f248d4`. The builder input digest is `2235a9d8189a7b81316f818f12a3bd278353e5ab30371beb3c0538f7d25115dd`. All recorded inputs matched at the check, and all eight output-file hashes were verified. This does **not** review or credit the planned 170-appraisal regeneration.

An independent standard-library reader, without importing the builder or its normalization module, reconstructed all **6,467 groups** and **85,255 alias cells** from raw responses. Every group ID/name, category value and type, factsheet entry, nutrientInfo array and ordered synonym was retained. Partition, URL, response file, hit/row indices and hit metadata also matched. The **592 JSON nulls**, category capitalization and the empty `dsld:2605:synonym:0` cell survive exactly. No source values, indices or IDs were lost or changed.

All alias IDs, roles, raw values, ordering and normalization/collision counters matched independently. The 69,708 normalized rows include one empty alias; the empty cell is preserved but correctly unsearchable. Normalization is a search index, not replacement of the original compound/form spelling.

All **6,800 name-crossref rows** match an independently reconstructed exact normalized-name join against the recorded inputs: 5,161 main catalog, 88 status-unverified, 1,480 FDA food and 71 curated-appraisal rows. There are no missing, extra or duplicate rows. Source alias-cell links match, and every row retains `EXACT_NORMALIZED_NAME_ONLY_NOT_CHEMICAL_EQUIVALENCE`. This validates a metadata join, not the chemical identity, formulation, regulatory status or efficacy of the linked source.

All group evidence states remain `UNKNOWN_NOT_ASSESSED_BY_INVENTORY`; all 7,362 corpus-record rows say `NONE_FROM_AUTOMATED_SCREEN`. All 41,496 mention rows retain `LEXICAL_DISCOVERY_CONTEXT_UNREVIEWED`, including the **1,077 unflagged** rows. The pinned summary reports 371 groups and 774 records with unflagged matches. Automatic reading, primary-paper and appraisal credits are zero.

The inspected matcher uses overlapping whole-token sequences and separates fields. Normalized crossrefs require the full name, while raw strings preserve form/mixture detail. Hyphen/common/short/shared flags are incomplete ambiguity aids. **Unflagged does not mean a specific administered compound or beneficial effect.** No scientific roles were assigned, no efficacy was inherited from a name link, and full literature-match recall was not independently rerun.
