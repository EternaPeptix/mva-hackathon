# Independent catalog crosscheck: Broad Repurposing Hub

21 September 2026. **This widens discovery coverage without adding an unverified
count of new medicines.** The main census remains 15,357 catalog/supplement
entities. Broad entries overlap that census and include preclinical compounds,
salts, mixtures and dated development-status annotations.

The [official portal](https://repo-hub.broadinstitute.org/repurposing) supplies
drug and sample annotation downloads dated **18 August 2025**. We independently
parsed both, compared exact normalized names and full InChIKeys with the saved
ChEMBL/DrugCentral census, and screened every preferred source name against the
same 7,362-record literature corpus. No supplier-alias expansion, fuzzy identity
matching or chemical-form equivalence was assumed.

| Measured coverage | Result and interpretation |
| --- | --- |
| Drug annotation entries | 7,540 distinct preferred source names; 7,448 normalized forms. Normalization collisions remain flagged. |
| Sample annotations | 22,612 rows across 7,522 actual names and 15,673 nonblank Broad IDs. Repeated IDs and 109 rows lacking an ID prevent calling these unique physical samples. |
| Missing sample groups | 18 drug entries have no sample annotation rows; no structure is invented. |
| No current name/full-key candidate | 2,713 source entries. This is a failed cross-reference, not 2,713 new or eligible medicines. Other mappings also retain collisions, conflicting keys and shared structures. |
| Preferred-name literature screen | 544 source names have 2,148 lexical record mentions; all 7,540 names receive an explicit state. A missing match remains UNKNOWN. |
| Mentions among the 2,713 unmapped entries | 127 names have mentions; 39 have a direct gene co-mention. Seventeen of those pass the stated lexical-specificity policy, which still supplies no efficacy or identity validation. |

**A catalog gap is not necessarily a scientific-review gap.** Unmapped names
with gene co-mentions include CAPE, hesperadin, monastrol, reversine and
S-trityl-L-cysteine, which already have appraisals outside the clinical catalog.
They do not become new discoveries when a second source names them. Conversely,
a compound may occur as a toxicant, assay control, background citation or tumor
killing agent. Every mention remains a discovery pointer until source-level
review establishes the actual intervention and outcome.

## Corrections found by independent reproduction

The sample file contains three unquoted newlines inside vendor-name fields:
two CK-101 rows and one triptorelin row. The parser repairs only those specified
records, preserves each field and checks an independent strict TSV parse. It
does not reinterpret continuation fragments as new compounds.

An early summary inflated distinct sample names to 7,540 by inserting empty
groups during lookups. The corrected value is **7,522**; the 18 empty groups
remain explicit. The original cross-reference assignments were independently
reproduced rather than silently merged.

The literature cleaner also required correction: a generic angle-bracket pattern
could erase text between scientific less-than and greater-than signs. The final
version strips a bounded set of formatting tags while preserving inequality
text. This recovered **38 mention rows**, increasing 536 matched names/2,110
mentions to 544/2,148. These are lexical recall improvements, not newly measured
drug benefits. Original abstracts and earlier reading ledgers are unchanged.

## Inspect and reproduce

- [Source rows and candidate cross-references](broad-hub-screen.tsv): all names,
  name/structure candidate IDs, full source keys, identity flags and search state.
- [Mention locations](broad-hub-mentions.tsv): source identifiers, titles,
  source-record hashes and field-specific token spans; no raw abstract text.
- [Coverage, input hashes and source notices](broad-hub-summary.json).
- [Offline screen](broad_hub_screen.py): optional prototype comparisons are not
  required to reproduce the result from source downloads and the saved corpus.

Obtain the official [drug annotations](https://repo-hub.broadinstitute.org/public/data/repo-drug-annotation-20250818.txt)
and [sample annotations](https://repo-hub.broadinstitute.org/public/data/repo-sample-annotation-20250818.txt).
Use the exact literature snapshot identified in the summary; a later retrieval
can differ. Run from this directory, writing derived outputs outside the checkout:

```sh
python3 broad_hub_screen.py \
  --source-dir /path/to/broad-downloads \
  --catalog-dir . \
  --corpus /path/to/literature-gene-unrestricted.json \
  --title-abstract-corpus /path/to/literature-gene-title-abstract.json \
  --output-dir /path/to/broad-replay
```

The source-only replay and independent matching checks validate parsing,
identity bookkeeping and search locations. They do not validate pharmacology,
current approval, complete worldwide coverage or comparative MVA efficacy.

## Attribution and source notices

Broad Institute Drug Repurposing Hub; Corsello et al., *Nature Medicine* 23,
405–408 (2017), [doi:10.1038/nm.4306](https://doi.org/10.1038/nm.4306).
The current portal identifies compound metadata as
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) for research, with
clinical-treatment/commercial-marketing restrictions in its notice. The dated
TXT banners separately state noncommercial use. Both source-specific notices
are retained in the summary. This is a research-derived crosscheck; the team's
license does not override third-party terms. ChEMBL/DrugCentral mappings retain
the source terms in the [census attribution](README.md#attribution-and-source-terms).
