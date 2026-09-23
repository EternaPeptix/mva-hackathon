# Research probes beyond clinical-development catalogs

**23 September 2026.** The [SGC collection](https://www.thesgc.org/chemical-probes)
and [donated collection](https://www.thesgc.org/donated-chemical-probes) add a
source-bounded inventory of research reagents. This broadens discovery beyond
marketed medicines. Individual stock, supplier eligibility, clinical use and
product identity were not established by reading these tables.

## What was enumerated

| Layer | Count and meaning |
|---|---|
| Main SGC table | 250 source rows, 249 distinct probe labels |
| Donated table | 141 source rows, 140 labels; all are also in the main table |
| Union | 391 preserved source rows, **249 probe labels**, not 391 molecules |
| Negative controls | 230 records scoped to their parent probe; 226 distinct nonempty control-cell strings |
| Clinical-catalog name candidates | 12 probe labels; exact name cross-references only, not structure-confirmed identities |
| Curated-appraisal name candidates | Nine probe labels, linked to the current ledger; not nine new efficacy findings |
| Saved-corpus retrieval | Seven mention rows in seven source/ID records; six probe labels, no negative-control hits |

`BAY 1816032` appears twice on each page with conflicting protein-family labels
(Kinase and 2OG); all rows are retained. Its target is **BUB1**, a different
protein from BUB1B/BUBR1. No source target cell literally names BUB1B; that says
nothing about indirect effects or unlisted activity.

The 237 probe labels without an exact clinical-catalog name candidate are
unresolved, not 237 proven new molecules. Even the matched generic label
`CRTH2 antagonist` requires identity qualification. No count is added to the
existing 15,357-entity catalog denominator.

## What the seven abstract readings change

Every complete saved abstract for these seven hits was read. The
[reading overlay](hit-review.json) preserves source hashes and decisions,
separately from the generated [retrieval rows](mentions.tsv).

| Hit | Actual evidence and current decision |
|---|---|
| **GSK343, PMID 27226494** | The subsequent [primary qualification](../ezh2-oocyte-qualification.md) separates genetic EZH2/BubR1 stabilization from drug methylation effects; admits GSK343 and DZNep as mechanism probes. |
| **MS023, PMID 39300069** | [Primary qualification](../ms023-qualification.md) retains ccRCC benefit and actual BUB1B RNA/protein decrease, with normal-cell and selectivity limits. |
| **GSK591, PMID 42121254** | [Primary qualification](../prmt5-mettl3-qualification.md) separates GSK591-only animal effects from cell/PDO combination effects and preserves the second-PDO null comparisons. |
| **Bromosporine, PMID 41663394** | [Primary qualification](../bromosporine-qualification.md) retains monotherapy/combination benefit but does not establish selective catalytic PCAF action, formal synergy or BUBR1 restoration. |
| **JQ1, PMID 28061448** | Abstract reports enhanced tumor-cell effects with volasertib; enantiomer unspecified. Do not silently assign this to (+)-JQ1. |
| **A-485, PMID 42211570** | Review article, not a new primary experiment; follow its underlying sources. |
| **BSP, PMID 21526452** | **False bromosporine match:** BSP denotes bisulfite PCR. The actual tested agent is now separately [qualified as lunasin](../lunasin-qualification.md), a BUB1B-transcript research reserve. This does not establish a dietary-soy effect. |

These readings add **zero** to accepted outer-corpus completion counts because
this is a separate, potentially overlapping pass. The primary appraisals are
counted once in the main decision ledger, not once per table or alias.

## Retrieval limits and replay

The [stdlib parser](../sgc_probe_inventory.py) checks the exact source table
headers, retains duplicates and controls, and validates catalog references.
Identity matching preserves punctuation, code separators and stereochemistry.
The positive and negative JQ1 enantiomers remain separate; bare JQ1 is flagged
as unspecified. The explicitly listed BSP acronym stays ambiguous. LP99's
`2S,3R-enantiomer` control is a contextual description and is not searched as a
global molecule. Only the explicit Homer control pair is split into two names.

The screen covers the 7,362 saved source/ID records, titles and abstracts only,
using literal preferred labels and a few explicitly documented aliases. It
misses alternate spellings, unlisted synonyms and compounds named only in full
text. **GSK620 illustrates this limitation:** its separately reviewed GBM
paper supplies actual drug evidence even though this narrow name screen did
not retrieve it. A zero-hit label remains unassessed, not ineffective.

Run from the repo root with the externally retained source files:

```sh
python3 workflow/track2-drug-census/sgc_probe_inventory.py \
  --main-html /path/sgc-probes.html \
  --donated-html /path/sgc-donated.html \
  --corpus /path/literature-gene-unrestricted.json \
  --out /path/replayed-sgc-inventory
```

[summary.json](summary.json) records source, script, decoder and output hashes.
The raw HTML pages and licensed article payloads remain external. Fresh source
downloads can change the enumeration; they are not expected to reproduce these
bytes unless their recorded hashes match. Regression tests cover enantiomer
separation, code/gene boundaries, intact inequality text, duplicate annotations,
contextual controls, changed headers and dangling catalog references.

This source collection covers research probes, not all purchasable compounds,
food constituents, industrial chemicals or every medicine worldwide. Those
remaining spaces and the unread literature require further work.

The current replay against **222 assessments** preserves all four SGC output tables and **nine labels with curated-appraisal name candidates**. The three common-compound additions add no probe-inventory matches or efficacy credit.
