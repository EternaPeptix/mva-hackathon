# Public drug census and BUB1B evidence screen

Read [the decision report](decision-report.md) for new findings and the current
publication-hold status.

**Snapshot: 21 September 2026, America/Phoenix (retrieval timestamps are UTC).**
This is an enumerated, source-bounded census and retrieval screen. It does not
estimate clinical efficacy or establish that the search contains every medicine
worldwide. The manually reviewed decisions are separate from database coverage.
Missing annotations and unmatched names remain **UNKNOWN**.

## Measured scope

- **16,784 ChEMBL clinical-development molecule records**, including all **364**
  withdrawal-flagged records returned by the separate withdrawal query. These
  reduce to **14,659 parent groups**, preserving every source molecule ID.
- Union with **2,497 DrugCentral IDs** present in at least one FDA, EMA or PMDA
  approval download adds **223** unmatched approved-list records: **14,882 main
  catalog rows**. The crosswalk preserves **151 name/structure conflicts** rather
  than silently combining them. Some unresolved rows can represent the same drug;
  this is not a certified count of unique active ingredients.
- **475 additional DrugCentral records have unverified eligibility/status** and
  remain in a separate table. They include possible medicines outside these three
  regulatory lists. They are included in the literature and interaction screens.
  No chemical or database membership alone was converted to an approved medicine.
- Thus **15,357 catalog/supplement entities** receive an explicit literature
  screening state. The main catalog carries **156,295 source aliases**. The
  literature matcher considers **92,792 normalized aliases**, including the
  status-unverified supplement.
- The retrieval scanned **all 7,561 ChEMBL mechanism records**, **all 18,552 ChEMBL
  targets** and **all 19,378 DrugCentral interaction records**. There are **433**
  DrugCentral interaction rows for status-unverified supplemental entities.
  Database-wide target rows are screened before retaining only referenced targets
  in the compact derived target table.
- The literature corpus contains **7,362 Europe PMC records**, with abstracts for
  **7,191**. A separate gene-in-title/abstract query contains **1,460 records**.
  The lexical screen finds **503 entities** with a sufficiently specific alias
  co-mention (495 main-catalog, eight supplemental), producing **1,779 rows**;
  **401 rows** also contain the gene name in title/abstract. The other **14,854
  entities remain UNKNOWN**. **6,160 ambiguous co-mention rows** are retained
  separately. These are retrieval counts, not counts of successful drug studies.

The source versions and every downloaded response's URL, SHA-256 checksum,
retrieval time and available HTTP modification date are in `provenance.json`.
Machine-readable totals, including annotation coverage, are in `coverage.json`
and `literature-coverage.json`. Do not substitute the above snapshot for current
product-specific regulatory review. `approval-overrides.json` and
`fda_freshness.py`, maintained separately, check selected later approval updates.
The FDA supplement retains **40 2026 novel-approval products**: 36 have direct
name cross-references, two combination products have separately matched
components, and two new products (daraxonrasib and milsaperidone) remain
unmatched to this older catalog and receive a complete-corpus alias screen.
These are supplemental product records, not additions to the **15,357**
parent/supplement screening denominator; product/component equivalence must
not be assumed from a name cross-reference.

## Inclusion, identity and status rules

`catalog.py` takes the union of ChEMBL `max_phase >= 1` and
`withdrawn_flag = true`. Phase 4 is an approval **record**, not proof that a
product is currently marketed, appropriate for this child or authorized in a
particular jurisdiction. Phases 1–3 can include discontinued programs; phase
alone does not show current trial activity. Veterinary flags are retained: the
catalog is deliberately broader than human-only medicines.

ChEMBL `molecule_hierarchy.parent_chembl_id` groups related forms; every original
ID and source record remains available. This is not an independently verified
active-moiety or formulation equivalence. Cross-source identity uses a unique
exact **full InChIKey**, or a unique exact normalized **preferred name** when
known structure keys do not conflict. No fuzzy brand-name matching or automatic
salt equivalence is used. Ambiguous and conflicting matches remain distinct.

The DrugCentral regulator files contain 2,331 FDA, 456 EMA and 435 PMDA rows.
They are **unioned individually**. The site's `FDA+EMA+PMDA` file is an
intersection according to its published SQL and is deliberately not used as a
union. DrugCentral's source exports have different dates: 2023 SDF, 2021
SMILES/interaction data and unversioned regulator CSVs. The source site currently
labels its full SQL dump 2023; that 1.40 GB compressed dump was not downloaded.
The compact sources provided usable identifiers, names and interaction data.
Sixteen interaction-export IDs absent from the structure/approval exports are
preserved: seven join clinical parents and nine enter the unverified supplement.

## What the target screen establishes

The live reviewed UniProt record confirms **BUB1B / BUBR1**, **O60566**,
**ENSG00000156970**, human taxon 9606. ChEMBL target **CHEMBL4295998** maps to
that accession. BUB1 (O43683) is a different protein. The optional neighborhood
contains only same-organism interactors explicitly annotated in the UniProt
BUB1B record; it is neither a complete biological network nor a therapeutic model.

There were no curated clinical BUB1B mechanism rows in the retrieved ChEMBL and
DrugCentral mechanism sources. The all-activity lookup found **24 BUB1B assay
records**, of which **two** map to census clinical drugs:

- Dasatinib: **residual activity >75% at 1 μM**, not an IC50 or evidence of useful
  BUBR1 restoration.
- Ceritinib: **IC50 13,080 nM**, an inhibitory biochemical result, not mutant rescue.

The underlying assay descriptions state unknown organism origin despite the
human target mapping. The full source rows preserve relations, units, descriptions,
document IDs and validity fields in `bub1b-clinical-drug-activities.tsv`.
The neighborhood's nicotinamide/suramin SIRT2 interaction annotations do not
establish a beneficial direction. Experimental compounds outside the clinical
filter, such as BRT-1, require separate appraisal; zero clinical database hits
cannot erase them or imply that BUB1B is undruggable.

## Literature matching and manual appraisal

`literature_match.py` scans **titles and available abstracts**, not the full text
of every article. The unrestricted Europe PMC query can retrieve records that
mention the gene only in references/full text; an explicit field marks whether
the title/abstract itself contains a gene literal. Punctuation and case are
normalized, and token boundaries prevent substring hits inside unrelated words.

Short/code-ambiguous aliases, aliases mapping to multiple entities, common
biology/chemistry terms, and trade-only or untyped synonyms are assigned to a
separate ambiguity table. This prevents ordinary-word brands such as *Search*,
*Equal* and *Success* from becoming automatic drug evidence. The ambiguity policy
is explicit and inspectable, not a complete linguistic disambiguation model.
Even specific generic-name hits can describe controls, assay reagents, harmful
exposures, another disease, endogenous substances or background discussion.
No hit means only **no specific catalog alias in the retrieved title/abstract**.

`decisions.json` contains the separately authored source appraisals.
`catalog-appraisals.tsv` joins their compound aliases to every screened entity
without inventing efficacy estimates. Multi-compound route groups and identity
collisions are flagged for review, and their disposition is not silently
inherited by all matches. An exact alias join does not strengthen the evidence
in the original manual decision. Root-owned `appraisal.md`, `literature-review.md`
and the report route table supply the biological assessment.

## Bounded completion audit

[`literature_triage.py`](literature_triage.py) gives all 1,779 specific co-mention
rows an explicit triage state, including the **401 gene-literal rows across 279
source records**. It uses transparent abstract-context routing rules plus
separately marked manual inspection of restorative, same-gene, injury and
identity-confusion cases. It does **not** claim that all papers received manual
full-text review. Remaining ambiguity is recorded as unresolved; titles alone
are not used to infer intervention outcomes. Preprints, patents and published
articles remain separately identifiable, not automatically independent evidence.

A final bounded manual pass read all **28 initially unresolved gene-literal rows**
across **21 source records**. Every record had an abstract; none was left unread
because of a missing abstract. Their compound roles resolve to conditioning or
injury exposures, culture/assay/background uses, computational nomination,
cytotoxic perturbation, an engineered genetic rescue, or related paclitaxel-risk
patent claims. No additional restorative intervention emerged from this set.
This resolves the abstract role, **not** full-text efficacy or clinical validity:
those remain unestablished. The patent abstracts in particular do not establish
validation of their risk-prediction claims. Duplicate AGR/MED study records and
related patents are explicitly identified.

The table distinguishes oncology cytotoxicity/resistance, computational or
biomarker associations, experimental cell-cycle perturbations, assay reagents,
other-lesion interventions, clinical care for MVA complications, false lexical
matches and unresolved roles. Specific corrections include Silver–Russell
syndrome versus silver, a generic *cofactor* versus folitixorin, decitabine versus
its deoxycytidine subphrase, CAPE versus caffeic acid, and genetic/model-dependent
rescue versus the reagent used to run that experiment.

The finite pass flags icariin (aged WT mice), vitamin K2 (neonatal anesthetic
injury; subtype unresolved), CAPE (a distinct ester with reduced BubR1 during
neurodifferentiation), NAC (Sod1-deficient oocytes) and griseofulvin (reported
Mad1/BubR1 binding with cancer-cell killing) for precise primary appraisal.
These annotations establish neither actual-mutant muscle rescue nor a new
successful medication. The 503-entity broad list was additionally filtered for
aging/muscle/protective contexts; the resulting audit table includes unrelated
lesions and injury exposures, which are not promoted by beneficial language.

## Files

| File | Role |
|---|---|
| `catalog.tsv` | Every main clinical/approved parent or unmatched DC record, with IDs, status and annotation counts |
| `chembl-clinical-records.tsv` | Every original ChEMBL molecule/formulation record and its status |
| `aliases.tsv` | Source-qualified main-catalog names and synonyms |
| `drugcentral-status-unverified.tsv` | Every supplemental DC record, aliases and explicit status uncertainty |
| `drugcentral-crosswalk.tsv` | Exact identity matching and unresolved conflicts |
| `chembl-mechanisms.tsv`, `chembl-targets.tsv` | Catalog-associated mechanisms and compact referenced target identities |
| `drugcentral-interactions.tsv` | All downloaded interaction rows, with main/supplement scope |
| `bub1b-target-neighborhood.tsv` | Direct-target or annotated-interactor matches; no efficacy inference |
| `bub1b-clinical-drug-activities.tsv` | The two clinical-drug assay records, preserving original units/relations |
| `literature-mentions.tsv` | Specific alias co-mentions for manual triage |
| `literature-ambiguous-mentions.tsv`, `literature-alias-policy.tsv` | Ambiguous hits and the applied screening rules |
| `literature-drug-coverage.tsv` | Every main/supplement entity, including explicit UNKNOWN rows |
| `catalog-appraisals.tsv`, `decision-alias-joins.tsv` | Traceable manual-decision joins and identity ambiguity |
| `literature-triage.tsv`, `literature-triage-summary.json` | Abstract-context routing of all specific co-mentions; manual restorative/identity checks explicitly distinguished from rule-based context |
| `literature-restorative-context-audit.tsv` | Finite additional audit of broad aging/muscle/protective-language contexts; not a rescue list |

Array/object cells are JSON within tab-separated files. Blank fields are missing
source values, not zero. Files are UTF-8. Stable IDs should be used for joins.

## Reproduce

Python 3.9 or newer. Catalog, matching, triage and verification use the standard
library; the separately reconstructed literature fetcher also requires `requests`. Run from the repository root. A new
cache directory obtains a fresh snapshot; an existing cache intentionally reuses
its checksummed source responses. Do not mix source versions by partially editing
the cache. The script checks page totals/identities and ChEMBL version stability.

```sh
python3 workflow/track2-drug-census/fetch_literature.py \
  --out /path/to/new-literature-snapshot --date 2026-09-21
python3 workflow/track2-drug-census/fetch_literature.py \
  --verify-only --out /path/to/existing-literature-snapshot
python3 workflow/track2-drug-census/catalog.py --cache /path/to/catalog-cache
python3 workflow/track2-drug-census/catalog.py --offline --cache /path/to/catalog-cache
python3 workflow/track2-drug-census/literature_match.py \
  --corpus /path/to/literature-gene-unrestricted.json \
  --title-abstract-corpus /path/to/literature-gene-title-abstract.json
python3 workflow/track2-drug-census/literature_triage.py \
  --corpus /path/to/literature-gene-unrestricted.json
python3 workflow/track2-drug-census/fda_freshness.py \
  --corpus /path/to/literature-gene-unrestricted.json
# Optional FDA replay input: --html /path/to/saved-fda-page.html
python3 workflow/track2-drug-census/target_lookup.py /path/to/target-lookups.json
python3 workflow/track2-drug-census/verify_census.py
```

The local source cache for this run is outside the repository at
`/Users/jeweled/mva-hackathon/track2-drug-census-work/catalog-cache`; literature
corpus snapshots are alongside it. The source queries are recorded in
`literature-coverage.json` and
[`fetch_literature.py`](fetch_literature.py), reconstructed from the original API
calls. The fetcher checks complete cursor pagination and unique record counts,
preserves existing snapshots, and records the retrieval recipe. Europe PMC
queries intentionally have no date restriction; the independent PubMed search
uses the stated publication-date upper bound. A refreshed live index may return
different counts, so the retained snapshot/checksums remain the dated evidence.
Verification of the retained literature snapshots passed (1,460 / 7,362 Europe
PMC records and 1,379 PubMed IDs). A fresh end-to-end network replay also
returned the same counts, plus 74 additional PubMed XML records, in the external
`track2-drug-census-work/reproduction-check-20260921` directory. An offline
catalog rebuild succeeded. Raw caches and run logs are excluded from Git.
All work used public metadata on the local Mac; no Spark compute or experiments
were used.

## Attribution and source terms

- [ChEMBL API](https://www.ebi.ac.uk/chembl/api/data/docs), **ChEMBL 37**, released
  1 May 2026. [Source license: CC BY-SA 3.0](https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/general-questions).
- [DrugCentral public downloads](https://drugcentral.org/download), from the
  University of New Mexico DrugCentral team. [Source license: CC BY-SA 4.0](https://drugcentral.org/privacy).
- [UniProt O60566](https://www.uniprot.org/uniprotkb/O60566/entry), entry version
  234, annotation updated 2 September 2026. [Source license: CC BY 4.0](https://www.uniprot.org/help/license).
- [Europe PMC REST API](https://europepmc.org/RestfulWebService), queried by the
  literature workflow. Article titles/abstracts retain their individual source
  rights; this screen redistributes identifiers, titles and match metadata rather
  than claiming ownership of article text.

Retain source attributions and applicable share-alike terms when redistributing
these database-derived tables. These source terms do not purport to relicense
unrelated repository content.
