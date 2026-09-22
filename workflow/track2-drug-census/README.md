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
The extended [`bub1b-all-target-activities.tsv`](bub1b-all-target-activities.tsv)
retains all **24 rows, 15 molecule identifiers and seven documents**, including
experimental compounds. The [source review](target-assay-review.md) verifies a
four-row residual-activity/inhibition error, with two compound/concentration
transpositions, and retains TG101209's thermal signal as a bounded probe.
Annotations are preserved; corrections do not overwrite source data.
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

## Coverage-gap follow-up

Two recovered blind spots are explicit additions to the original screen:

- [`pubmed_supplement.py`](pubmed_supplement.py) now screens the 74 fetched
  PubMed XML records. There are **48 additional source PMIDs**, yielding **7,410
  source records** in the union. Four are preprints with journal updates already
  in the base corpus; these are not four additional independent studies. The
  supplement has 19 missing abstracts and supplies one previously missing
  abstract for an existing PMID. Six specific and 30 ambiguous alias rows are
  retained. [Three full-text mechanism appraisals](pubmed-supplement-review.md)
  distinguish tumor killing/resistance and developmental disruption from rescue.
- [`gene_record_screen.py`](gene_record_screen.py) inventories **all 1,460 gene
  title/abstract records without requiring a catalog drug-name hit**. A bounded
  reviewer pass inspected 594 prioritized titles plus all missing-abstract titles
  (609 unique titles), read 21 full abstracts, and verified seven primary
  methods/results passage sets plus one additional primary abstract. These are
  this pass's counts, not additive independent studies or reviews of all papers.
  The [appraisal](catalog-independent-review.md) recovers administered melatonin,
  folic acid and experimental interventions that the alias gate missed. None
  currently establishes a reason to displace the niacin/apitegromab comparison.

The record inventory preserves every unadjudicated entry. Machine priority,
manual title inspection, complete abstract reading and primary-source verification
are separate fields. The fixed manual ledger is bound to the exact reviewed
snapshot: a changed input is rejected rather than inheriting a historical review.
The 28 missing abstracts, unprioritized records, full-text-only interventions and
unindexed work remain limitations. This is a sensitivity pass, not a completed
systematic review of all interventions. The original 7,362-record alias tables
retain their original denominator; supplemental tables do not silently rewrite it.

## Complete saved-abstract follow-up

The [new review](abstract-review-findings.md) reads the gene-query records in
source order without an alias gate. All **1,432 available saved abstracts** have
been read. That pass read nine missing-abstract bodies and left 19 title-only;
a [separate recovery pass](missing-records-followup.md) subsequently read nine
more short bodies and one newly recovered patent abstract. Among the 28
originally missing abstracts, current depth is therefore **18 complete bodies,
one new abstract, two partial texts and seven title-only records**.
The [original coverage summary](gene-abstract-review-summary.json)
preserves that pass's actual completed reading; the
[record table](gene-abstract-review.tsv) separates complete saved abstracts,
accessible bodies, title-only records and unread entries. Three AI reviewers
read disjoint batches; this is not three independent reviews of each paper.
Original assessments are in `abstract-review-notes.json`. The compiler validates
source identities and exact abstract hashes; it does not generate judgments.
Supplemental recovery assessments are in `missing-records-followup.json`; they
retain original identities and depths separately. The [primary lead qualification](lead-qualification.md)
raises Aurora-B balancing for functional comparison and qualifies geraniin,
TG101209, KP372-1 and okadaic acid without claiming clinical efficacy.
Source truncations, duplicate versions and unanswered full-text questions remain
explicit. Complete saved-string reading is not an exhaustive worldwide efficacy
review. Geraniin and other newly retained routes have explicit comparative
decisions in the linked findings.

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
| `bub1b-all-target-activities.tsv`, `target-assay-source-checks.json` | All target rows plus separately recorded original-workbook corrections |
| `literature-mentions.tsv` | Specific alias co-mentions for manual triage |
| `literature-ambiguous-mentions.tsv`, `literature-alias-policy.tsv` | Ambiguous hits and the applied screening rules |
| `literature-drug-coverage.tsv` | Every main/supplement entity, including explicit UNKNOWN rows |
| `catalog-appraisals.tsv`, `decision-alias-joins.tsv` | Traceable manual-decision joins and identity ambiguity |
| `literature-triage.tsv`, `literature-triage-summary.json` | Abstract-context routing of all specific co-mentions; manual restorative/identity checks explicitly distinguished from rule-based context |
| `pubmed-supplement-records.tsv`, `pubmed-supplement-summary.json` | All supplemental PMID identities, version links, missing abstracts and union counts |
| `pubmed-supplement-mentions.tsv`, `pubmed-supplement-ambiguous.tsv` | Supplement screened with the original alias policy |
| `gene-record-screen.tsv`, `gene-record-screen-summary.json` | Every gene-query record with separate priority, title/abstract review and verification states |
| `abstract-review-notes.json`, `gene-abstract-review.tsv`, `gene-abstract-review-summary.json` | Authored assessments and source-bound compilation of the deeper reading; no raw abstracts |
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
python3 workflow/track2-drug-census/pubmed_supplement.py \
  --corpus /path/to/literature-gene-unrestricted.json \
  --pubmed-xml /path/to/literature-pubmed-additional.xml \
  --pubmed-search /path/to/literature-pubmed-esearch.json \
  --title-abstract-corpus /path/to/literature-gene-title-abstract.json
# Replay of the reviewed snapshot; changed snapshots require fresh manual review.
python3 workflow/track2-drug-census/gene_record_screen.py \
  --input /path/to/literature-gene-title-abstract.json
python3 workflow/track2-drug-census/compile_abstract_review.py \
  --corpus /path/to/literature-gene-title-abstract.json
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


## Beyond the gene-in-abstract boundary

The [additional review](outer-corpus-findings.md) reads 180 selected
complete abstracts in two nonoverlapping passes (60 plus 120) from the 5,902 records outside the 1,460-record gene-title/
abstract query. That outer set has 5,759 saved abstracts; the unselected records
remain unreviewed by this pass. Selection uses muscle, progeroid and chromosome
fidelity contexts without a clinical-catalog alias requirement. The
[selection ledger](outer-corpus-selection.json) records exact IDs and source
hashes; [original assessments](outer-corpus-review.json) distinguish drug,
genetic, harmful-challenge and unspecified interventions. The
[419-row discovery ranking](outer-corpus-candidate-ranking.json) is metadata
triage, not 419 reviewed experiments. The overlap audit refers to the named
pre-integration commit, preserving its historical meaning. The
[next selection](outer-next120-selection.json) and
[next assessments](outer-next120-review.json) preserve the subsequent fixed-order
120-record pass. These are records, not globally deduplicated studies; 5,722 outer
records remain outside both selections.

These records add concrete routes including cortical-tension modulation,
progerin-specific chaetocin and hesperetin/CISD2. Primary follow-ups are labeled
separately from abstract assessments. No source reading is converted into a
numerical MVA efficacy estimate. `verify_census.py` reconciles saved selection,
identity, depth and coverage claims; it does not validate biological efficacy or
independently reproduce human/animal experiments.
