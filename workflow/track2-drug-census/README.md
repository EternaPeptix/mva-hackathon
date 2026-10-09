# Public drug census and BUB1B evidence screen

Read [the decision report](decision-report.md) for new findings and the current
publication-hold status.

[PR1 reuse](pr1-reuse.md) records which work is already incorporated, adds
bounded UMK57 control findings, and preserves the later source corrections.

The [food-substance inventory](food-substances/README.md) retains 3,971 FDA rows with explicit source-status flags. The current scientific comparisons are in the [decision report](decision-report.md).

The dated [EPA CPDat v4.0 and v4.1 consumer-product overlays](cpdat-v4.1/README.md) add a product-formulation chemical discovery layer. Their name matches are triage only and do not establish BUB1B efficacy or exposure. Source-qualified findings include high-concentration ibuprofen/naproxen BUB1B-transcript reduction in nonhuman-primate Sertoli cells, an MnBP reduction in pooled human granulosa cells, and an ozone-associated Bub1b reduction in neonatal mouse lung; none is a restoration or MVA-benefit result. Separate [quercetin](quercetin-bubr1b-qualification.md) and [resveratrol](resveratrol-bubr1b-qualification.md) public-array checks add consumer-compound evidence: quercetin gives small mixed-direction HepG2 shifts, while all three resveratrol-treated A549 arrays are below a single control. A [lunasin](lunasin-qualification.md) reanalysis reproduces its published 2.23-fold RWPE-1 BUB1B RNA value from a likely matching GEO matrix and finds a smaller RWPE-2 shift; GEO's citation is missing, the matrix uses one array per condition, and neither result measures BUBR1 function.

The comparison covers **310 appraisals**. [Creatine monohydrate](creatine-inflammatory-myopathy-qualification.md) now has a randomized functional benefit in adults with stable inflammatory myopathy receiving shared exercise; separate aging strength nulls and pediatric results remain part of its assessment. [Sildenafil](sildenafil-dystrophinopathy-qualification.md) did not demonstrate the primary benefits sought in two dystrophinopathy trials, and one produced an imprecise cardiac safety signal. A source-specific [poplar propolis extract](propolis-qualification.md) adds an adverse-direction BUB1B transcript observation at a cytotoxic lymphoma-cell exposure, not a restoration or muscle-benefit result. These findings strengthen compound-level comparison without establishing MVA benefit or changing the provisional nominees.

The [4 October common-compound extension](common-compound-extension.md) adds L-tyrosine, riboflavin, rutin, MK-4 and MK-7, with preparation-specific positives and counterevidence. The independently qualified [D2/inulin/HMB trials](batch11-common-compound-followup.md) add two appraisals and amend calcium-HMB without counting it twice. Controlled adult function signals remain visible alongside dose, analysis and task-test limits; none demonstrates BUB1B/MVA rescue or changes the provisional nominees. The separate outcome-first screen now accepts 2,280 exact source/ID records, with 53,937 unscreened; [batch 17](outcome-first-batch17-screen.md) retains all available saved abstracts, incomplete originals, nulls and alternatives under its completed source-specific review.

The separate [NIH DSLD ingredient inventory](dsld-ingredients/README.md) preserves **6,467 source groups, 78,788 synonym cells and 85,255 total name/synonym cells**. Its **41,496 lexical matches, including 1,077 unflagged matches, provide no efficacy assessment** and confer no reading or appraisal credit. A separate first pass annotates six rows across four previously read abstracts, without new unique credit. Groups, label terms and matches are not unique drugs, verified chemical equivalents or additions to the main entity denominator. API version 9.5.0 (May 2026) and index-name suffix `2025-09-25n` are distinct source descriptors, not current-market or record-update guarantees.

The separate [full-label prevalence layer](dsld-label-prevalence/README.md) processes a January 2026 NIH DSLD export: **214,780 label records**, including **121,959 marked on-market**, with exact Supplement Facts ingredient spellings retained separately. It identifies consumer-label breadth for common compounds but does not estimate sales, actual exposure or efficacy; generic “Niacin” entries do not identify the formulation used in the clinical niacin study.

The [research-probe inventory](sgc-probes/README.md) adds 249 SGC probe labels
and 230 parent-scoped negative-control records as a separate discovery layer.
Seven complete matched abstracts were read; one BSP match is a PCR abbreviation.
The [GSK620 follow-up](gsk620-gbm-qualification.md) credits actual mouse cancer
survival, while the [EZH2 study](ezh2-oocyte-qualification.md) separates genetic
BUBR1 stabilization from drug inhibition. Those primary follow-ups added four appraisals, not 249 demonstrated treatments,
and no outer-corpus reading credit.

The [food-name review](food-substances/common-abstract-findings.md) retains 212 available abstracts across 213 source IDs. [Batch15](batch15-revalidation-findings.md) adds the final 532 saved strings: **3,598 accepted revalidations, 5,759 combined outer source/ID records and zero pending in the original selection**. This closes that queue, not the broader search or primary-source follow-up. The earlier 2,161 annotations were not independently re-audited here. Root checked 16 source/annotation pairs in the final batch; versions are not independent studies and the 143 initially empty fields retain separate recovery statuses.

**Catalog snapshot: 21 September 2026; corrected text screen: 22 September 2026,
America/Phoenix (retrieval timestamps are UTC).**
This is an enumerated, source-bounded census and retrieval screen. It does not
estimate clinical efficacy or establish that the search contains every medicine
worldwide. The manually reviewed decisions are separate from database coverage.
Missing annotations and unmatched names remain **UNKNOWN**.

The [47-question follow-up reconciliation](batch15-followup-reconciliation.md) links component-level answers to the actual qualified sources. Related compounds, different experiments and source-name matches do not close a question; missing experiments remain explicit.

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
  The [corrected lexical screen](literature-v2/summary.json) finds **511 entities**
  with a sufficiently specific alias co-mention (503 main-catalog, eight
  supplemental), producing **1,811 rows**; **401 rows** also contain the gene
  name in title/abstract. The other **14,846 entities remain UNKNOWN**.
  **6,195 ambiguous co-mention rows** are retained
  separately. These are retrieval counts, not counts of successful drug studies.

The source versions and every downloaded response's URL, SHA-256 checksum,
retrieval time and available HTTP modification date are in `provenance.json`.
Machine-readable totals, including annotation coverage, are in `coverage.json`
and `literature-v2/summary.json`; `literature-coverage.json` is the historical
uncorrected screen. Do not substitute the above snapshot for current
product-specific regulatory review. `approval-overrides.json` and
`fda_freshness.py`, maintained separately, check selected later approval updates.
The historical FDA supplement retains **40 2026 novel-approval products**: 36 have direct
name cross-references, two combination products have separately matched
components, and two new products (daraxonrasib and milsaperidone) remain
unmatched to this older catalog and receive a complete-corpus alias screen.
These are supplemental product records, not additions to the **15,357**
parent/supplement screening denominator; product/component equivalence must
not be assumed from a name cross-reference.

The [4 October FDA refresh](fda-approval-refresh-20261004.md) retains **45 products**: five new rows all have existing catalog cross-references. It adds no catalog entities or efficacy appraisals. The frozen40 annotations remain unchanged; bounded new-label sections and pending primary-paper questions are separate from accepted eligibility screens. Use `fda_freshness.py --output <dated-file>` to preserve the historical layer.

A separate [Broad Repurposing Hub crosscheck](broad-hub-review.md) screens all
7,540 preferred source names against the same corpus, retaining 2,148 lexical
mention rows and unresolved identity/status flags. These are overlapping source
entries, not additions to the 15,357 main denominator. The independent parser
and name screen preserve missing sample groups and scientific inequality text.

The direct [EMA inventory crosscheck](ema-coverage-review.md) adds a separate
layer of **1,568 currently authorised human product rows**. Of these, **180**
have no exact whole-field name cross-reference to the saved catalog layers.
These are unresolved product/identity questions, not 180 new unique medicines.
The 15,357-entity denominator is unchanged, and the EMA layer performs no
BUB1B efficacy assessment. Initial authorisation, later decisions and positive
opinions remain distinct.

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

[`literature_rescreen.py`](literature_rescreen.py) scans **titles and available abstracts**, not the full text
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


The current `literature-v2/` layer preserves mathematical inequalities that the
historical `literature_match.py` tag stripper could delete. The saved alias policy
is unchanged. All 7,362 retained records are reprocessed, gaining 32 specific and
35 ambiguous catalog-to-record matches with no losses. An independent replay
using Python's HTML parser and a separate token matcher reproduced both complete
match sets. The 32 recovered specific rows use the [source-role audit](decoder-source-role-review.md);
other manual annotations are explicitly inherited rather than newly verified.
The 96-row keyword routing table is not the complete primary-follow-up queue;
separate source-role follow-ups remain visible. No manual-reading credit is added.

Replay into a new directory, with the original externally retained corpora:

```sh
python3 workflow/track2-drug-census/literature_rescreen.py \
  --corpus /path/literature-gene-unrestricted.json \
  --title-abstract-corpus /path/literature-gene-title-abstract.json \
  --out /path/new-literature-v2
python3 workflow/track2-drug-census/verify_rescreen.py \
  --corpus /path/literature-gene-unrestricted.json \
  --title-abstract-corpus /path/literature-gene-title-abstract.json
```

The first command refuses nonempty output directories or changed corpus versions.
Input bytes are captured before processing, and concurrent changes cause failure.
The second checks the committed layer against a fresh source replay. Neither
command certifies scientific interpretation. Historical tables and hashes stay
available for earlier evidence chains.

`decisions.json` contains the separately authored source appraisals.
`catalog-appraisals.tsv` joins their compound aliases to every screened entity
without inventing efficacy estimates. Multi-compound route groups and identity
collisions are flagged for review, and their disposition is not silently
inherited by all matches. An exact alias join does not strengthen the evidence
in the original manual decision. Root-owned `appraisal.md`, `literature-review.md`
and the report route table supply the biological assessment.

## Historical 1,779-row triage

The historical [`literature_triage.py`](literature_triage.py) gave 1,779 specific co-mention
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

The finite pass flags icariin (aged mice; genotype unverified), vitamin K2 (neonatal anesthetic
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
- [Broad Repurposing Hub](https://repo-hub.broadinstitute.org/repurposing),
  annotations dated 18 August 2025; Corsello et al.,
  [doi:10.1038/nm.4306](https://doi.org/10.1038/nm.4306). The
  [crosscheck attribution](broad-hub-review.md#attribution-and-source-notices)
  preserves both the current portal and dated file notices.

Retain source attributions and applicable share-alike terms when redistributing
these database-derived tables. These source terms do not purport to relicense
unrelated repository content.


## Beyond the gene-in-abstract boundary

The [earlier review](outer-corpus-findings.md) read 433 selected
complete abstracts in three nonoverlapping passes (60 + 120 + 253) from the 5,902 records outside the 1,460-record gene-title/
abstract query. That outer set has 5,759 saved abstracts; the unselected records
remain unreviewed by this pass. Selection uses muscle, progeroid and chromosome
fidelity contexts without a clinical-catalog alias requirement. The
[selection ledger](outer-corpus-selection.json) records exact IDs and source
hashes; [original assessments](outer-corpus-review.json) distinguish drug,
genetic, harmful-challenge and unspecified interventions. The
[419-row discovery ranking](outer-corpus-candidate-ranking.json) supplied metadata
triage. Every ranked record now has a complete-abstract assessment, but these
include reviews and genetic experiments; they are not 419 drug-efficacy tests. The overlap audit refers to the named
pre-integration commit, preserving its historical meaning. The
[next selection](outer-next120-selection.json) and
[next assessments](outer-next120-review.json) preserve the subsequent fixed-order
120-record pass. The [completion selection](outer-ranked-completion-selection.json)
and [assessments](outer-ranked-completion-review.json) cover the remaining 253
ranked records. These are records, not globally deduplicated studies; **5,469**
outer records remained outside those three selections.

The subsequent [unfiltered review](unfiltered-findings.md) reads another 384
complete saved abstracts under a fixed SHA-256 ordering without a biology,
date, drug-name or favorable-result filter. Those four selections totaled **817**
outer source/ID records, leaving **4,942 originally saved abstracts unread** at that stage.
The original recovery ledger records all **143** missing outer abstracts:
18 recovered abstracts and 53 short bodies read (including 14 availability
notices and one keyword index), 11 long bodies unread, five empty retrieved
bodies and 56 HTTP-500 retrieval failures. Missing access is not negative
evidence. `verify_unfiltered.py` checks coverage and, with the original two
retrieval JSON files, replays the fixed selection and every abstract hash.

These records add concrete routes including cortical-tension modulation,
progerin-specific chaetocin and hesperetin/CISD2. Primary follow-ups are labeled
separately from abstract assessments. No source reading is converted into a
numerical MVA efficacy estimate. `verify_census.py` reconciles saved selection,
identity, depth and coverage claims; it does not validate biological efficacy or
independently reproduce human/animal experiments.

The [intervention/endpoint qualification](intervention-qualification.md) checks
all fourteen primary follow-ups from the final ranked selection. It resolves
drug identities, genetic-versus-pharmacological effects, positive organ outcomes
and control harm; it does not change the historical abstract-reading ledgers or
claim full-text review of all ranked records. At that historical checkpoint, the decision ledger contained
152 curated route/compound appraisals, including research tools and a combination
RNA program rather than 152 approved medicines.

The [unfiltered primary comparison](unfiltered-primary-qualification.md) resolves nine prioritized questions, adding six route/compound appraisals and updating existing NAC and nicotinamide entries. Collagen matrices remain a culture method; sFRP3 resolves to an already-qualified genetic study. Neither enlarges the medicine catalog. These qualifications preserve favorable functional results and source-specific counterevidence; they do not close the remaining search.

The [next unfiltered pass](next576-findings.md) adds 576 complete saved abstracts, bringing that stage to **1,393 source/ID records** with **4,366 originally saved abstracts unread**. It flags 62 source questions; these are not completed primary appraisals. `verify_next576.py` independently replays the fixed selection and content hashes when supplied the original retrievals. The [primary comparison](next576-primary-qualification.md) adds ursolic acid, MG132 and DHMMF at explicit evidence depths. It also revises the existing bortezomib appraisal for combination-specific contracture benefit and growth/mortality tradeoffs, without changing the 96-entry total. Six conference abstract texts and one primary letter form a separate recovery overlay; no entire conference book is claimed as read.

The [intervention-identity follow-up](intervention-identity-followup.md) adds three route/compound appraisals and revises the existing STLC route, bringing that stage of the ledger to **99**. That primary follow-up did not add abstract readings; those stage counts are historical.


The [NQO1/SIRT2 follow-up](nqo1-sirt2-followup.md) upgrades NMN's reference
rationale with reported mitotic-completion benefit and adds an adverse dicoumarol
appraisal. The [EMA identity follow-up](ema-identity-followup.md) qualifies three
additional products without finding a supported MVA nomination. Together these
bring the ledger to **103 route/compound appraisals** before the subsequent mebendazole qualification. The recovered
37 names/identifiers and explicit typography variants remain unmatched; zero
retrievals are not negative efficacy. Neither follow-up changes the catalog
entity denominator or the historical 180-product EMA mismatch table.


The [additional 256-entry reading](next256-findings.md) brought that stage to
**1,649 source/ID records**, leaving **4,110 originally saved
abstracts unread**. Sixteen flags retain specific source questions; they are not
sixteen qualified new medicines. Every saved string was read, including an
explicitly documented upstream fragment. `verify_next256.py` replays the exact
selection, exclusions and content hashes against the original retrievals.


[Mebendazole primary qualification](mebendazole-qualification.md) adds a real
other-tissue protein-availability lead and brought that stage of the ledger to **104
appraisals**. Its liver outcomes and missing factorial control are explicit;
no BUBR1 restoration or muscle efficacy is claimed. Further lines remain open.


The [next 512-entry reading](next512-findings.md) brought outer coverage at that stage
to **2,161 source/ID records**, leaving **3,598 originally saved abstracts
unread**. Two readers completed 256 disjoint records each at indices 1217-1728.
It flags 115 source questions; those are not 115 qualified medicines. Four
appraisals follow from its primary work — proTAME, Apcin, fisetin, and dasatinib
plus quercetin — with revisions to the existing Eg5-balancing and navitoclax
entries, bringing that stage to **108 route/compound appraisals**.
`verify_next512.py` replays the exact selection, exclusions and content hashes
against the original retrievals; it passes with the source selection re-derived
from the original corpus. The senotherapeutic reserves are conditional and
combination-specific; neither displaces the niacin/apitegromab comparison.


## Source revalidation after a failed completion audit

The [quality follow-up](review-quality-followup.md) rejects an unreliable 3,598-record completion claim. The initial 186 literal revalidations plus eight 128-record batches reached **1,210 at the batch08 checkpoint**. The [explicit admission](food-admission-findings.md) of 64 already-read food-overlay records reached a historical checkpoint of **1,274 accepted / 3,435 outer / 2,324 pending**. [Batch09](batch09-revalidation-findings.md) then added 128 records at the historical **1,402 accepted / 3,563 outer / 2,196 pending** checkpoint. [Batch10](batch10-revalidation-findings.md) reached the historical **1,530 accepted / 3,691 outer / 2,068 pending** checkpoint. [Batch11](batch11-revalidation-findings.md) reached the historical **1,786 accepted / 3,947 outer / 1,812 pending** checkpoint. [Batch12](batch12-revalidation-findings.md) added 256, reaching the historical **2,042 accepted / 4,203 outer / 1,556 pending** checkpoint. The current original-queue totals are given at the top of this document. This reconciles existing literal reading evidence; it is not 64 newly read papers or independent studies. The preceding 2,161 annotations were not independently re-audited here. Four initial strings are non-abstract items; originally empty fields and recovered bodies retain separate depths. Reading a source does not endorse its claims.

The historical [first](next128-revalidation-findings.md), [second](batch02-revalidation-findings.md), [third](batch03-revalidation-findings.md), [fourth](batch04-revalidation-findings.md), [fifth](batch05-revalidation-findings.md), [sixth](batch06-revalidation-findings.md) and [seventh](batch07-revalidation-findings.md) batches flagged 52, 49, 31, 34, 27, 22 and 14 questions respectively. [Batch08](batch08-revalidation-findings.md) flags **18 further questions**, not 18 drugs. Eight predetermined root cross-reader checks required no substantive changes; no error-rate estimate follows. Saved metadata retain one retracted publication, two erratum-linked records, two preprint records, four linked-preprint notices and five commentary links. Notices were not read merely because their links were retained; no retraction reason is inferred.

[UMK57](umk57-qualification.md) adds actual acute and post-division fidelity
evidence with limited durability, bringing that stage to 109 appraisals.
The subsequent [six-paper qualification](protein-and-tissue-qualification.md)
adds seven compound/tool appraisals and updates metformin, bringing that stage's
ledger to **116**. These are not 116 approved drugs or independent studies. It
retains protein-specific stabilization, ex vivo progenitor support and an
uncontrolled clinical neuromuscular response, while correcting overbroad organ
recovery claims. None establishes MVA efficacy or changes the nominated pair.
The [vascular/enzyme qualification](vascular-enzyme-qualification.md) then adds
six appraisals from three primary studies, bringing that stage to **122**.
Anti-KLK8 antibody earns controlled short-term survival/barrier evidence, K858
earns vascular-lesion evidence, and alglucosidase alfa earns biochemical and
BUB1B-transcript evidence in Pompe myocytes. Drug, genetic and control arms
remain distinct; these findings do not change the nominated pair.
The [STAT3-compound follow-up](stat3-compound-qualification.md) adds SC-43 and
garcinone D, for **124** appraisals at that stage. Its mouse survival, cell protein and
proliferation findings remain specific to the tested models; no MVA efficacy is
established.
A separate
[text-decoder correction](decoder-alias-check.md) recovers 32 catalog-to-record
co-mentions across 11 abstracts under the unchanged alias policy. The current
`literature-v2/` layer propagates them through per-entity coverage and triage;
historical tables are retained. Neither layer is an efficacy classification.
`read_abstracts.py` avoids text-stripping, and `verify_revalidation.py` reconciles
the corrected subset and pending queue without certifying the interpretations.

The [Apcin liver study](apcin-liver-qualification.md) adds actual-drug
pretreatment survival and injury evidence to the existing Apcin entry. Its
survival timing, replication and target-attribution limits remain explicit;
this adds no duplicate compound or abstract-reading count.


The [renal](renal-growth-qualification.md) and
[reprogramming](reprogramming-qualification.md) qualifications add four distinct
compound entries and update existing rapamycin and mirdametinib assessments,
bringing that stage to **128**. Positive renal function and cocktail-level
pluripotency outcomes remain attributed to their actual interventions. These
are not four additional independent efficacy experiments or demonstrated MVA
responses. The 15,357-entity retrieval denominator is unchanged.


The [MAPK injury](mapk-injury-qualification.md) and
[SIRT2 substrate](sirt2-substrate-qualification.md) qualifications add five
compound entries, bringing that stage of the ledger to **133** while preserving all
128 preceding decision objects. FR180204/SP600125 injury protection, the
SB203580 morphology/actin null, NPD11033 protein/ATP findings and SirReal2's
separate ATP result retain their model and endpoint boundaries. All five lack
an exact alias match in the saved clinical catalog; they are research-compound
appraisals, not additions to the 15,357-entity retrieval denominator. No MVA
efficacy or change of nominee is established.


The [Aquamin](aquamin-qualification.md) and [KIF18A](kif18a-compound-qualification.md)
qualifications added five entries, bringing that stage to **138**, with all
133 preceding objects unchanged. Aquamin Soluble is one characterized mixture;
its calcium comparator remains explicit without nominating each constituent.
AM-0277, AM-1882, AM-5308 and AM-9022 retain separate compound-specific results.
All five lack exact aliases in the saved clinical catalog, so the 15,357-entity
retrieval denominator is unchanged. Neither source establishes MVA efficacy or
changes the provisional research pair.


The [corneal](corneal-qualification.md) and [cardiomyocyte](cardiomyocyte-qualification.md)
qualifications add four entries and extend the existing ascorbate-phosphate
assessment, bringing that stage of the ledger to **142**. 2-DG and Stattic retain actual
mouse tissue benefits; rotenone is an adverse cell control. FGF4 plus AA retains
combination-level optical function, while the subsequent hypoxia experiment
measures injury. These are two primary studies, not four independent positive
experiments. Only 2-DG has an exact alias join (CHEMBL1880349); FGF4, Stattic and
rotenone do not match this saved catalog. Neither absence nor a clinical-development
record establishes current approval or efficacy. The 15,357-entity denominator
and provisional research comparison are unchanged.


The [unesbulin qualification](unesbulin-qualification.md) adds one entry, for
**143** appraisals at that stage. It preserves randomized leiomyosarcoma combination
PFS benefit and separate glioblastoma-model survival evidence, with direct
tubulin binding and unestablished normal-tissue selectivity. Its catalog alias
joins CHEMBL4594353; that does not establish approval or MVA benefit. The
[Greatwall qualification](greatwall-qualification.md) extends the existing
okadaic-acid entry: 500 nM restores mitotic-entry timing, without demonstrating
faithful daughter cells. Its older 20 nM completed-mitosis benefit remains
intact. All 109 PR1 IDs remain; 102 decision objects are unchanged and seven have
qualified updates, including CAPE. The reader, its test and original 3,598-record selection
remain unchanged. No source qualification adds abstract-reading credit.

## BET and neural-stem-cell follow-up

**Historical BET/neural checkpoint:** nine appraisals brought the ledger to **152**. All 143
preceding assessment objects were unchanged at that stage. The
[BET review](bet-macrophage-qualification.md) distinguishes seven compounds
from one shared macrophage study; [the sFRP3 review](sfrp3-direction-qualification.md)
adds a recombinant protein and levetiracetam. AAV gene delivery is kept
separate and does not create another compound entry.

| Unelected route | Result retained | Why it does not currently displace the nominees |
| --- | --- | --- |
| ABBV-744 | Controlled mouse plaque reduction during diet challenge. | No matching MVA tissue mechanism, BUBR1 protein or useful-muscle result. Some cytokine comparisons are nonsignificant. |
| Molibresib, (+)-JQ1, apabetalone, GSK046, GSK620, MZ1 | Compound-specific cell inflammatory/lipid findings, including conditional nulls. | Six comparators share one study; no animal benefit may be borrowed from ABBV-744. MZ1 is a degrader. Apabetalone's separate cardiovascular trial missed its primary significance threshold. |
| Recombinant sFRP3 | Reduced excess stem-cell division in APP mice. | Long-term pool/behavioral findings used AAV. Hypomorphic BubR1 models start with deficient division and support the opposite genetic direction. |
| Levetiracetam | Frzb RNA increase among selected seizure responders. | Nonresponders excluded, inconsistent sample count, and no new BUBR1/muscle or all-treated functional result. |

These historical primary follow-ups supplied **zero additional abstract-reading credit**. Current coverage is recorded in the source-revalidation section above. No biological replication is claimed; niacin/restoration and apitegromab/compensation remain provisional under the same functional and safety requirements.

ABBV-744, molibresib, apabetalone and levetiracetam join existing catalog
entities by exact normalized aliases. The other five do not; none of these
joins establishes approval, efficacy or historical formulation equivalence.
The 15,357-entity retrieval denominator is unchanged.

## Unresolved protein-restoration source pointers

A separate bounded protein-restoration search identified two quercetin RNA-source pointers ([PMID15309432](https://pubmed.ncbi.nlm.nih.gov/15309432/), [PMID21632981](https://pubmed.ncbi.nlm.nih.gov/21632981/)). Their complete primary abstracts describe cancer-cell experiments but do not name BUB1B. Lawful publisher and institutional routes did not yield the native expression rows; both BUB1B directions remain unknown. These pointers add no assessed monotherapy entry or reading credit, and quercetin alone does not inherit dasatinib-plus-quercetin outcomes. This limited search is not an exhaustive negative result.


The [food/BCAA/Shoden follow-up](common-preparation-counterevidence.md) adds three bounded preparation appraisals: Camembert physical-function nulls, shared-training rather than added BCAA benefit, and controlled Shoden symptom nulls despite hormone changes. Two mains remain inaccessible; none adds an objective-function positive or changes the nominees.

The [three further preparation comparisons](common-preparation-functional-followup.md) retain shatavari training/contractility and MFGM-mixture mobility positives alongside strength/pain nulls and analysis limits. The protein-meal walking signal remains abstract-bounded, with grip decline in both arms and relative hsCRP increase. These bring the ledger to 303 without changing the nominees.

The [MFGM fraction qualification](mfgm-function-qualification.md) admits one preparation-family appraisal, bringing coverage at that checkpoint to **304**. It retains a small controlled primary strength benefit, agility/conduction positives and comparative strength nulls under shared exercise. Fraction identity, carryover, subsets and unknown participant overlap prevent constituent attribution or four-independent-trial counting; nominations remain unchanged.

The [batch 14 primary follow-ups](batch14-primary-followups.md) add one tannase-treated green-tea preparation and amend existing riboflavin/CoQ10 appraisals. Complete and partial source access remain distinct; narrow exercise signals and nulls do not establish BUB1B/MVA benefit. Ledger coverage at that checkpoint was **305 appraisals**.

The source-specific [Risulose isomaltulose running appraisal](isomaltulose-running-qualification.md) retains a controlled running null alongside metabolic positives. Nonzero symptoms, missing results and printed statistical/unit/timing conflicts remain explicit. It adds one named-preparation appraisal with empty aliases, no eligibility increment and no nominee change; BUB1B rescue and superiority to the current leaders remain unestablished.

The [source evidence and review-depth register](evidence-stage-register/README.md) preserves the complete older 309-appraisal snapshot and an additive link to the current 310th appraisal. Original reading claims, source scopes, preparation/contrast fields and every recorded gap remain visible; byte checks do not establish efficacy or complete primary reading. Validated restoration, function and leader-competitiveness denominators remain UNKNOWN.
