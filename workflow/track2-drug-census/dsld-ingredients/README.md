# NIH supplement ingredient groups: a separate discovery source

**Retrieved 23 September 2026.** The [NIH DSLD API](https://dsldapi.od.nih.gov/)
adds ingredient groups and related label terms beyond the FDA food-additive
inventory. All 27 documented A–Z/Other query responses match their exact reported
totals. This is a source-defined group enumeration, not a list of unique
molecules, current products or demonstrated treatments.

| Preserved or screened item | Count |
|---|---:|
| Source group IDs / distinct exact group names | 6,467 / 6,080 |
| Source name cells / synonym cells | 6,467 / 78,788 |
| Nonempty normalized aliases | 69,707 |
| Aliases shared across multiple groups | 6,228 |
| Saved literature records screened | 7,362 |
| Group–alias–record mention rows | 41,496 |
| Groups / source records with a mention | 1,143 / 6,578 |
| Unflagged mention rows / groups / records | 1,077 / 371 / 774 |
| Exact-name-linked groups: main catalog / status-unverified / FDA / appraisal | 1,121 / 43 / 545 / 65 |
| Added reading or efficacy credit | 0 |

**The automated inventory marks every lexical match context-unreviewed.** A separate [initial role review](first-role-findings.md) now annotates six matching rows across four previously read abstracts, without new unique reading or appraisal credit. Unflagged means only that
none of the specified ambiguity rules applies; it does not establish a specific
intervention, biological benefit or safety. Generic terms dominate many matches.
A no-match group remains unknown. The layer does not add its overlapping group
counts to the existing 15,357-entity catalog.

## What is retained

[groups.tsv](groups.tsv) preserves source IDs, names, categories, source positions,
factsheet references, empty nutrient-information arrays and cross-references.
[alias-cells.tsv](alias-cells.tsv) retains all 85,255 original name/synonym cells,
including the empty synonym. [normalized-aliases.tsv](normalized-aliases.tsv)
keeps collisions visible. Matching uses the corrected text decoder and Unicode
alphanumeric token boundaries with NFKC/case normalization; punctuation and
hyphens separate tokens. All overlapping phrases are found. A component can
therefore match inside a hyphenated derivative, which receives an explicit flag.
No chemical/formulation equivalence follows from a match or a source synonym.

[mentions.tsv](mentions.tsv), [source-records.tsv](source-records.tsv) and
[name-crossrefs.tsv](name-crossrefs.tsv) retain source identity and explicit
unknown states without raw abstract text. [The policy](screen-policy.json)
describes common, short, shared and adjoining-hyphen flags; flags are not
exclusions. Factsheet links include topic searches and receive no reading credit.

All 6,467 nutrient-information lists are empty. Category arrays contain 592
**JSON null elements**, with case variants preserved; they are not the text
"null". Repeated group names span distinct IDs and cannot be silently merged.
[Metadata anomalies](metadata-anomalies.json) record these distinctions.

The API reports version 9.5.0 (May 2026) and index names containing
`2025-09-25n`; the saved static API documentation is version 9.4.0. These strings
are not per-record update dates or a current-market check. The separate queries
are not a transactional database snapshot. [Retrieval metadata](source-retrieval-manifest.json)
bind their exact source bytes, which remain in the external cache.

## Reproduction and checks

From the repo root, using the frozen source cache and corpus:

```sh
python3 workflow/track2-drug-census/build_dsld_inventory.py \
  --source-dir /path/supplement-ingredients-20260923 \
  --catalog-dir workflow/track2-drug-census \
  --corpus /path/literature-gene-unrestricted.json \
  --out /path/dsld-replay --expected-appraisals 222
```

The builder stages a complete immutable run, rechecks inputs and atomically
switches the output's `current` symlink. Resolve it once and check the
[summary hashes](summary.json). Outputs are copied here from the completed run.
The command is a snapshot replay, not a live worldwide inventory update.

The [builder checks](inventory-check.md) passed fifteen data-handling fixtures
and a second identical replay. A different regex matcher agreed on 5,295
alias/source/field hits for ten fixed aliases; it did not independently verify
all 69,707 aliases. A [separate source/table check](independent-source-check.md)
round-tripped every group and source cell and reconciled all name cross-references
at the **168-appraisal checkpoint**. The [parent reconciliation](parent-reconciliation.json)
then verified the 172-appraisal replay: all source fields and six other outputs
remain unchanged. The original 170-appraisal reconciliation is retained. Two
additional name links point to the 10-HDA appraisal, for 6,805 cross-reference
rows and 65 appraisal-linked groups. The earlier independent review remains
scoped to its 168-appraisal checkpoint.

No literature interpretation, new compound appraisal, current availability,
clinical status or MVA efficacy is established by this inventory. The
[decision report](../decision-report.md) records separately qualified evidence.

## Qualified source follow-ups

Primary reviews now distinguish [CAPE](../cape-neurological-qualification.md),
[Gemini-Cur](../gemini-cur-qualification.md) and [10-HDA](../10hda-qualification.md).
These refine or add curated appraisals without new abstract-reading credit.
The first six-row role review remains an immutable account of its earlier
abstract-only scope; subsequent primary evidence is linked separately.

**Source-name inconsistency:** groups 3291 and 6680 both display
“10-Hydroxy-2-Decanoic Acid.” Their matching synonym is respectively
“10-hydroxy-2-decenoic acid” or “10-HDA.” These lexical links do not establish
saturated/unsaturated chemical equivalence or confer the study’s results on
a marketed ingredient. The exact source cells remain unchanged. CAPE likewise
does not confer its findings on the parent caffeic-acid group.

The preceding **174-appraisal replay** preserves all eight output tables byte-for-byte against the 172-entry snapshot: **6,805 name cross-references and 65 appraisal-linked ingredient groups** remain. Both previous parent reconciliations stay historical. Neither new Claramine nor Nutlin-3a appraisal creates an ingredient-name match or extra reading credit.

The preceding **177-appraisal replay** adds five exact name links to selenomethionine across three source groups, for **6,810 cross-references and 68 appraisal-linked groups**. Source fields and six other outputs remain unchanged. Group147 is reported as Selenium but matches through a seleno-L-methionine synonym; this does not equate all selenium forms or products. Repeated Selenomethionine groups2384/7901 are source IDs, not distinct molecules. Historical parent reconciliations and independent source reviews remain preserved; no extra reading or efficacy credit follows.

The preceding **179-appraisal replay** adds two capsaicin name links, for **6,812 cross-references and 70 appraisal-linked groups**. Groups 3951 (Capsaicin) and 3807 (Capsaisin) share that alias; the original spelling and synonym relation are preserved without inferring product or dietary equivalence. The tested glucosamine/chondroitin combination is not equated to either ingredient alone. All source fields and six other outputs remain unchanged, with no new reading or efficacy credit.

The preceding **181-appraisal replay** preserves all eight output tables byte-for-byte: **6,812 cross-references and 70 appraisal-linked groups** remain. Neither tested Potentilla preparation creates a new exact name join. Ingredient or species names are not assigned the mixture/extract findings; no new source, reading or efficacy credit follows.

The preceding **184-appraisal replay** adds two resveratrol name cross-references in one existing ingredient group, for **6,814 cross-references and 71 appraisal-linked groups**. All source fields and six other output tables remain unchanged. The two topical preparations create no new exact join. Name relations do not establish chemical/formulation equivalence or efficacy; the original independent inventory review remains historical.

The preceding **194-appraisal replay** preserves all eight output tables, 6,814 name cross-references and 71 appraisal-linked groups. The ten new source appraisals add no NIH name relation or automatic reading/efficacy credit.

The preceding **197-appraisal replay** preserves all eight output tables, 6,814 name cross-references and 71 appraisal-linked groups. The three experimental-probe assessments create no new NIH name relation or automatic reading/efficacy credit.

The preceding **200-appraisal replay** adds five name cross-references in five existing groups, for **6,819 cross-references and 76 appraisal-linked groups**. Curcumin-name relations include the Turmeric group, and alpha-tocopherol names occur in two vitamin E groups. All five retain shared-alias flags; neither whole turmeric nor another vitamin E preparation inherits the experimental result. All source fields and six other output tables remain unchanged, with no additional reading or efficacy credit.

The preceding **202-appraisal replay** adds one name cross-reference for diclofenac in existing group 8775, for **6,820 cross-references and 77 appraisal-linked groups**. Source fields and six other output tables are unchanged. Ingredient-database presence and exact normalized naming do not establish a supplement preparation, chemical equivalence, eligibility, exposure or efficacy. This layer adds no reading credit.

The preceding **203-appraisal replay** preserves all eight output tables: **6,820 name cross-references and 77 appraisal-linked groups**. The exact Yiqi Huayu Jiedu preparation creates no new name relation, and its predicted constituent list grants no ingredient-level efficacy or reading credit.

The preceding **207-appraisal replay** preserves all eight output tables: **6,820 name cross-references and 77 appraisal-linked groups**. The all-trans retinoic-acid assessment creates no additional NIH name relation. A molecular-name match elsewhere in the catalog does not make retinol, other retinoids or an ingredient preparation inherit its culture result.

The preceding **212-appraisal replay** adds one normalized-name relation from existing Geranylgeraniol group 9833 to the source-labelled trans-geranylgeraniol assessment: **6,821 cross-references, 92 curated relations and 78 appraisal-linked groups**. Only its appraisal-name field and the cross-reference table change; all source fields and six other output tables remain unchanged. The source’s incomplete double-bond geometry and unshown high-dose result remain unresolved. Database synonymy does not establish an equivalent supplement, exposure, current availability or efficacy.

The preceding **216-appraisal replay** preserves all eight output tables, including **6,821 name cross-references, 92 curated relations and 78 appraisal-linked groups**. The CENP-E/Aurora additions and source refinements introduce no new ingredient-name relation or efficacy credit.

The preceding **219-appraisal replay** adds two L-ascorbic-acid synonym relations in existing Vitamin C groups 10 and 3428: **6,823 name cross-references, 94 curated relations and 80 appraisal-linked groups**. Only the two appraisal-name cells and their relations change; all source fields and six other output tables remain unchanged. Shared synonyms do not establish chemical preparation, retail-product or efficacy equivalence.

The current **222-appraisal replay** adds six CoQ10/spermidine name relations across existing groups 7126, 7127, 7134 and 7793: **6,829 name cross-references, 100 curated relations and 84 appraisal-linked groups**. Group 7793 is source-labeled Polyamine (unspecified); its spermidine synonym does not confer class-wide benefit. No mitoquinone-mesylate relation is forced from an unspecified source salt. All source fields and six other output tables remain unchanged.
