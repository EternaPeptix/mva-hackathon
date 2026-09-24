# Reusing PR1 without losing source corrections

**22 September 2026.** The work in [PR1 revision 26b30dd](https://github.com/EternaPeptix/mva-hackathon/commit/26b30ddca828e7fb36e888f0a012bfced4e759e1)
is reusable. Most of it was already present in the evidence branch. This
integration reconciles that revision with the later source corrections; it does
not approve publication or declare the drug search complete.

| PR1 contribution | Disposition |
|---|---|
| Niacin/apitegromab proposal, experimental plan, current pitch and earlier drug comparisons | Already inherited through the common ancestor, `84a1a01`. Preserved. |
| Drug decisions | Before this integration, **107 of PR1's 109 decoded decision objects were equal** to the current branch. The two substantive differences were metformin and UMK57. All 109 decision IDs remain represented; seven later appraisals brought the integration checkpoint to 116. |
| Metformin | Keep the later [primary qualification](protein-and-tissue-qualification.md), which distinguishes TET2 stabilization, functional tumor outcomes, exposure and engineered-cell limits. |
| UMK57 | Retain the actual acute and post-division fidelity benefits, durability limits and kinetochore BUBR1 measurement. Add PR1's RPE-1/BJ controls with the narrower interpretation below. |
| Exact-string reader and regression test | Both files are byte-identical between PR1 and the pre-integration branch. No duplicate implementation is needed. |
| 3,598-record selection | Every selection record is identical to the retained [revalidation queue](outer-revalidation-selection.json). The search investment is preserved. |
| Completion annotations and zero-unread claim | Do not admit these as verified coverage. Specific errors remain in the linked revision, and record/hash reconciliation does not establish faithful interpretation. The original draft remains accessible in PR1's immutable history. |
| Earlier 512-record findings | Adopt the wording that identifies 108 appraisals as the count **at that stage**, rather than the current count. |

## What the added UMK57 controls establish

[Orr et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC5130158/) Figure 2B reports
no statistically significant change in lagging-anaphase frequency in
non-transformed, hTERT-immortalized RPE-1 and BJ cells at 100 nM. Both start with
low error rates. This usefully limits transfer from the positive cancer-cell
results. It does not establish inactivity in all non-cancer cells, and neither
line models the child's BUB1B genotype. The patient's attachment-stability defect
has not been measured; calling it the opposite of the cancer-cell defect is
unsupported. The vinblastine comparison supports a different mechanism from
direct tubulin targeting; it is not a direct UMK57 uptake/permeability assay.
The [UMK57 appraisal](umk57-qualification.md) and decision ledger now state these
distinctions. No nomination changes follow from them.

## Why the completion claim is not restored

Direct re-reading of the original, hash-matched strings still contradicts three
missing-name assertions in PR1's `source_abstract_quality` fields:

| Source | Name present in the saved abstract |
|---|---|
| [PMID 29044175](https://pubmed.ncbi.nlm.nih.gov/29044175/) | rhGAA, including the dose-dependent rescue sentence |
| [PMID 39166113](https://pubmed.ncbi.nlm.nih.gov/39166113/) | K858, in the administered-inhibitor result |
| [PMID 24905233](https://pubmed.ncbi.nlm.nih.gov/24905233/) | piceatannol, in the phosphoproteomics method |

These are concrete counterexamples, not an estimated error rate. The source
strings are intact; the cause of all annotation discrepancies has not been
established. The separate mathematical-inequality decoder bug has been
[reproduced and repaired](decoder-alias-check.md), but it is not proven to explain
every manual error. Likewise, the 143 initially empty abstract fields have
separate recovery records and cannot all be called never read.

At that integration checkpoint, coverage was **2,347 outer source/ID records**, including 186
direct revalidations; **3,412 selected strings still required revalidation**. The
earlier 2,161 annotations have not been independently re-audited here. New work
in progress receives no credit until checked and integrated. Further source
questions remain, so the user's conditional publication/submission gate is
still open. The [machine-readable comparison](pr1-reuse.json) records the exact
revision, identical components and source checks.

Subsequent qualified additions are tracked in the [current census](README.md).
The checkpoint counts above are historical, not the current coverage claim.

The later [Apcin qualification](apcin-liver-qualification.md) changes one more
inherited assessment by crediting an actual-drug mouse organ-protection result.
At checkpoint `4975530`, all 109 PR1 IDs remained; 106 decoded decision objects
were unchanged, with metformin, UMK57 and Apcin updated. Subsequent primary
qualifications also update the existing rapamycin and mirdametinib entries,
while preserving all 109 PR1 IDs. The reader/test and all 3,598 selection
records remain identical. The unsupported completion verifier is retained in
Git history rather than used for current coverage.

## Preservation checkpoints, 23 September 2026

At the 170-appraisal checkpoint, the exact requested PR1 revision was an ancestor of this branch. All **109**
assessment IDs remained, with **103 decoded objects unchanged**. Six had later
source qualifications: metformin, UMK57, Apcin, rapamycin, mirdametinib and
okadaic acid. The reader and its test remain byte-identical, and all **3,598**
selection records remain identical. No further merge or duplicate implementation
is needed to use those contributions.

The inherited 143-entry ledger is preserved intact by the BET/sFRP3 follow-up,
which adds nine appraisals for **152**. Reading coverage remains **1,082 literal
revalidations**, **3,243 outer source/ID records** and **2,516 pending strings**.
Those are source/version records, not independent studies. The proposal and
search infrastructure are reusable; the unsupported all-read/exhaustive-negative
conclusion remains excluded. Both PRs remain private drafts and unmerged.


The preceding research-probe follow-up brought the ledger to **156**: four new
appraisals and refinements to GSK620 and ABBV-744, with the other 150 prior
objects unchanged. The PR1 preservation counts above are unchanged. The
separate 249-label SGC inventory and its seven actual abstract readings do not
increase outer-corpus completion counts or prove new efficacy.

The earlier food/epigenetic follow-up brought the ledger to **161**, adding five
appraisals and a new source to (+)-JQ1 while preserving the other 155 preceding
objects. PR1 identity/reader/selection preservation is unchanged. The 3,971 FDA
food source rows are a separate inventory, not unique-drug or efficacy credit.

The preceding follow-up brought the ledger to **170**: taurine and arachidonic acid are added, EPA receives a source refinement, and the other 167 preceding objects are unchanged. [Batch09](batch09-revalidation-findings.md) adds 128 saved-source revalidations after the food 64 reconciliation, for **1,402 accepted /3,563 outer source IDs /2,196 pending strings**. The separate NIH ingredient-group inventory receives no reading or efficacy credit. PR1 preservation remains 109 IDs/103 unchanged objects, the same reader/test and 3,598 selected records.

The preceding **172-appraisal** update preserves all 109 PR1 IDs, with **102 decoded objects unchanged** and seven qualified updates. [CAPE](cape-neurological-qualification.md) is the additional refinement: actual neurological benefits remain, while the BUBR1-associated measurement is correctly identified as RNA, not protein. Two new appraisals cover [Gemini-Cur](gemini-cur-qualification.md) and [10-HDA](10hda-qualification.md); the other 169 preceding objects are unchanged. [Batch10](batch10-revalidation-findings.md) yields **1,530 accepted /3,691 outer records /2,068 pending strings**. The original reader/test, all 3,598 selection records and the exclusion of unsupported completion claims are unchanged.

The preceding **174-appraisal** update adds [Claramine](claramine-qualification.md) and [Nutlin-3a](nutlin-qualification.md), preserving all 172 previous objects and the same 109 PR1 IDs/102 unchanged objects. [Batch11](batch11-revalidation-findings.md) brings totals to **1,786 accepted /3,947 outer /1,812 pending**. Reader/test and all 3,598 selected records remain identical; no exhaustive-negative conclusion is restored.

The preceding **177-appraisal** update adds separate selenomethionine, methylseleninic acid and pectolinarigenin assessments, preserving all 174 preceding objects. All 109 PR1 IDs/102 unchanged objects, the reader/test and all 3,598 original selection objects remain. [Batch12](batch12-revalidation-findings.md) brings totals to **2,042 accepted /4,203 outer /1,556 pending**. No exhaustive-negative conclusion is restored.

The preceding **179-appraisal** update adds capsaicin and the tested glucosamine/chondroitin combination, preserving all 177 previous objects. All 109 PR1 IDs, 102 unchanged objects, the reader/test and all 3,598 selected records remain. [Batch13](batch13-revalidation-findings.md) reaches **2,554 accepted / 4,715 outer / 1,044 pending** after one root sample cohort-wording correction. Neither the new source reading nor measured benefits in other settings restore an exhaustive-negative conclusion.

The preceding **181-appraisal** update adds two preparation-specific Potentilla assessments, preserving all 179 preceding objects. PR1 retention remains 109 IDs/102 unchanged objects, byte-identical reader/test and all 3,598 selected records. [Batch14](batch14-revalidation-findings.md) reaches **3,066 accepted / 5,227 outer / 532 pending** after 16 root spot-checks with no substantive correction. The unsupported exhaustive-negative conclusion remains excluded.

The preceding **184-appraisal** update adds resveratrol and two distinct topical Potentilla preparations. All 181 previous objects, all 109 PR1 IDs/102 unchanged PR1 objects, the reader/test and original selection remain. Abstract counts are unchanged: **3,066 accepted / 5,227 outer / 532 pending**. The next batch is incomplete and gains no completion credit. The newly found rat resveratrol paper is absent from the frozen corpus; its full-text access remains unresolved.

The preceding **194-appraisal** update adds ten assessments and refines MG132 with a separate oocyte source. All 184 prior IDs remain; 183 prior objects are unchanged. All 109 PR1 IDs remain, with 101 unchanged objects and eight qualified refinements. Reader/test and all original selection objects remain identical. [Batch15](batch15-revalidation-findings.md) closes the original queue at **3,598 accepted / 5,759 combined outer / zero selected strings pending**. This does not validate the rejected completion annotations, independently re-audit the earlier 2,161 records, or close primary follow-up and broader discovery.

The preceding **197-appraisal** update adds PD166285, BAY-320 and BAY-524 as distinct experimental comparisons. All 194 prior decision objects remain unchanged, as do all 109 PR1 IDs/101 unchanged PR1 objects, reader/test and original selection. Assay and exposure descriptions are corrected or expanded in source qualifications. Original-queue counts remain 3,598 / 5,759 / zero pending; source-question reconciliation confers no new reading or efficacy credit.

The preceding **200-appraisal** update adds ordinary curcumin, dehydrozingerone and alpha-tocopherol as distinct source-qualified comparisons. All 197 prior objects and all 109 PR1 IDs/101 unchanged PR1 objects remain, along with the reader/test and original selection. Original-queue counts remain 3,598 / 5,759 / zero selected strings pending. Name inventories and source-question updates add no unique reading or efficacy credit.

The preceding **202-appraisal** update adds piperlongumine and diclofenac sodium and refines ordinary curcumin with its original animal evidence. All 200 prior IDs remain, with 199 objects unchanged. All 109 PR1 IDs/101 unchanged PR1 objects, reader/test, catalog/aliases and 3,598 selection records remain preserved. No new abstract-reading credit or exhaustive-negative conclusion follows.

The preceding **203-appraisal** update adds the tested Yiqi Huayu Jiedu preparation and preserves all 202 preceding objects. All 109 PR1 IDs/101 unchanged PR1 objects, reader/test, catalog/aliases and the original selection remain. Database ingredient and gene lists confer no added efficacy or reading credit.

The preceding **206-appraisal** update retains all 203 prior IDs, refines lithium chloride and adds three actually administered Figure 4 comparators. All 109 PR1 IDs/100 unchanged PR1 objects (nine refined, including lithium chloride), reader/test, catalog/aliases and original selection remain. Detection-qualified probe signals and reporter phenotypes add no MVA efficacy or global-completion claim.

The preceding **207-appraisal** update adds all-trans retinoic acid as a bounded culture-barrier comparator, preserving all 206 prior objects. All 109 PR1 IDs/100 unchanged PR1 objects and nine source refinements remain, along with the reader/test, catalog/aliases and original selection. A cross-model gene-set association is not converted into drug-induced BUB1B restoration.

The preceding **212-appraisal** update adds four source-specific CENP-F comparisons and the exact Kalanchoe KT-NB fraction. All 207 previous objects and all 109 PR1 IDs/100 unchanged PR1 objects remain, along with the reader/test, catalog/aliases and original selection. Nine earlier PR1 refinements are unchanged. No new outer-reading credit or exhaustive-negative conclusion follows.

The preceding **216-appraisal** update retains all 212 prior IDs, with 208 objects unchanged and four separate source refinements, and adds four CENP-E/Aurora comparisons. All 109 PR1 IDs remain; 97 objects are unchanged and twelve have qualified updates. GSK-923295, Aurora balancing and Eg5 balancing are the three additional PR1 refinements; MG-132 was already among the prior nine. Reader/test, catalog/aliases and the original selection remain unchanged. No preprint reading or global-completion claim follows.

The preceding **219-appraisal** update adds urolithin A, putrescine dihydrochloride and L-ascorbic acid, preserving every one of the 216 preceding objects. All 109 PR1 IDs remain (97 unchanged and twelve previously qualified refinements). Reader/test, catalog/aliases, original selection and all 47 historical source-question objects remain unchanged. These human-function and preclinical comparators confer no added outer-reading or global-completion claim.

The preceding **222-appraisal** update adds CoQ10, spermidine and mitoquinone, while refining only vitamin K2 among the 219 prior objects. All 109 PR1 IDs remain: 96 unchanged and thirteen qualified refinements, now including vitamin K2’s direct RNA and tissue-immunofluorescence evidence. The reader/test, catalog/aliases, original selection and all 47 historical source questions remain unchanged.

The preceding follow-up retains **222 assessments** and refines resveratrol, icariin and TG101209; the other 219 objects remain unchanged. All 109 PR1 IDs remain: 94 unchanged and fifteen source refinements, now including icariin and TG101209. Reader/test, catalog/aliases, original selection and the 47 historical source questions remain unchanged. These source-specific refinements do not restore an exhaustive-negative conclusion.

The preceding update adds three intervention assessments for **225 total**, preserving all 222 previous decision objects. All 109 PR1 IDs remain: 94 unchanged and fifteen previously qualified refinements. Creatine and glutamine share one trial source; the separately qualified citrulline–metformin combination does not modify existing metformin/AICAR evidence. Reader/test, aliases, original selection and all 47 historical source questions remain unchanged.

The preceding update adds L-leucine and calcium-HMB for **227 assessments**, preserving all 225 previous objects. All 109 PR1 IDs remain: 94 unchanged and fifteen previously qualified refinements. Reader/test, aliases, original selection and all 47 historical source questions remain unchanged. The new source qualifications retain both functional positives and nulls, with no comparative MVA-efficacy or original-queue completion credit.

The preceding update adds one L-carnitine/levocarnitine assessment for **228 total**, preserving all 227 previous objects. PR1 retention remains 109 IDs, 94 unchanged objects and fifteen previously qualified refinements. Reader/test, aliases, original selection and all 47 source questions remain unchanged. Functional endurance is not converted into BUBR1 restoration or comparative MVA efficacy.

The preceding update adds GlyNAC and caffeine for **230 assessments**, preserving all 228 previous objects including NAC prior-28. PR1 retention remains 109 IDs/94 unchanged objects/fifteen earlier refinements. Reader/test, aliases, original selection and all 47 historical source questions remain unchanged. Functional trajectories and selected-biopsy associations do not establish BUBR1 restoration or comparative MVA efficacy.

The preceding update retains all 230 prior IDs, refines only NR and adds three distinct comparisons for **233 assessments**. PR1 retains all 109 IDs, with 93 unchanged objects and sixteen source refinements including NR. Reader/test, original selection and all 47 historical questions remain unchanged. Randomized functional evidence is credited without transferring preparation identity or declaring MVA superiority.

The preceding update retains all 233 prior IDs, refines only NR with NEM6 evidence and adds sodium-AKG and source-labelled epicatechin for **235 assessments**. All 109 PR1 IDs remain (93 unchanged; sixteen previously refined). Reader/test, original selection and all 47 historical questions remain unchanged. The report now presents a compact comparison table; its previous detailed block is preserved in the extended appendix.

The preceding update adds beta-alanine, apigenin and beta-lapachone for **238 assessments**, preserving all 235 preceding objects. All 109 PR1 IDs remain (93 unchanged; sixteen previously refined). Reader/test, original selection, all 47 historical source questions and the preserved comparison appendix remain unchanged.

The current update adds source-labelled betaine and nitrate-containing beetroot juice for **240 assessments**, preserving all 238 preceding objects. PR1 remains 109 IDs/93 unchanged/sixteen earlier refinements. Reader/test, original selection, all 47 historical source questions and the preserved comparison appendix remain unchanged.
