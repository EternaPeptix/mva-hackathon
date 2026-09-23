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

The current **184-appraisal** update adds resveratrol and two distinct topical Potentilla preparations. All 181 previous objects, all 109 PR1 IDs/102 unchanged PR1 objects, the reader/test and original selection remain. Abstract counts are unchanged: **3,066 accepted / 5,227 outer / 532 pending**. The next batch is incomplete and gains no completion credit. The newly found rat resveratrol paper is absent from the frozen corpus; its full-text access remains unresolved.
