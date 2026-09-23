# Reusing PR1 without losing source corrections

**22 September 2026.** The work in [PR1 revision 26b30dd](https://github.com/EternaPeptix/mva-hackathon/commit/26b30ddca828e7fb36e888f0a012bfced4e759e1)
is reusable. Most of it was already present in the evidence branch. This
integration reconciles that revision with the later source corrections; it does
not approve publication or declare the drug search complete.

| PR1 contribution | Disposition |
|---|---|
| Niacin/apitegromab proposal, experimental plan, current pitch and earlier drug comparisons | Already inherited through the common ancestor, `84a1a01`. Preserved. |
| Drug decisions | Before this integration, **107 of PR1's 109 decoded decision objects were equal** to the current branch. The two substantive differences were metformin and UMK57. All 109 decision IDs remain represented; seven later appraisals bring the current total to 116. |
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

The committed coverage remains **2,347 outer source/ID records**, including 186
direct revalidations; **3,412 selected strings require revalidation**. The
earlier 2,161 annotations have not been independently re-audited here. New work
in progress receives no credit until checked and integrated. Further source
questions remain, so the user's conditional publication/submission gate is
still open. The [machine-readable comparison](pr1-reuse.json) records the exact
revision, identical components and source checks.
