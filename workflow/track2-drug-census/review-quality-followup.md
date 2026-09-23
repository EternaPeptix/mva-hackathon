# Source revalidation and new research leads

**22 September 2026.** The proposed completion of the remaining 3,598-record
outer review is **not accepted**. An audit found annotations describing missing
drug names and results that are present in the saved source. The accompanying
“no further leads” and zero-unread conclusions therefore cannot trigger
publication or submission.

The [corrected subset](outer-revalidation-review.json) contains **186 distinct
saved strings**, re-read directly: a 154-record muscle/force sensitivity screen,
a 21-record purposive quality audit, and 11 integrating-reviewer reads, with
overlap removed, plus six unique queue records from an additional 11-record
decoder-role audit. It retains **34 primary-source follow-up questions**. Of these
186 strings, 182 are abstracts and four are a placeholder, a background fragment,
a correction notice and a conference notice. Reading those four does not make
them usable scientific abstracts. Reviews, preprints and versions are also not
independent drug experiments.

Together with the preserved earlier 2,161 entries, the current ledgers represent
**2,347 outer source/ID records**. **3,412 selected saved strings still require
source revalidation.** This audit did not independently revalidate the earlier
2,161 annotations. The 143 originally empty abstract fields are a separate
category: [later recovery](unfiltered-findings.md) includes recovered abstracts
and bodies, so they must not all be described as unread.

## What failed, and what changed

The rejected selection and annotation hashes reconcile, but that proves identity,
not faithful reading. The original annotation files do not identify their readers
or bind the text snapshot at reading time; the current chunk inputs were modified
after the annotations. The cause of the discrepancies remains unknown. The
original drafts are preserved outside Git rather than silently relabeled as
verified evidence.

| Source | Error in the rejected annotation | Literal-source correction |
|---|---|---|
| [PMID29044175](https://pubmed.ncbi.nlm.nih.gov/29044175/) | Claimed the enzyme name disappeared from the rescue result. | rhGAA is named and dose-dependently reduced glycogen accumulation in Pompe patient-derived myocytes. This is biochemical rescue, not measured force or MVA treatment. |
| [PMID24905233](https://pubmed.ncbi.nlm.nih.gov/24905233/) | Treated a Syk inhibitor as unnamed. | Piceatannol is explicitly named. Recovering its identity does not establish a beneficial intervention. |
| [PMID39166113](https://pubmed.ncbi.nlm.nih.gov/39166113/) | Claimed the KIF11 inhibitor name was absent. | K858 is named; retain its vascular-model effect at that scope. |
| [PMID27077114](https://pubmed.ncbi.nlm.nih.gov/27077114/) | Claimed numerical outcomes and inclusion criteria were missing. | Those values and criteria are present. The paper is prognostic classification, with no administered-drug benefit. |
| [PMID30022161](https://pubmed.ncbi.nlm.nih.gov/30022161/) | Claimed the TET2 phosphorylation residue was unspecified. | Serine 99 is explicit; metformin-mediated TET2 stabilization is a real source claim requiring drug-arm qualification. |

These are purposively checked examples, not an estimated error rate. Corrected
records preserve retraction flags, preprint/version relationships and unread
linked errata. The screening context also incorrectly treated trans phase and
Asn1002Lys instability as established. They remain unmeasured. Neither niacin
nor apitegromab is an MVA-proven benchmark against which challengers must first
demonstrate superiority. Segregation percentages can be functional outcomes;
they should be interpreted alongside absolute output, not dismissed categorically.

## Scientific consequences

The [UMK57 primary appraisal](umk57-qualification.md) adds one entry, bringing
the ledger to **109 route/compound appraisals**. Four existing entries gain
qualified updates. The number includes research tools and combinations, not
109 approved medicines. UMK57 has actual acute and post-division fidelity
benefit in cancer cells, with limited durability and unestablished MVA transfer.

| Further line | Evidence retained | Next question |
|---|---|---|
| CX-4945 | Administered CK2 inhibition raises CDCA3 protein and platinum sensitivity in NSCLC, [PMID34050247](https://pubmed.ncbi.nlm.nih.gov/34050247/). | Quantify substrate specificity, exposure and normal-cell effects before investigating any BUBR1 transfer. |
| Wy-14643 | The intact abstract reports PPARalpha activation, accelerated hepatocyte proliferation and improved post-hepatectomy recovery, [PMID42019507](https://pubmed.ncbi.nlm.nih.gov/42019507/). | Identify recovery endpoints, regimen and controls. Do not transfer this tool-compound result to a clinical PPARalpha agonist. |
| Exogenous bFGF; separately collagenase | Patient-derived marrow-stromal support improved ex vivo, [PMID24014875](https://pubmed.ncbi.nlm.nih.gov/24014875/). | Separate treatment arms and quantify functional support, delivery and context-specific tumor risks. |
| Adrenergic agonists | A report describes beneficial responses in three AGRN-CMS siblings, [PMID28937031](https://pubmed.ncbi.nlm.nih.gov/28937031/). | Identify agents and measured responses. A different neuromuscular lesion is a transfer question, not evidence that the MVA child's weakness must be purely myopathic. |
| Muscle repair and senolysis references | The 154-record sensitivity pass retains 22 questions, including 10 previously unflagged drug, RNA, genetic or review-tracing questions. | Trace actual compounds, age/tissue windows and functional outcomes. Ten questions are not ten new therapies. |

No new primary skeletal-muscle force/grip gain was identified in that 154-record
sensitivity subset. This is a bounded observation, not a claim that none exists
in the unvalidated records or the wider literature. These leads are insufficient
to change the nominees today and sufficient to keep follow-up open.

The [subsequent six-paper primary qualification](protein-and-tissue-qualification.md)
resolves the intervention and endpoint questions for CX-4945, metformin,
Wy-14643, bFGF/collagenase, ephedrine/salbutamol and eplerenone. The first five
papers correspond to five of the 34 original source questions; eplerenone came
from the separate decoder-role audit. UMK57 was qualified earlier. The [status table](primary-question-status.json) links six of the 34 questions
to primary qualifications; the remaining 28 flags are not reconciled as closed
in this correction. Experimental and translational questions also remain for
the qualified sources.

CX-4945 and metformin earn protein-specific pharmacological evidence; bFGF and
collagenase earn ex vivo progenitor-support evidence. The AGRN family report
earns an uncontrolled clinical response after sequential ephedrine/salbutamol,
without assigning the full walking-distance change exclusively to either drug.
Wy-14643's treated arm shows regeneration markers, not demonstrated liver-mass
or functional recovery; eplerenone's cycling effect has unresolved physiological
consequences. The ledger now contains 116 compound/route appraisals. These
qualifications add no new reading credit to the 186-record abstract subset.

## Automated text-screen repair

A separate [decoder comparison](decoder-alias-check.md) reproduced a concrete
bug: arbitrary HTML-tag stripping deletes prose between a mathematical less-than
sign and a later greater-than sign. Holding the saved alias policy fixed,
restricted tag handling recovers **32 catalog-to-record matches across 11
abstracts**, with no losses; the detailed match count changes from 1,779 to
1,811. These remain lexical co-mentions, not efficacy findings. Five direct-gene
flags also reappear. This does not establish the cause of the separate manual
annotation errors.

The original tables remain dated historical outputs. The
[correction overlay](decoder-alias-delta.json) records additional IDs and aliases;
[their roles have been checked](decoder-source-role-review.md): most are predictions,
associations or assay/background uses. Eplerenone has a cellular effect requiring
primary follow-up, now completed at the bounded depth described above.
The [corrected screen](literature-v2/summary.json) now propagates the recovered
matches through all 15,357 entity states and downstream routing. It has 1,811
specific rows, 6,195 ambiguous rows, 511 entities with a specific match and 14,846
without one. The 32 recovered specific rows retain source-checked compound roles;
105 other manual/lexical annotations are clearly inherited. The 96-row keyword
shortlist is not a complete follow-up queue. An independent full replay using a
different parser/matcher reproduced the complete match sets. Historical tables
remain intact. The portable
comparison checks all original records, both decoder implementations and the
same alias policy; it does not change the source cache.

For future manual reading, `read_abstracts.py` emits the exact saved string with
its hash and an end-of-batch marker, with no HTML stripping. Its regression check
preserves inequalities, drug names and markup. `verify_revalidation.py` checks
subset identities, declared depths and the pending count, and can replay the
selection from the original corpus. Passing these checks is not proof of the
scientific interpretation. Raw abstracts and article figures remain outside Git.
