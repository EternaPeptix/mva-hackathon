# Public-release boundary review

Reviewed 21 September 2026 against the official event material and the 497 tracked files at private v28 commit `45bb5c47cd1c2f92893ff6f8965051428965bfb3`. **The repository remains private.** This document records content decisions; it does not authorize a visibility change, submission, outreach or deletion.

The [file inventory](PUBLIC-RELEASE-INVENTORY.json) pins every baseline Git blob and assigns an explicit decision: 493 keep, three generalize/update, one publication-license block. It also identifies twenty selected text candidates from the prepared package. Those are proposals for the current private review, not a recursive copy or certification of future edits. Changes made after v28 must receive their own final review.

## Controlling distinctions

Sage's privacy officer distinguishes reusable findings from material carrying the child's genome. The clarification permits HPO terms, limited candidate lists, gene/pathway/mechanism work and code/reports/pitches without prohibited payloads. Genomic files, their slices/indexes, genotype-scale intermediates, pasted variant blocks and raw-genomic trained artifacts remain excluded and subject to deletion. A file's label does not decide its category; its contents and reconstructive scope do. [Organizer clarification](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/2).

The dataset card's CC BY 4.0 metadata does not remove the gated access conditions. [Official dataset page](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data).

The rules require deletion of controlled data within 30 days of close and confirmation to the organizers, prohibit family/MVA Society recontact, and limit identifying disclosures. They require an acknowledgement and the official dataset citation for public communications. The embargo runs from close until the organizing team's summary/preprint and addresses journal manuscript submission; code/models/permitted derived outputs may be shared at any time. Conference abstracts/posters require prior written approval. These provisions do not authorize resharing controlled genomic payloads. [Official rules, data handling, embargo and acknowledgements](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/rules.py).

## Decisions for the actual source archive

“Keep” means eligible content under this review, subject to licensing, correct attribution and a separately authorized future release. It does not mean anonymous, generic, clinically validated or authorized for execution. Hashes and dataset identifiers are interpreted below by function; the organizer did not issue a categorical ruling about each metadata format.

| Category and exact examples | Decision | Grounds and boundary |
|---|---|---|
| HPO inputs: `workflow/configs/proband01-exomiser.yml`, `revel-only-strict.yml`, `revel-only-permissive.yml`, `lirical-subject.json` | **Keep**, explicitly labelled case-derived inputs | HPO retention is expressly permitted. Do not claim these files contain no case-derived information. Their small term sets are not genotype tables; the associated sample label does not make them fictional. |
| HPO syntax in `workflow/scripts/resume_studio.py`, `workflow/tests/test_prioritize.py`, `test_summarize_lirical.py` | **Keep** with method/test scope | The term identifiers are permitted vocabulary/input selection. Tests and the launcher are not evidence of phenotype certainty. |
| Sample/input labels in case configurations and historical identity guards | **Keep** as declared technical identities | The challenge's `PROBAND01` identifier is required by its own Track 1 handler; the dataset page itself shows the WGS filename family. These labels are not a newly discovered person identity. Do not describe the archive as anonymized. New portable examples may use parameters, while exact historical sources and qualification bindings remain preserved. |
| Bounded allele/region selectors: `bub1b-read-phasing.json`, `bub1b-four-lane-phasing.json`, `bub1b-noncoding-review.json`, `joint-haplotype-review.json`, `nabp2-clipped-read-review.json`, `serpina1-cluster-review.json` | **Keep**, bounded research selection | A small selected allele/region hypothesis is a finding or method selector. It must not be mistaken for confirmed phase, diagnosis, a full native record export or a negative call. No rule requires removing every gene or candidate coordinate. |
| Public/reference gene selections: `copy-number-review.json`, `coding-depth-v2.json`, fixed secondary-panel configurations | **Keep** with selection scope | Gene lists and annotation targets do not contain the child's base-by-base genotypes. Retaining a fixed target list improves reproducibility; it does not establish measured coverage. |
| Exact private-input SHA256, counts, expected paths and runtime/source hashes in configs and snapshots | **Keep** as provenance, with nonportable/private-input labels | A digest/pointer does not supply or grant access to its referent. This is our content-based assessment, not an explicit organizer hash exemption. No raw records, reconstructive genotype tables or access credentials may accompany it. Do not execute a historical `ready` flag merely because it is archived. |
| Gene-named documents: `workflow/docs/ADAMTS1-NORMALIZED-OVERLAY.md`, `CLIPPED-READ-REVIEW.md`, `JOINT-HAPLOTYPE-REVIEW.md` | **Keep** | These describe methods and selected technical questions. Gene names are not inherently prohibited; descriptions must retain their stated inference limits. |
| Seven JSON files matching literal GT/AD/PL examples in the inventory | **Keep**, public/fictional provenance retained | They are the four HG002 small-variant/representation receipts plus mitochondrial, ploidy-boundary and structural-CRAM fictional controls. Their public/fictional records must not be relabelled as patient validation. Similar actual-case receipts outside this selection are not covered. |
| `README.md`, `REPRODUCIBILITY.md`, `SOURCE-SNAPSHOT.json` | **Generalize/update**, preserving baseline | Reconcile current versus historical scope, disclose the case-derived selectors/findings, and regenerate the reviewed inventory. This is not a mandate to destroy or silently sanitize qualified historical source. |
| `LICENSE` and any claim that MIT alone states all submission rights | **Block that public-license claim until reconciled** | The event overview explicitly describes submissions, code and reports under CC BY 4.0. State how the event submission grant relates to the existing MIT grant for owned code and separate third-party rights. This review does not require removing MIT, relicensing others' work, or altering the private baseline. [Official overview](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/about.py). |
| FASTQ/VCF/BCF/BAM/CRAM/index content; raw read names/sequences; genotype-scale streams, caches and pasted variant blocks; full raw evidence directories | **Omit** from public selection | A reformat, renamed file or embedding does not turn genomic source data into a report. None is selected by this inventory. Mixed or newly added artifacts are blocked until reviewed by content. |

Local machine/account paths are present in 103 baseline files and are nonportable metadata; no event rule specifically forbids them. They are not secrets or proof of genomic disclosure by themselves. Parameterization can improve a new reusable interface, but silently replacing paths inside frozen sources would invalidate their recorded identities. The inventory records affected lines rather than rewriting them.

This is a static file/category review informed by existing provenance. It is not a comprehensive historical data-flow, provider-contract or secret audit, and a pattern count is not itself proof of privacy. No prohibited current payload was identified by the reviewed categories; that is narrower than proving arbitrary future files safe.

## Prepared findings may be included in this private revision

The prepared package's following members are eligible: `jeweled_track1_report_submission.md`, `jeweled_track2_report_submission.md`, `jeweled_methods_form_answers.md`, `jeweled_evidence-reconciled-v1.csv` and its small provenance file, `jeweled_track2_pitch_reconciled.md`, `PITCH-SONIC-REVIEW-2026-09-20.md`, `PROVENANCE.md`, and the twelve enumerated text appendices. Each exact source hash appears in the inventory. The appendix set includes mechanisms, exposure/assay plans and curated aggregate depth/population/secondary/structural findings; it is not the underlying per-base or per-variant evidence.

Root may select these bounded reports/findings for review and explicitly replace the earlier source-only description. Preserve local links, scientific qualifiers and evidence omissions. The package README should not automatically replace the repository README; media remains a separately managed reviewed artifact. Do not copy the whole artifact tree. Public findings remain subject to the event acknowledgement, attribution, privacy and license conditions.

## Verified delivery facts

The current published deadline is **24 October 2026, 23:59 UTC**, with dates subject to organizer updates. [Official timeline](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/about.py). Recheck before a future timed release; do not substitute the later winners date for submission close.

Submission occurs in the [official Hugging Face Space](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026), through the two submission tabs after account sign-in. Track 1 takes a GitHub URL, PDF/Markdown report and CSV, plus optional display name; its handler accepts `PROBAND01`. [Track 1 form source](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/submit_track1.py). Track 2 takes a GitHub URL, pitch-video URL and PDF/Markdown report, with optional display name and notes; it currently permits three submissions, with the latest reviewed. [Track 2 form source](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/submit_track2.py).

The track instructions permit a private repository during the event and require it public afterward. That is distinct from the manuscript embargo. No sign-in, submission, visibility change, data download, outreach or deletion was performed for this review.

## Final revision delta — 21 September 2026

The inventory's `final_delta` reviews **36 exact payload files**: 21 report files, seven public-replay files, and eight root/index/public-method/status files. All **497 baseline decisions are unchanged**. Later corrections resolve the earlier documentation issues only for the newly pinned bytes.

| Revision content | Decision and boundary |
|---|---|
| Four current report/script texts, revised one-row CSV and small provenance, two report navigation/provenance files | **Keep.** Selected findings, methods, research prediction, proposed experiments and disclosures are permitted output categories. The technical sample label and nominated pair remain explicitly case-derived; no claim of anonymity, raw-record inclusion or clinical validation is made. |
| Twelve dated appendices and current causal-interpretation appendix | **Keep.** Bounded candidate findings, aggregate coverage/mixture summaries, public resource/benchmark results and prospective methods. Per-marker pileups, large interval arrays and native genotype streams remain omitted. Current interpretation is distinguished from dated method detail. |
| Seven public-replay files, reproduction addition and public operational-status receipt | **Keep.** Inputs are public HG002/reference/truth/dependencies. The receipt distinguishes the completed full Mac replay, bounded Linux smoke and final source hardening. Runtime qualification comes from the recorded execution; this privacy assessment does not add it. |
| `LICENSE` | **Keep; baseline block resolved.** The additional event submission grant states CC BY 4.0 for owned submitted material alongside MIT for code, excludes controlled data/third-party resources and acknowledges bounded findings. |
| Root README, CHANGELOG, component index, public benchmark method and saved public review | **Keep.** Current scope, history and public-versus-case evidence are distinguished; raw genomic exclusions remain. |
| `SOURCE-SNAPSHOT.json` | **Regenerated.** Its baseline pin and historical source mappings are preserved; current payload identities are verified in the release fields. The containing Git commit binds the generated snapshot without a recursive self-hash. |

Each payload has an individual path, size, mode, SHA-256, category, decision and grounds. The two audit documents and generated snapshot are bound by the final Git commit without recursive self-hashes. Later additions or edits need a bounded delta review. No source sanitization, private-data deletion, new scientific control suite or public visibility change was performed.

### Attribution and verified project reference

The exact required acknowledgement appears in [report provenance](../../reports/PROVENANCE.md#acknowledgement), linked by the root README. The revised [LICENSE](../../LICENSE) records the event grant without relicensing third-party resources.

A bounded unauthenticated lookup verified **[Rare Disease, Real Kid. The 2026 MVA Hackathon — syn76251147](https://www.synapse.org/Synapse:syn76251147/wiki/642892)** and the [Hugging Face dataset resource](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data). The sole public Synapse wiki confirms 24 October 2026 close and the open-output/controlled-data distinction. It has no dedicated formatted dataset citation or DOI. Cite the verified project reference now and recheck the dataset reference supplied there at publication, as the rules require. Do not invent a DOI or impose one as an additional requirement. These links do not make the underlying case data public.

The revision remains **private and unsubmitted**. Vimeo is on hold at the user's request; competition submission is pending and Cartesia settings are unverified. This content review does not authorize public visibility, waive the organizer terms or turn limited findings into clinical conclusions.

### Scientific-review wording corrections

The later adversarial review corrects the first-anaphase endpoint label, the conditional
compound-heterozygote scoring example and the justification for the EPCR revision.
The same selected findings remain; no new case records, inputs or experimental results
are added. The inventory updates the affected document identities. The CSV and drug
ranking are unchanged; reviewer agreement does not establish probability calibration.
Vimeo work is now on hold at the user's request pending review and possible recalibration.

The pitch document now also includes a proposed narration explaining the existing N1002K
abundance test and drug ranking. This adds no case findings or data category. The approved
narration remains exact, and the movie/captions are unchanged. The text proposal is not a
recording, approval, upload or new scientific result; its bytes and updated navigation are
covered by the inventory's pitch follow-up.

The measured-depth addition uses only the already selected public HG002 BAM and GIAB
regions. It adds aggregate depth sums/counts, explicit filters and a reproduction example
to existing files. No per-base output or additional case finding is included; this does not
establish completion of the separate full calling/comparison replay.

The later terminal record establishes full public Mac replay completion, exact primary
count/audit agreement and measured resource use. Only public aggregate results and
execution identities are added to the existing status/documentation. Final launcher
safeguards and limited Linux checks retain their separate verification scope.

### Earlier Track 2 lead and existing model import

The earlier private revision selected everolimus and added an evidence pipeline and
unelected-route comparison. The new `workflow/track2-model/` contains public Ensembl
responses, reference-based conditional products for the already nominated pair,
assumed cell-population results, source and figures. These are bounded candidate
findings and public-reference/hypothetical computations, not patient transcript
sequences, raw reads, alignments or genome-scale genotypes. The corresponding
inventory delta records this content decision; the repository remains private.

A fresh Mac Studio run of the packaged launcher reconstructed the transcript from
the retained reference and reproduced all scientific population outputs exactly,
excluding execution metadata. This validates replay, not pathway activation or
drug efficacy. The current reports and unrecorded pitch text adopt the new lead;
Track 1 content, the prediction CSV and approved narration remain unchanged.

### Public mouse transcriptome follow-up

`workflow/track2-rna/` adds source and a fresh-download replay of deposited
GSE134781 mouse counts, plus complete aggregate gene/program results, figures and
the executed environment. No case data are added. MSigDB-derived membership and
overlap retain attributed CC BY 4.0 terms. The absence of positive mTORC1 transcriptional
corroboration and scheduling counterevidence appear in the proposal, methods and
unrecorded pitch. The inventory refresh covers this public-data addition; the
repository remains private, and media remains on hold.

### Final route comparison

The current proposal selects everolimus as the lead drug hypothesis and niacin
as the primary comparator, tested separately. Arimoclomol and dasatinib are
conditional reserves. Selection compares causal relevance with actual drug
implementation: prioritizing the former favors niacin, whereas the latter favors
everolimus. This supports the two-route set without asserting efficacy superiority.
The RNA null does not rebut published phosphorylation; state-specific BUBR1 loss
and underreported rapamycin nonbenefit remain concerns to test. The alternatives
table, evidence pipeline, methods, navigation and proposed pitch are synchronized.
No new case data, laboratory response or media is added. Existing computational
outputs, Track 1 science and the approved recording remain unchanged. The
repository remains private; media and competition submission remain on hold.

The subsequent readability and clinical-context revision uses the already included
muscle-atrophy/failure-to-thrive findings, adds an explicit everolimus/arimoclomol
comparison, and summarizes four public RMS clinical sources in the governing
pipeline. This adds no new case payload or treatment recommendation. The original
oncology appendix remains dated method detail; the approved narration is unchanged.

The subsequent decision-rule and export revision specifies the four prospective
drug-by-cohort claims and their joint benefit/preservation requirements. A bounded
exporter reads only the three committed report/methods files and existing prediction
CSV, converting local report links into commit-pinned repository URLs. It does not
include raw data or expand the selected findings. Its outputs remain private review
drafts; file preparation does not authorize submission, media work or public release.

### Track 1 evidence packaging

The current Track 1 revision adds one case-evidence appendix and one selected
aggregate JSON. They recover existing prioritizer/predictor results, support for
the two nominated alleles, BUB1B locus counts, the ten-gene coding comparison,
depth and phase. Public-resource hashes are retained; raw/native variant records,
read names, case-input hashes and genome-scale outputs are excluded. This is the
same limited findings category permitted for reports, with private repository
handling retained. The map makes omitted case prerequisites explicit. GIAB truth
versions are separated. The prediction CSV, classifications, EPCR, Track 2 science
and approved narration are unchanged. No new analysis, public release, media
upload or competition submission is performed.

### Track 1 delivery record

Track 1 model 1 was sent once and has a verified public automated score of
100 rank points and F-max 1.000. The exact combined report upload is retained
alongside the unchanged CSV. This combines already reviewed Track 1 prose and
methods, adds the existing acknowledgement, and updates status wording; no new
case-data category is added. Final durable report-save verification remains
pending after loss of the browser confirmation view. Repository visibility and
Track 2/media status are unchanged. See the current report index for the precise
verification boundary; historical no-submission statements describe earlier states.

### Track 2 research review and new pitch

The 21 September review adds a route-decision appendix and a plain-language main
report, with mechanism-specific reasons to retain, reserve or decline alternatives.
The selected pair remains everolimus and niacin; no biological experiment or
additional case-data category is introduced. The new Skylar pitch text supersedes
the historical Lauren narrative for the proposed submission. New media are prepared
locally for human review; binary media are not included in this source revision.
Cartesia Free-plan status is verified, while its training/retention settings remain
unverified. Track 1 upload bytes and its methods prefix are preserved. Public upload,
Track 2 submission and later host cleanup remain separate pending actions.

### Private repository review package

At the user's request, the repository now includes the unchanged Skylar MP4,
editable PPTX, SRT captions and source storyboard, plus the clean-state AI rescore
and a review landing page. Media contain the same curated report-level findings,
conceptual artwork and synthetic narration; they add no raw case data or new
biological measurement. A separate packaging check found no local-path or
credential-like content in the media/review payloads. These are private review
copies, with human approval, public video publication and Track 2 submission still
pending. The rescore identifies its assessed revision and is not an official score.

### Niacin evidence correction after human review

The main report, route table, evidence pipeline and Track 2 methods now disclose
Pirinen's treated-control strength improvements and distinguish muscle NAD
engagement from disease-specific rescue. The proposed matched treatment controls
and interpretation of genotype-by-treatment effects are explicit. This is a
primary-literature interpretation correction, without new case data, biological
experiments or a change to the drug pair. The prior clean-state review remains
unaltered and is identified as predating the correction. Slide 4 in the deck and
video replaces “five adults” with “five patients,” reflecting the participant aged
17 in Table S1. Track 1, narration, captions and spoken storyboard are unchanged;
private human review and the submission pause continue.

### Carrier-blot and phase corrections after human review

The governing reports, methods and pipeline now describe Sieben Figure 5C as an
unquantified carrier-mouse blot, without total-substrate denominators or inferential
statistics. Confirmed cis retires the proposed compound-heterozygous explanation;
a separate BUB1B rationale requires demonstrated additional trans-allele or dominant
causality. The living Track 1 report and causal appendix make this condition explicit;
the exact submitted upload and CSV remain unchanged. Slides 2 and 3, narration,
caption timing and pitch text are synchronized. No new case data, phase result or
biological experiment is added. Private human review and the submission pause continue.

### Public-entry and pitch comparison follow-up

The 21 September comparison assessed six public Track 2 repositories and four
YouTube pitches using available transcripts and selected frames. A fifth,
age-restricted pitch was assessed from metadata only. This bounded review did
not execute competitors' code, validate their biology or establish a competition
ranking. The JAK/secretory-suppression reserve explicitly attributes the public
proposal that prompted it and separates its cited primary evidence. The current
revision also states the first drug-free model-qualification milestone and its
unsecured resource dependencies, and simplifies judge navigation. No new
scientific or efficacy novelty is claimed; the selected drug pair, biological
results, current media, exact Track 1 upload and prediction CSV are unchanged.

The [current Skylar video](https://vimeo.com/1228990028) is public with English
captions. The user conditionally authorized Track 2 submission after assessment,
and the requested public-entry/pitch comparison is complete. Track 2 portal
submission remains **pending**; this comparison is not a submission-readiness
certification. The repository remains private. These current delivery facts
supersede the historical media-hold and approval-pending statements above.

### Niacin-first revision and bounded reusable demonstrations

The current reports and revised local pitch place niacin first on mechanism
match, with everolimus as a downstream probe. This is a research ordering, not
measured efficacy. The offline [explorer](../track2-model/explorer/README.md)
packages 192 retained aggregate-population comparisons; it adds no simulation,
individual-lineage model, patient forecast or drug-response result. The
[transfer demonstration](../track2-transfer/README.md) combines the existing
bounded BUB1B hypothesis with a synthetic public TRIP13 genotype. Its five cached
records are public reference sequences/annotations and one public ClinVar record.
Sequence consequences and curated route gates do not establish a second patient,
biological replication, automatic variant discovery or generalization accuracy.
The 19 added files are permitted reference, synthetic-result, source or method
content; no additional controlled case payload is included.

The revised MP4, deck, captions and storyboard remain local review artifacts.
The existing [Vimeo video](https://vimeo.com/1228990028) has the superseded drug
ordering and awaits replacement. Track 2 remains unsubmitted and the repository
private. The exact submitted Track 1 report and CSV are unchanged. Current
release identities cover these additions and revisions; historical audit entries
remain records of their earlier states.

### Revision-specific independent review and acknowledgement

The [niacin-first review](../../reports/reviews/track2-niacin-first-20260921.md)
is retained with its assessed artifact hashes and primary-check limits. Its
77/100 is an unofficial, revision-specific judgment, not an efficacy result or
a certified prediction of judging. Only the review text is imported; reviewer-local
reproduction outputs are not included. The standalone Track 2 report now carries
the exact acknowledgement from provenance. Its previously assessed scientific
body is byte-identical, and the other four pinned review artifacts still match.
Navigation links identify this review separately from the earlier assessment.
Public-video replacement and Track 2 submission remain pending; the repository
is private and the exact Track 1 upload and CSV remain unchanged.

The report exporter now accepts bracket-only genotype notation inside code while still rejecting inline/reference-style Markdown links there; this syntax correction changes no scientific, media or CSV content.

### 21 September 2026: revised Vimeo publication verified

At the user's request, [Vimeo 1228990028](https://vimeo.com/1228990028) now
serves V2 of `reports/media/jeweled_track2_skylar.mp4`: 8,699,221 bytes,
180 seconds, 1920×1080, SHA-256
`363827525c61b122141005fbf4440daec2035d4e8d1011d0e528c68055aa2e1f`.
The owner interface reported “Optimization complete”. The saved title is
“Preserving muscle in MVA: niacin first, everolimus as probe | jeweled Track 2”;
the description was updated and English closed captions were replaced with the
current SRT and enabled. Signed-out private-window playback at 00:52 showed the
niacin lead slide and the matching subtitle “Separate experiments connect NAD
and SIRT2 to BUBR1 stability.” This verifies visual playback and caption
correspondence, not a full listening review. The former recording is superseded
publicly at the same URL and remains in private version history.

Active delivery documents and current file identities now reflect this event.
Historical assessments and audit entries retain their original wording; their
artifact identities describe the versions they assessed. Scientific claims,
local MP4/SRT/deck bytes and exact Track 1 uploads are unchanged. Track 2 remains
unsubmitted and the repository remains private.


### 21 September 2026: broad public drug census and new leads

The public drug census adds source-qualified clinical/regulatory metadata,
reproducible retrieval/matching code, screening states and selected primary-source
appraisals. Raw bulk caches and full article abstracts are kept outside Git.
ChEMBL, DrugCentral and UniProt attribution/source-license exceptions are retained.
The additions contain no new controlled genomic payload or patient experiment.

The comparison identifies new lines, including apitegromab as a serious downstream
muscle competitor, AZD1152 same-gene pharmacology with a sarcopenia null, and FRZB
and PAK2 genetic target opportunities. Missing efficacy remains unknown, not an
experimental rejection. Catalog coverage and lexical triage do not certify an
exhaustive worldwide inventory or full-text review of every co-mentioned medicine.
The existing two-drug report/video are marked provisional pending incorporation;
no public repository release or Track 2 submission occurred. Exact Track 1 uploads,
CSV, local media and historical reviewer identities are preserved.


### 21 September 2026: apitegromab integration and Vimeo replacement

The reports, prospective experiment, editable deck and narration now compare
niacin restoration with apitegromab muscle compensation, tested separately.
Everolimus is a conditional mechanism probe. No MVA efficacy or new biological
measurement is claimed. An independent check found no material scientific error
in the revised pitch; historical rubric scores do not score this new proposal.

At the user's explicit request, Vimeo V3 replaced V2 at the existing
[public URL](https://vimeo.com/1228990028). The current movie is 180 seconds,
1920×1080 at 25 fps, with 166.64 seconds of speech and 13.36 seconds of transitions.
All eight slides were visually reviewed; decoded frames, audio and caption text
passed technical integrity checks. Vimeo reported optimization complete; the
new title, description and enabled English captions were saved. Signed-out
playback at 01:19 showed the new apitegromab slide and matching subtitle. No
independent full listening review is claimed. Earlier versions remain in private
history and no longer appear at the public link. Exact Track 1 uploads and
Track 1 methods are unchanged. The repository remains private and Track 2
remains unsubmitted. The current file identities are recorded in the inventory.


### 21 September 2026: recovered literature-screen gaps

Supplemental PubMed metadata and a catalog-independent gene-record inventory
close identifiable retrieval/alias gaps without asserting exhaustive efficacy
review. Seven additional route appraisals and a deeper NAC/ascorbate review retain
model-specific outcomes and unknown actual-allele muscle benefit. Four explicit
preprint/journal links prevent treating added records as independent studies.
New files contain public identifiers/titles, derived screening metadata and
original methods/appraisal prose; abstracts and full texts remain outside Git.
The reviewed fixed input is required before historical manual-review flags may
be reproduced. Independent count/code and selected primary-source checks passed.
The elected proposal, video and exact Track 1 uploads are unchanged.

### 21 September 2026: complete saved-abstract reading and target-assay source checks

The census now retains authored assessments for all 1,432 available saved
abstracts in the 1,460-record gene query. Nine records without abstracts have
accessible-body reviews; 19 remain at title/metadata depth. Source-truncated
thesis text, version duplicates and unresolved full-text questions remain
explicit. Three AI reviewers read disjoint batches, not each paper three times.
The six new route appraisals bring the curated total to 56; geraniin and
TG101209 warrant bounded fidelity/stability qualification, without changing the
niacin/apitegromab nominations.

All 24 retrieved target activity rows are preserved, including experimental
compounds. Two independent readings of the original JAK3-study supplement
confirm four residual-activity values mislabeled inhibition and two swapped
compound/concentration assignments. A separate interpretation layer records
those findings and available nested concentration metadata. Thermal resistance
is not treated as demonstrated usable protein or checkpoint rescue.

The additions contain public metadata, original assessment prose and methods;
raw literature bodies and restricted case files remain outside the package.
Source licenses still apply. New leads and unresolved efficacy preclude an
exhaustive-completion claim; the repository remains private and Track 2
submission remains pending. Current report, movie and exact Track 1 uploads
are unchanged by this literature follow-up.
### 21 September 2026: primary lead qualification and missing-text recovery

The public-source comparison now contains 59 curated route/compound appraisals.
Primary-source checks raise Aurora-B balancing as a functional-compensation
experiment, with arrested-cell attachment endpoints distinct from completed
accurate divisions. Geraniin retains a cellular result but an unresolved intact
parent-exposure bridge; TG101209 remains a thermal probe with competing
pharmacology. KP372-1 and okadaic acid receive contextual/reagent assessments.
Niacin/apitegromab remain the research nominations, not proven superior drugs.

A separate recovery layer covers all 19 formerly title-only access investigations:
nine complete short bodies, one patent abstract, two partial texts and seven
title-only records. Original corpus/review depths remain unchanged. New payloads
contain original assessments and source links, not downloaded abstracts, full
manuscripts, figure images or controlled data. Counts and current report/methods
text were reconciled. Track 1 upload artifacts and pitch media are unchanged.
Repository visibility remains private and Track 2 remains unsubmitted; further
leads and unresolved evidence prevent asserting the user's no-further-leads
condition. These are literature judgments, not drug experiments or clinical advice.


### 21 September 2026: outer-corpus routes and primary muscle endpoints

Sixty selected complete abstracts outside the gene-title/abstract query now have
original assessments, exact selection/identity hashes and discovery-only metadata.
This is 60 of 5,902 outer records, not full outer-corpus or all-drug efficacy
coverage. Primary follow-up elevates TM5484/PAI-1 and hesperetin/CISD2 as
experimental muscle-function alternatives, preserving prevention versus later
rescue, composite motor versus direct force, dose units and genetic/drug
separation. PLK1's patient-drug observation remains explicitly data not shown;
Bletilla identity/mediation and chaetocin model-transfer limits remain visible.
DEVD-CHO is qualified as a degradation-mechanism reagent, not productive rescue.
The ledger now contains 68 route/compound appraisals.

An independent audit reconciled all 60 identities/hashes and the 419-record
metadata selection; it independently read 19 complete abstracts, not all primary
papers. Source reviewers checked the PLK1, hesperetin and PAI-1 integrations.
The final wording avoids unproved cohort independence or genetic-versus-drug
survival superiority. Current census reconciliation and local-link checks cover
the changed package. Historical source provenance and baseline inventory are
preserved; raw abstracts/manuscripts stay outside Git. The nominated pair,
Track 1 uploads and pitch media are unchanged. This remains a private draft;
new leads do not satisfy the conditional publication/Track 2 submission request.

### 21 September 2026: translation checks and another 120 abstracts

The comparison now records 74 curated appraisals. Primary qualification separates
human-stage TM5614 from the muscle compound TM5484, preserves observed adverse
events and the Japanese unapproved designation, and distinguishes chromosome
rescue from parthenogenesis, control-cell harm, drug–siRNA dependence and mixture
effects. Niacin/apitegromab remain nominations without demonstrated MVA efficacy.

Another 120 fixed-order complete abstracts have original assessments, bringing
the two outer selections to 186 distinct records; 5,722 outer records remain
outside those selections. An independent reviewer verified all 120 source
identities/hashes and exclusions, and independently reread 15 abstracts including
all nine follow-ups. Original gene/60-record assessments are preserved. Source
reviewers checked the primary integration; corrections distinguish nonsignificant
from absent effects and explicitly unused from unreported randomization/blinding.

New files contain original judgments and public metadata, not raw papers,
abstracts, figures or controlled case data. Historical source fields and baseline
inventory remain intact. Track 1 upload artifacts and pitch media are unchanged.
This is a private draft update, not publication, Track 2 submission, a new rubric
score or biological validation. The continued leads leave the user's conditional
publication/submission gate unmet.

### 21 September 2026: ranked completion and compound-specific qualification

The additional 253 complete saved abstracts finish the stored 419-record ranking.
Three disjoint selections now contain 433 outer source/ID records, leaving 5,469
outside those selections. An independent reviewer checked every source identity,
raw abstract hash, selection/exclusion and merge, then reread 29 complete
abstracts (a fixed sample plus all 14 new follow-ups). The historical 60/120
ledgers remain unchanged. This is bounded abstract coverage, not 433 independent
studies or a complete worldwide drug-efficacy review.

The ledger contains 79 curated appraisals. New primary interpretation separates
PZ15227 from DT2216, preserves platelet toxicity and the suspended pediatric
combination trial without single-agent causal attribution, and credits actual
navitoclax organ-function and contextual tumor findings. Gedatolisib gains a
conditional cell-size/fidelity rationale. Okadaic acid's partial completed-mitosis
result, SB-431542's duration-dependent growth penalty, EPZ-5676's replication gap
and the nucleoside mixture's BUB1 model are explicit. Independent scientific
checks led to a correction: different nucleoside preparations/endpoints cannot
be treated as a controlled same-mixture comparison.

The three new files contain original assessments, selection metadata and public
source identifiers. Raw primary texts, figures and external review notes remain
outside Git. Prior source fields and baseline release inventory are preserved.
The exact Track 1 uploads, movie, captions, slides and narration are unchanged.
A stale methods sentence was corrected to acknowledge the completed revised
movie. No new rubric score, wet-lab result, public release or Track 2 submission
is claimed. Further primary leads remain; niacin/apitegromab remain hypotheses.

### 21 September 2026: intervention and endpoint qualification

The fourteen follow-ups from the final ranked selection now have targeted
primary-source assessments. One new comparison document resolves named drugs,
genetic interventions, organ-function results and control harm. Eight additional
route/compound appraisals bring the ledger to 87; this includes research tools
and a combination RNA program, not 87 approved medicines. The myeloma summary
is corrected to demonstrated MAD2 vector delivery, without CDC20-restoration
attribution. The existing abstract/selection ledgers remain unchanged.

The new document contains original interpretation and public source links.
Primary bodies, figures and external reviewer notes remain outside Git. Source
identity fields and the baseline release inventory are preserved. Track 1 upload
artifacts and current pitch media remain unchanged. This is a private research
update; it supplies no new biological validation, rubric score, publication or
Track 2 submission. The alternatives do not establish the conditional absence
of further leads.

### 21 September 2026: unfiltered literature and independent catalog crosscheck

A fixed selection adds 384 complete saved abstracts without biological, drug-name
or favorable-result filtering. The four outer selections now contain 817 distinct
source/ID records; 4,942 originally saved abstracts remain unread. All 143 missing
outer records have recovery states. Seventy-one recovered text units were read,
but only 53 contain some assessed scientific content; notices, metadata, unread
long bodies and failed access are explicitly distinguished. An independent audit
replayed source identities, hashes and exclusions and corrected one inferred
drug annotation whose modality was unspecified in the abstract.

Primary and figure-level review qualifies actual nicotinic-acid protection in
toxin-exposed oocytes. The human result is partial and lacks a displayed direct
treatment-contrast significance test. SIRT1 evidence and toxicant-induced BUBR1
recruitment do not establish SIRT2-mediated drug restoration or MVA muscle rescue.
The report, pipeline and existing niacin appraisal retain those boundaries.

The Broad Repurposing Hub crosscheck preserves 7,540 overlapping source entries
separately from the main census. Independent parsing corrects sample-name counts
and three malformed source line breaks. A text-cleaner correction preserves
scientific inequalities and recovers 38 lexical mentions. Source-only replay,
independent matching and compact-table reconciliation passed. Public metadata
and source notices accompany the reproducible script; raw downloads, abstracts,
primary bodies, images and external reviewer work remain outside Git.

Historical source fields, baseline inventory, Track 1 uploads and current media
are preserved. This private update supplies no new efficacy experiment or rubric
score. Further lines remain open; publication and Track 2 submission are pending.

### 21 September 2026: qualification of nine unfiltered-search questions

One source-linked comparison preserves functional tumor/immune, secretion and
osteogenic findings while separating markers, selected survivors, genetic
interventions and culture methods. Six route/compound appraisals bring the
ledger to 93; existing NAC and nicotinamide entries are updated. Collagen
culture support and the already-appraised sFRP3 study add no medicine entity.
Independent source checking and competitor analysis support the provisional
niacin/apitegromab comparison, with no superiority or MVA-efficacy claim.
Historical reading ledgers and catalog denominators are preserved. Raw sources
and external reviews remain outside Git. This private update does not satisfy
the condition for publication or Track 2 submission; further lines remain open.

### 21 September 2026: additional unfiltered review and functional competitors

The next 576 complete saved abstracts bring outer coverage to 1,393 source/ID
records, leaving 4,366 originally saved abstracts unread. Six embedded conference
abstract texts and one primary letter are separately scoped recoveries. Original
assessments preserve source identity, exclusions and reading depth; raw source
text and images remain outside Git.

Three added appraisals bring the ledger to 96. Ursolic acid receives credit for
mouse force and mixed controlled human evidence. MG132 receives distinct
alignment and ex vivo anti-proteolysis credit; neither proves completed faithful
division or contractile benefit. DHMMF phosphorylation accompanies tumor-cell
damage. The existing bortezomib appraisal now credits combination-specific
contracture prevention, with residual mortality, reused controls and reduced
normally innervated muscle area/volume. Independent source and competitor review retain the current research
comparison without declaring clinical superiority. Historical ledgers, baseline
inventory, Track 1 uploads and current media are preserved. Further scientific
questions remain, so publication and Track 2 submission are pending.


### 21 September 2026: intervention identity and direct EMA product coverage

Three new appraisals bring the ledger to 99: distinct NO donors, GSK-923295,
and an access-limited flavopiridol behavioral report. The existing STLC route
now distinguishes partial frog-embryo structural benefit from reduced human-cell
bipolarity and control harm. Source reading depths, nested observation counts,
an assay-label inconsistency and unavailable primary details remain explicit.
No abstract-reading records or experimental efficacy outcomes were added.

The separate EMA supplement retains 1,568 currently authorised human product
rows from an official 2,746-record feed; 180 whole-field names lack exact
cross-references. It does not add 180 medicines or change the 15,357 main/supplement
denominator. Original metadata preserves IDs, status, initial/later dates,
attribution and an unexplained source-clock discrepancy. Raw feed/page responses
stay outside Git; the included code and compact metadata support source replay.
Each product explicitly has unassessed BUB1B efficacy in this layer. Independent
reconstruction confirms source counts, examples and combination traps.

The update preserves historical source fields, baseline inventory, Track 1
uploads and current pitch media. It is a private draft, with no new rubric score,
clinical sign-off, biological experiment or Track 2 submission. Further literature
and product-identity work remains; the publication/submission condition is unmet.


## Further literature and EMA identity review

The latest update adds 256 completed saved-string readings, bringing outer
coverage to 1,649 source/ID records with 4,110 original saved abstracts unread.
One upstream fragment remains explicitly incomplete. Original selections and
reviewer-file hashes are preserved; a proTAME wording refinement distinguishes
SAC-independent metaphase arrest from checkpoint override.

NMN's reported NQO1-model mitotic-completion benefit strengthens its reference
role. Dicoumarol, three independently identified EMA products, and mebendazole
add five appraisals, bringing the ledger to 104. The mebendazole source supplies
protein-availability and liver-pathology evidence but lacks the full knockout
control design needed for its strongest causal claim. The niacin/apitegromab
research comparison remains unchanged.

The 37 recovered EMA aliases/identifiers and explicit typography variants yield
no saved catalog/corpus cross-references. These results neither merge substances
nor modify the historical 180-product mismatch snapshot or 15,357-entity census.
Every product's MVA efficacy remains unknown. Raw papers, source feeds, figures,
abstracts and reviewer working files remain outside Git. Added payloads are
original interpretations, public-source naming metadata and reproducible methods.
Existing source/provenance fields, exact Track 1 uploads and pitch media are
preserved. Independent source checks and source-selection/identity replay passed;
they do not prove clinical efficacy. Repository publication and Track 2 submission
remain pending because further leads and unresolved coverage remain.


## Source-revalidation correction — 22 September 2026

An unreliable 3,598-record completion draft is excluded from this correction.
The added subset has 186 literal-source re-readings, with original source IDs,
hashes, publication metadata and evidence-depth distinctions; 3,412 queue
entries remain pending. A separate decoder overlay recovers 32 catalog-record
co-mentions under the same policy. Original abstract/full text and figure files
remain outside Git. The new UMK57 appraisal and bounded updates to four prior
entries bring the ledger to 109 routes/compounds, not 109 approved medicines.

The selection, source revalidation and decoder replay were checked independently
of efficacy. Existing Track 1 uploads, pitch media, controlled-data exclusions and
source licenses remain intact. The earlier pipeline and inventory snapshots are
preserved as historical provenance. No visibility change, public release or
Track 2 submission follows from this private corrective draft.

## Corrected text replay and six-paper qualification — 22 September 2026

A separate current lexical layer replays all 7,362 original source records under
unchanged alias policy, preserving scientific inequalities. It contains 1,811
specific and 6,195 ambiguous catalog-record matches: 511 of 15,357 entities have
a specific match; 14,846 remain without one. Recovered source roles propagate
into triage while inherited annotations remain labeled. An independent parser
and matcher reproduce both complete match sets. No extra manual-reading credit,
efficacy estimate or worldwide-completeness claim follows.

Six primary papers add seven compound/tool appraisals and update metformin,
leaving 116 entries. Shared papers and sequential treatments are not independent
replications. Positive protein, progenitor-support and reported locomotor results
are preserved, alongside exposure, attribution, model and functional limits.
Niacin/apitegromab remain provisional nominations; no MVA efficacy is established.
The 3,412-string revalidation queue and remaining source questions stay open.

## Further source readings and research reagents, 22 September 2026

The [128-record addition](../track2-drug-census/next128-revalidation-findings.md)
preserves the original 186 accepted revalidations and adds 128 separately read
abstracts in the existing order. Its ledger contains original interpretation,
public bibliographic/version metadata and external-source hashes; raw abstracts
and reading logs remain external. The verifier distinguishes local identity
checks from replay requiring the original retrievals and external evidence.

The [three-paper comparison](../track2-drug-census/vascular-enzyme-qualification.md)
adds six intervention appraisals. It includes research-reagent supply identities,
controlled mouse organ-protection findings, vascular-lesion results and a Pompe
enzyme-replacement transcript response. Drug and genetic interventions remain
separate; no MVA efficacy or comparative clinical probability is inferred.
The resulting 122 appraisals and 2,475 outer records leave 3,284 source strings
pending. No raw literature or new controlled case payload is added. Track 1
upload artifacts and current pitch media are unchanged; this remains a private
draft update without publication or submission.

## PR1 reconciliation, 22 September 2026

The [PR1 reuse comparison](../track2-drug-census/pr1-reuse.md) reconciles revision
`26b30dd` with the later corrected branch at `5a11cfa`. Of PR1's 109 decisions,
107 were already identical; the reader/test files and all selection records
were also already retained. Bounded RPE-1/BJ UMK57 controls and historical-stage
wording are incorporated. The 116-entry ledger, primary qualifications,
corrected text screen, accepted reading counts and pending queue remain intact.
The unaccepted completion annotations and summary remain preserved in PR1's
history, without being admitted as validated coverage in this branch.

The new payload consists of original interpretation and public-source comparison
metadata. It contains no source abstract/full text, figure or controlled case
payload. Track 1 and pitch media are unchanged. This is private branch
integration, not a public release or Track 2 submission.

Added material is original interpretation, public bibliographic/naming metadata,
source hashes and screening code. Raw abstracts, articles and figures remain
external. Historical input files, source provenance, Track 1 uploads and current
pitch media are preserved. Repository visibility and submission status are
unchanged: this is a private draft update.

## Further fixed-order readings and compound qualification, 22 September 2026

The [second 128-record batch](../track2-drug-census/batch02-revalidation-findings.md)
preserves the prior 314 revalidated records and adds 127 saved abstracts plus one
conference abstract. Cumulative revalidation reaches 442 strings, with 2,603
outer source/ID records and 3,156 pending strings. A sample comprised ten
cross-reader checks and six self-rereads; two wording corrections are explicit.
The linked erratum was read separately and changes author names. None of these
checks establishes exhaustive scientific coverage or adds independent studies.

The [compound comparison](../track2-drug-census/stat3-compound-qualification.md)
adds SC-43 and garcinone D, bringing the decision ledger to 124. Their actual
mouse survival and cellular phenotypes are credited within the measured models;
they do not demonstrate BUBR1 restoration or MVA muscle benefit. The updated
[Apcin appraisal](../track2-drug-census/apcin-liver-qualification.md) retains a
mouse pretreatment liver-injury/survival signal, including its timepoint and
mechanism limits, without creating a duplicate compound entry.

New files contain original assessments, public bibliographic metadata, hashes
and verification code. Raw papers, abstracts, figures and reading logs remain
outside the repository. The original release inventory and source-provenance
fields are preserved; current artifact hashes are refreshed. Track 1 uploads
and pitch media remain unchanged. This private draft update authorizes neither
publication nor submission.

## Renal and reprogramming qualification, 22 September 2026

The [third 128-record batch](../track2-drug-census/batch03-revalidation-findings.md)
preserves the preceding 442 and adds 124 saved abstracts plus four conference
abstracts. Accepted revalidations reach 570, with 2,731 outer source/ID records
and 3,028 selected strings pending. Eight independent cross-reader checks found
no substantive corrections; this does not estimate an annotation error rate.
Two linked retraction notices were read separately, with zero extra coverage;
two other linked errata remain open.

The [renal](../track2-drug-census/renal-growth-qualification.md) and
[reprogramming](../track2-drug-census/reprogramming-qualification.md) assessments
retain actual functional positives while separating exposures, interventions,
models and reading depth. Four new compound entries and updates to rapamycin
and mirdametinib bring the ledger to 128 appraisals. All 109 PR1 candidate IDs
remain, with later source corrections preserved. Its historical reuse counts
are explicitly labeled as a checkpoint.

Payloads are original assessments, bibliographic metadata and verification
code; raw literature and reading logs remain external. Targeted independent
checks corrected the delayed calorie-restriction endpoint wording and verified
the reprogramming timing and identity boundaries. The original inventory,
historical source fields, Track 1 uploads and pitch media are preserved. This
private draft update is not a publication or submission approval.
