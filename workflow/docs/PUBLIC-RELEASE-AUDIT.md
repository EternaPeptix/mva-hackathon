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


## Fourth revalidation batch and injury/SIRT2 qualification

The [fourth 128-record batch](../track2-drug-census/batch04-revalidation-findings.md)
preserves the preceding 570 records and raises accepted revalidation to 698,
with 2,859 outer source/ID records and 2,900 selected strings pending. It contains
124 abstracts, one saved fragment, two commentaries and one preprint abstract.
Eight independent sample checks required no substantive changes. One retraction
notice was read separately, without additional reading credit; four linked
errata remain pending. Original part annotations and earlier ledgers are
preserved. Source/log replay and structural checks are separate from scientific
interpretation and primary-paper qualification.

The [injury](../track2-drug-census/mapk-injury-qualification.md) and
[SIRT2](../track2-drug-census/sirt2-substrate-qualification.md) assessments add
five research-compound entries for 133 appraisals. They preserve transient
injury protection, adverse checkpoint/protein directions and distinct active
comparators. All 128 preceding decision objects and all 109 PR1 candidate IDs
remain. Neither assay protection nor ATP reduction establishes MVA benefit.

Only original assessments, public bibliographic metadata and verification code
are added. Source article text, images, supplementary spreadsheets and reading
logs remain in the external cache. Track 1 uploads and pitch/media are preserved.
This is a private draft update, not permission to publish or submit Track 2.


## Colonoid and KIF18A qualification — 22 September 2026

This private update adds 128 fixed-order source revalidations, bringing accepted
revalidations to 826 and outer ledgers to 2,987; 2,772 selected strings remain
pending. Thirty-two complete literal-reading outputs and eight independent
sample checks support the declared reading depth. Original part annotations and
all previous ledgers are preserved. Structural/full-source replay and altered
fixtures passed; they do not certify scientific interpretation.

Five compound-specific entries bring the ledger to 138, with all 133 preceding
objects unchanged. Aquamin Soluble earns colonoid differentiation evidence;
its calcium comparator and mixture identity remain explicit. Four KIF18A
compounds retain separate cellular/animal tumor results and limited normal-cell
controls. Neither primary source establishes BUBR1 restoration or MVA muscle
benefit. KIF18A original figure images were unavailable; the complete Methods and
selected legends received an independent audit. Aquamin main/supplement images
were inspected, with targeted independent checks of the donor context, mineral
table and Ki67/antibody panels. Raw articles, supplements and images stay outside
Git. These follow-ups add no abstract-reading credit or catalog entities.

Historical source fields, inventory decisions, Track 1 uploads and pitch media
remain intact. No publication, merge, new rubric score or Track 2 submission is
authorized by this package check. The provisional research pair remains
unchanged and unproven in MVA; open source questions and coverage remain visible.


## Sixth literal-source batch and tissue-model qualifications, 22 September 2026

This update retains the prior 826 revalidations and adds 128 exact queued source
strings, for 954 accepted revalidations, 3,115 outer source/ID records and 2,644
pending strings. Both readers' original objects remain unchanged. Eight bounded
independent comparisons required no amendments. Full source/log replay, an
independent repository-only copy and four altered-evidence rejection checks pass;
these are provenance checks, not proof of biological conclusions. Two linked
errata remain unread, and 22 further source questions remain separately flagged.

The corneal and cardiomyocyte qualifications add four compound assessments and
extend the existing ascorbate-phosphate assessment, for 142 total. Corneal tissue
benefits and cardiac optical function remain positive evidence at their tested
scope. BUBR1 rescue, MVA muscle benefit and clinical equivalence are unestablished.
Both appraisals received bounded independent Methods/figure checks. Source text,
images, video files and supplemental data remain external; only original
assessments, bibliographic metadata, provenance hashes and verification code are
included. Track 1, current pitch media and all 109 PR1 decision IDs are preserved.
This private draft update does not approve publication, merger or submission.


## Seventh source batch and distinct tumor/timing outcomes, 22 September 2026

This private update preserves the preceding 954 revalidations and adds 128
fixed-order strings: 1,082 accepted revalidations, 3,243 outer source/ID records,
and 2,516 pending strings. Both readers' original objects are unchanged. Eight
bounded independent annotation comparisons agree; full source/log replay,
fresh-copy structural checks and four altered-evidence rejection fixtures pass.
Fourteen new source questions and nine unread publication links remain visible.
These checks establish identities and declared scope, not biological efficacy.

Unesbulin adds one appraisal, for 143, preserving randomized LMS combination
PFS and separate GBM-model survival benefits. Direct tubulin binding does not
establish selective BMI1 action or MVA rescue. The existing OA entry gains a
Greatwall-depletion mitotic-entry timing result while retaining its prior
completed-mitosis evidence. Both receive bounded independent source checks;
unavailable supplements and unassessed normal-tissue selectivity remain explicit.
Raw papers, abstracts, figures, registry payloads and reading logs remain outside
Git. Original appraisals, metadata and verification code are included.

All 109 PR1 IDs remain; 103 objects are unchanged and six have qualified updates.
Track 1 uploads, current pitch/media, prior source ledgers, historical snapshot
source fields and original inventory decisions are preserved. This update does
not authorize publication, merge or submission. The research pair remains
provisional and unproved in MVA.

## BET/sFRP3 primary follow-up and PR1 reuse, 22 September 2026

Nine original compound/reagent assessments and two primary qualifications bring
the ledger to 152. All 143 preceding objects are unchanged; all 109 PR1 IDs,
103 unchanged PR1 objects, its exact-string reader/test and the 3,598-record
selection are preserved. This update adds no abstract-reading credit and
claims no new independent reviewer. Existing coverage remains 1,082 literal
revalidations, 3,243 outer source/ID records and 2,516 pending strings.

The added material contains original interpretations, public bibliographic
metadata, source hashes and intervention/endpoint distinctions. Raw articles,
figures, supplementary files and external retrieval logs stay outside Git.
No controlled-case payload is added. Track 1 uploads, current pitch/media,
prior source ledgers, historical source identities and original inventory
rows remain intact. Neither compound evidence nor this package check approves
publication, merging or submission.


## Research-probe inventory and primary follow-ups, 23 September 2026

The added SGC layer preserves 391 public source rows as 249 probe labels,
230 parent-scoped control records and seven retrieval hits. Seven complete
abstracts were read in a separate overlay, including a false BSP/PCR match;
accepted outer-corpus totals are unchanged. The source tables are inventory,
not a count of unique drugs or qualified efficacy assessments.

Two primary qualifications credit GSK620's published mouse cancer-survival
benefit and distinguish EZH2 genetic stabilization from enzyme inhibition.
Four appraisals are added for 156 total; GSK620 and ABBV-744 are updated and
150 preceding objects remain identical. All PR1 IDs, its reader/test and
selection survive. New material is original assessment, public source-table
metadata and inventory/replay code. Raw source HTML, abstracts, article text,
figures and supplements stay external. No controlled-case payload is added.
Track 1, current pitch/media, original inventory decisions and historical
source identities are preserved. No publication, merge or submission is authorized.


## Food-substance and epigenetic primary follow-up, 23 September 2026

Five new appraisals bring the ledger to 161; (+)-JQ1 gains a source and the other
155 prior objects are unchanged. Lunasin's transcript signal, MS023's drug-specific
BUB1B reduction, bromosporine's cancer benefit/target ambiguity and GSK591/STM2457's
intervention/statistical boundaries are retained. Delegated source qualification
and a separate lunasin source check are bounded reviews, not whole-submission
clean-state scoring or biological replication. The 36-compound ccRCC screen has
no extra main-ledger appraisal credit.

The new FDA layer preserves all 3,971 retrieved rows, 36 explicit former-use
suffixes and 234 lexical hits. These public source-table facts, original
interpretations and parser/replay code add no controlled-case data. An independent
parser check reproduced the source/match coverage, found failed-rebuild and
ambiguity issues, and prompted fixes. Raw source HTML, abstracts, articles,
figure images, source workbooks and derived previews remain external.

Accepted outer-corpus coverage stays 1,082 literal revalidations and 3,243 source/ID
records, with 2,516 selected strings pending. PR1 reader/test/selection and all 109
IDs remain; 103 PR1 objects are unchanged. Track 1/current media, historical source
identities and original inventory decisions are preserved. No publication,
merge or Track 2 submission is authorized by this update.

## Food abstract and primary follow-up, 23 September 2026

Two primary qualifications bring the ledger to 163, preserving all 161 prior objects. Citral retains cancer-cell effects and the MRC-5 comparator while its BUB1B RNA change is not FDR-significant. Astaxanthin retains composite-culture positives, mixed markers and unresolved assay timing/compound attribution. Separate source checks corrected the timing description and preserved source limitations; this is not a whole-submission score or biological replication.

The food overlay records 67 saved abstracts behind 73 literal-policy rows, with raw text/logs external and no automatic outer credit. It preserves retraction, preprint/version and intervention-identity distinctions; 161 common-name rows remain unreviewed. Separate identity/data validation reconciles the exact source strings and overlaps. FDA source rows and 234 mention rows are preserved, with only astaxanthin/citral name cross-references added. The earlier independent parser review remains identified as its 161-entry snapshot.

Batch08 adds 128 disjoint literal source readings for 1,210 accepted revalidations and 3,371 outer source/ID records, leaving 2,388 selected strings pending. It retains one retracted publication, two erratum-linked records and two preprints; linked notices are not claimed as read. Eight root cross-reader comparisons required no substantive change. Source/log replay verifies identities and declarations, not scientific truth.

Original assessments, compact bibliographic/status metadata and verification code are the only new repository content. Raw abstracts, papers, figures, workbooks and previews remain external. No controlled-case payload is added. All 109 PR1 IDs and 103 unchanged objects, Track 1/current media, original inventory decisions and historical source identities are preserved. No publication, merge or Track 2 submission is authorized.

## Common-name food review and lipid comparison, 23 September 2026

Five new compound/role appraisals bring the ledger to 168, preserving all 163 preceding objects. BMS-309403 gains an actual-compound BUB1B transcript observation from twelve public normalized-expression values. This is a descriptive extraction, with source/timing/target limits and no inferred P value, protein or functional rescue. Oleic acid and free EPA retain distinct lipid-comparator roles; PA-d31 and bromate are stressor/adverse entries. The unavailable fatty-acid comment remains title/metadata-only. A duplicate Wy-14643 flag is linked to its existing marker-only qualification.

All 234 food query rows across 213 source IDs now have role annotations; 212 saved abstracts were read and one is missing. The common-name overlay retains its original annotations and publication/version flags. A separate reconciliation admits 64 already-read source records once, preserving all 1,210 prior acceptances and all 3,598 original selection objects: 1,274 accepted revalidations, 3,435 outer source IDs and 2,324 pending strings. Twelve predetermined root spot-checks required no substantive correction. Exact source/log replay and negative checks verify identities and declared gates, not scientific truth or independent reading of every record.

New payloads are original assessments, compact metadata, replay code, twelve public expression values/sample mappings and an original plot of those values. Raw abstracts, complete matrices, source articles and third-party figures remain external. No controlled-case payload is added. The FDA replay adds two exact name cross-references with source values/all 234 mentions preserved; the earlier independent parser audit stays historical. PR1 identities and source selection, Track 1, current pitch/media, original inventory decisions and historical source identities are preserved. This update does not authorize public release, merging or submission.

## Supplement ingredient discovery and metabolic qualifications, 23 September 2026

The ledger reaches 170: taurine and arachidonic acid are added; EPA gains a separate oocyte-source qualification, and 167 preceding objects remain unchanged. Actual taurine muscle electrophysiology is separated from transporter-knockout survival; source access limits remain explicit. ARA dose-specific culture injury and genetic metabolic rescues do not become nutrient-treatment benefit. The main report now presents measured reserves before a concise source/validation summary, with full audit detail linked.

The separate NIH DSLD layer preserves 6,467 source groups and 85,255 original name/synonym cells, including JSON null categories, repeated names and an empty synonym. All 27 saved partition totals reconcile. The 41,496 name matches are discovery records, not unique drugs or efficacy. Independent source/table checks apply to the 168-appraisal checkpoint; a separate parent replay against 170 preserves source fields and adds three name links. Six matching rows receive an initial role review across four previously read abstracts, without new unique, outer, primary or appraisal credit.

Batch09 adds 128 disjoint saved-string readings after the 64-food admission: 1,402 accepted revalidations, 3,563 outer source/ID records, 2,196 pending. Its retraction, preprint, erratum and lineage metadata remain distinct from biological endorsement. Eight predetermined root comparisons required no substantive change. Initially empty fields retain their separate recovery depths. Full source/log replay checks declarations and identity, not scientific truth.

New payloads are original assessments, public agency ingredient facts, compact bibliographic/status metadata and code. Raw abstracts, complete articles, third-party figures and controlled-case payloads remain outside the repo. Original inventory decisions, historical source identities, Track 1 and current media are preserved. No public release, merge or submission is authorized by this update.

## Common-compound assay qualification, 23 September 2026

The private evidence branch now has 172 appraisals: Gemini-Cur and 10-HDA are new; CAPE retains neurological-function positives with its BUB1B result corrected to RNA, not measured protein loss. All other 169 preceding assessments are unchanged. Source figures, tables, controls and reporting contradictions are distinguished in the linked qualifications; none establishes MVA efficacy or ingredient interchangeability.

Batch10 adds 128 exact saved-source revalidations after excluding all 1,402 prior acceptances. Coverage at that checkpoint was 1,530 accepted / 3,691 outer source-ID records / 2,068 pending strings. The root read and compared eight fixed samples; no full independent 128-record interpretation review or error-rate estimate is claimed. Prior reading artifacts and the 3,598-record queue are preserved.

The NIH replay adds only two lexical links for 10-HDA, with the decanoic/decenoic source-name inconsistency explicit. It adds no source groups, reading credit or chemical-equivalence claim. FDA and SGC data tables are unchanged; input identities are refreshed. The original independent inventory checks remain historical.

Ten new files contain original qualifications, compact selection/review metadata and verification code. No article body, raw abstract, third-party figure, full supplementary matrix or controlled case input is added. Track 1, current pitch/media, baseline release decisions and historical source identities are preserved. This is not public-release or Track 2 submission approval.

## Metabolic and tumor-associated function, 23 September 2026

The ledger now contains 174 appraisals. Claramine gains a distinct metabolic-function reserve based on receptor-processing/signaling experiments and a separate diabetic-mouse glycemic result. Nutlin-3a gains a tumor-associated physical-function reserve: movement and baseline-relative grip benefits remain credited without inferring normal-muscle restoration. The source-specific comparator, regimen and access limitations remain explicit; neither demonstrates BUB1B/MVA rescue. All 172 preceding appraisal objects are unchanged.

Batch11 adds 256 exact saved strings in 64 outputs, with 16 predetermined root source/annotation checks. Totals at that checkpoint were 1,786 accepted revalidations, 3,947 outer records and 1,812 pending. One retracted review is excluded as supporting evidence, and version/correction links do not themselves confer notice-reading credit. The full source/log replay checks identities and declarations, not scientific truth or independent review of all records.

The 174-appraisal inventory replay preserves every NIH, FDA and SGC output table. Input summaries change because the decision ledger changed; no new ingredient-name match or automatic efficacy/reading credit follows. Prior independent source reviews and parent reconciliations retain their historical scope.

Eight new files contain original qualifications, compact review/selection metadata and replay code. Raw abstracts, source papers, third-party figures and controlled case payloads remain external. Track 1, current media, original inventory decisions and historical source identities are preserved. Both PRs remain private drafts; no merge or Track 2 submission is authorized by this update.

## Selenium forms and flavone qualification, 23 September 2026

Three distinct entries bring the ledger to 177, preserving all 174 preceding appraisal objects. Selenomethionine and methylseleninic acid retain their relative transcriptional comparison and separate source-limited antimutagenesis lead. Neither relative RNA abundance nor phase redistribution establishes BUBR1 restoration. Conflicting GEO reference/time labels prevent an inferred rescue ratio; the adult SELECT prevention result retains its distinct intervention and population. Pectolinarigenin retains actual cancer-cell and xenograft positives, while the low-mass BUB1B-assigned spot does not establish total BUBR1 reduction. None demonstrates MVA efficacy or equivalence to dietary mixtures.

Batch12 adds 256 exact strings in 64 outputs. Sixteen predetermined root source/annotation checks required no substantive correction. Coverage at that checkpoint was 2,042 accepted revalidations, 4,203 outer records and 1,556 pending strings. Fifteen follow-up questions, three preprints and six correction-linked records retain their source roles; no saved retraction label occurs in this batch. Full source/log replay and rejection checks test identity and declarations, not scientific truth or independent reading of every record.

The NIH replay adds five selenomethionine name links across three source groups, with source fields and six other outputs unchanged. Group147 is generically named Selenium but matches a specific synonym; this is not equivalence across selenium forms. Historical parent reconciliations remain intact. FDA and SGC output tables are unchanged, and inventory work adds no reading or efficacy credit.

Eight new files contain original qualifications, compact source metadata and replay code. Raw articles, figures, complete GEO matrix/platform data and controlled case inputs remain external. A bounded independent selenium source check and selected root pectolinarigenin figure checks retain their actual scope. Track 1, current media/model/transfer, historical source identities and original inventory decisions are preserved. This private draft update is not approval to merge or submit Track 2.

## Common symptom and biomarker qualification, 23 September 2026

Two new appraisals bring the ledger to 179, preserving all 177 previous objects. Capsaicin retains local analgesic behavior with sensory denervation; the BUB1B transcript is injury-associated, not an established drug-response measurement. The tested glucosamine/chondroitin combination retains randomized human CRP reduction, with CRP confirmed as the registered primary endpoint. Its exploratory plasma antibody data do not establish intracellular BUBR1 restoration. Raw workbooks, articles and figure images remain external; selected source and arithmetic checks do not confer full independent scientific review.

Batch13 adds 512 exact saved strings in 128 outputs. Sixteen predetermined root source/annotation comparisons required one cohort-wording correction. Totals reach 2,554 accepted revalidations, 4,715 outer source/ID records and 1,044 pending strings. Retraction, expression-of-concern, correction and version warnings remain attached to their sources; linked notices do not gain reading credit from metadata. Full source/log checks verify identities and declarations, not efficacy or independent review of every record.

The NIH replay adds two capsaicin name links across source groups3951 and3807, retaining their Capsaicin/Capsaisin spellings and the chemical/product-equivalence limitation. All source fields and six other outputs are preserved. FDA and SGC data tables remain byte-identical; historical independent inventory checks and parent reconciliations retain their original scope.

Eight new files contain original assessments, compact selection/review metadata and verification code. Track 1, current media/model/transfer/RNA, all prior assessment objects, historical source identities and original inventory decisions remain preserved. This private draft update neither changes repository visibility nor authorizes merging or Track 2 submission.

## Potentilla preparation-specific protection, 23 September 2026

Two preparation-specific entries bring the ledger to 181, preserving all 179 preceding objects. The methanolic root extract retains a 2/10 versus 10/10 carcinoma observation at 260 days and partial premature chromosome-separation protection during carcinogen co-exposure. ECGU retains the separation phenotype, without inheriting the extract tumor endpoint or individual-component efficacy. Native figures, tables and supplement qualify reporting discrepancies and nonsignificant marker comparisons; BUBR1 and muscle function were not assayed. Root and bounded independent checks retain their actual source and arithmetic scope, without biological replication.

Batch14 adds 512 exact strings in 128 outputs, preserving all 2,554 prior acceptances. Sixteen root source/annotation checks required no substantive correction. Totals reach 3,066 accepted revalidations, 5,227 outer source/ID records and 532 pending strings. Retraction, correction and preprint relations remain source-status evidence, not notice-reading or independent-study credit. The new verifier splits JSONL on literal LF to preserve embedded Unicode paragraph separators; source/log bytes and earlier verifiers remain unchanged.

Fresh 181-appraisal name replays preserve every FDA, SGC and NIH output table. There are still 6,812 NIH name cross-references and 70 appraisal-linked groups; ingredient/species names are not assigned preparation-specific efficacy. Input summaries and parent reconciliations retain their historical/current distinction, with no automatic reading or efficacy credit.

Six new files contain original interpretation, compact source-selection/review metadata and verifier code. Raw articles, images, native supplementary material and raw-reading logs remain external. All original inventory decisions and historical source identities are retained. Track 1 and current media/model/transfer/RNA artifacts are preserved. This private draft update is not approval to merge, change visibility or submit Track 2.

## Resveratrol and topical preparation qualification, 23 September 2026

Three additions bring the ledger to 184, preserving all 181 preceding objects. Resveratrol retains human exercise-associated power and torque benefits with actual between-group interactions, smaller-pilot uncertainty and conflicting human evidence. Its rat BUBR1-staining source remains abstract-only. Two successive topical Potentilla extracts retain wound-repair results without transfer from the earlier oral preparations. No BUBR1-MVA efficacy, independent biological replication or clinical superiority is established.

An independent reader completed the wound source and a bounded second reading of the 2017 resveratrol main text. Account limits interrupted the final resveratrol write-up and parallel abstract work. Root completed its primary review and package checks; no completed independent package review is claimed for this update. The incomplete Batch15 remains external and gains no completion credit. Accepted counts remain 3,066 revalidations, 5,227 outer source/ID records and 532 pending strings. The rat resveratrol source is absent from the frozen corpus, demonstrating a discovery gap beyond queue completion.

Replayed FDA and SGC output tables remain unchanged. NIH gains two resveratrol name cross-references in one existing group, for 6,814 cross-references and 71 appraisal-linked groups. Every source field and six other NIH output tables remain unchanged. Name relations are not chemical/formulation equivalence or efficacy. Historical independent inventory reviews retain their original scope.

Four new files contain original source qualifications and compact provenance. Raw articles, images, registry data and unfinished reading logs remain external. Track 1, current pitch/media/model/transfer/RNA and protected experimental-plan bytes are preserved. This private draft update does not merge either PR, change visibility or submit Track 2.

## Original queue completion and drug-specific phenotype checks, 23 September 2026

The original 3,598-record queue is revalidated after the final 532 literal reads and 16 root source/annotation spot-checks. Combined outer coverage is 5,759 source/ID records, with zero pending in that selection. The earlier 2,161 annotations were not independently re-audited; primary questions, absent-source recovery depths and wider discovery remain distinct. This does not reinstate rejected annotations or an exhaustive-negative conclusion.

Ten new compound/combination assessments and a separate MG132 source refinement bring the ledger to 194. All 184 prior IDs remain and 183 objects are unchanged; all 109 PR1 IDs remain, with 101 unchanged. Native drug/normal-control oocyte panels retain partial improvements and stalled normal maturation. Original bladder-cancer tables distinguish large-tumor prevention from all-lesion clearance, and targeted original workbooks show lower tumor Bub1b RNA after erlotinib/gefitinib. No result establishes MVA rescue. Independent primary checks are bounded and documented; root alone inspected the workbook cells. The prior resveratrol independent write-up is now complete without a substantive change to its scientific claims.

Eight new files contain compact original interpretations, batch annotations, selection and verifier; raw article bodies, figures, workbooks and reading logs remain external. Cached FDA, SGC and NIH replays preserve all output tables. Track 1, methods prefix, pitch/media/model/transfer/RNA and protected experimental-plan bytes are unchanged. Both PRs remain private drafts, unmerged, and Track 2 unsubmitted.

## Progression, Bub1 inhibition and source-question reconciliation, 23 September 2026

The ledger reaches 197 with separate PD166285, BAY-320 and BAY-524 assessments; all 194 preceding objects and all original queue payloads remain unchanged. PD166285 produces progression effects with residual spindle defects. Formal paclitaxel synergy belongs to BAY-320, while severe errors in the RPE1 combination at 4 nM paclitaxel limit normal-cell-sparing claims. The earlier BAY biochemical assay conflation is corrected; BubR1 kinetochore localization is separated from protein restoration. Resveratrol and bladder-source exposure details are added from already-read primary material.

Five of the final batch's 47 historical questions have material direct qualified answers; 42 have none in the inspected package. Two study-level questions are answered within source limits and three retain source ambiguities. Matching names or related experiments do not close a question. No new abstract-queue or biological-replication credit follows. Original counts stay 3,598 accepted revalidations, 5,759 combined outer source/ID records and zero pending in that selection; wider discovery and primary questions remain open.

Six new files contain bounded original source qualifications and component-level reconciliation, with portable document links. Raw article bodies, figures, workbooks and reading logs remain external. FDA, SGC and NIH replay outputs are unchanged. Current PR1 preservation counters use the latest baseline; historical checkpoint objects remain preserved. Track 1, methods prefix and protected media/model/transfer/RNA/experimental plan bytes remain unchanged. This private draft update does not merge either PR, publish the repository or submit Track 2.

## Ordinary curcumin, structural comparator and alpha-tocopherol, 23 September 2026

Three additions bring the ledger to 200, preserving all 197 prior objects. Curcumin has measured cancer-cell arrest/death and Cdc27/APC/C-associated effects without definitive direct-target causality or normal functional rescue. Dehydrozingerone is a distinct tested weak-response comparator. Alpha-tocopherol has partial pretreatment/co-exposure MTT protection at both 100 and 200 micromolar; lower BUB1B RNA is not protein restoration. Source identities, assay controls, denominator limits and unmeasured endpoints remain explicit. Full-resolution source checks corrected an integrating-reader digit misread before any source-error flag was admitted.

Seven of the final batch's 47 questions now have material direct source-qualified answers: four original study-level questions answered within source limits and three retaining source ambiguities. Forty lack direct qualified answers in the package. Only the two newly addressed record objects changed; original questions, flags and frozen readings remain unchanged. Queue counts stay 3,598 accepted / 5,759 combined outer / zero original selection pending; broader research remains incomplete.

Four new files contain original bounded source qualifications and compact provenance. Raw scientific bodies, figures and reading logs remain external. FDA and SGC outputs are unchanged. NIH adds five explicitly ambiguous name links across existing groups, without changing source fields or six other tables. Alpha-tocopherol's main-catalog name join remains ambiguous between TOCOFERSOLAN and VITAMIN E, not a verified experimental-formulation match. Source appraisals confer no automatic chemical equivalence, product efficacy or new abstract credit.

Track 1, methods prefix, current pitch/media/model/transfer/RNA and protected experimental plan remain unchanged. The repository and both draft PRs remain private/unmerged, with Track 2 unsubmitted.

## Curcumin animal evidence, piperlongumine and diclofenac, 23 September 2026

The ledger reaches 202 assessments: two new compound entries and one refinement of ordinary curcumin, retaining all 200 preceding IDs and 199 unchanged objects. Curcumin gains reported mouse tumor and disease-endpoint lifespan benefits with an unresolved 1 mg/kg versus 1 g/kg dose conflict. Piperlongumine has complementary cellular anticancer outcomes but no direct/selective CDK1 or BUBR1 rescue evidence. Diclofenac retains low-micromolar pancreatic-cell combination signals without formal synergy, uniform superiority over single agents, or a measured normal-cell exposure bridge. Its author-reported midpoint has undefined normalization; 2 mM NAC is explicit only for the morphology panels.

The prior 2012 curcumin qualification is narrowly corrected: its linked tumor blots do not establish exact specimen reuse or an independent efficacy cohort. All other prior qualifications remain unchanged. Six new files contain original source interpretation and compact provenance; raw articles, figures and reading logs remain external. Root read all three scientific mains and selected original figures; separate independent checks retain their actual scope.

Only final-batch follow-up records n35, n52 and n225 change. Nine of 47 questions now have direct source-qualified information: five study-level questions answered within limits, four retaining unresolved source components, and 38 without a direct qualified answer. Original questions, flags and queue readings remain unchanged at 3,598 accepted / 5,759 combined outer / zero originally selected strings pending. Broader research remains incomplete.

FDA and SGC output tables remain unchanged. NIH gains one exact-name relation in an existing diclofenac group, with no chemical/formulation equivalence or reading credit. Track 1, its methods prefix, pitch/media/model/transfer/RNA and protected experiment plan remain unchanged. Repository visibility, draft PR state and Track 2 submission state are unchanged.

## Exact-mixture lung-injury qualification, 23 September 2026

The ledger reaches 203 assessments by adding the tested Yiqi Huayu Jiedu eight-herb decoction; all 202 earlier objects are unchanged. The actual animal regimen is seven-day pretreatment before LPS, with qualitative lung-histology and author-reported AKT1/pan-Ras findings. Respiratory function, survival and post-injury recovery are unmeasured. BUB1B is only a database disease-target input. Predicted constituents are not measured batch chemistry or isolated-ingredient efficacy. Representative blots cannot reconstruct the five-unit means; no demonstrated image/summary contradiction is inferred.

Two new files contain original interpretation and compact provenance; raw sources remain external. Root read the scientific main, all main tables and selected figures/pages; source and independent reviewers retain their separate, explicitly bounded supplementary scopes. Only final-batch record n156 changes. All original questions and other 46 objects remain. Ten of 47 questions have direct qualified information: six study-level answers within limits, four with unresolved source components and 37 without direct qualified answers. Outer-reading totals remain 3,598 accepted / 5,759 combined outer / zero originally selected strings pending.

A separate targeted search leaves two quercetin native-RNA pointers unresolved; no assessed monotherapy object, protein-restoration direction or new reading credit is added. All FDA, SGC and NIH output tables remain unchanged. Track 1, methods prefix, protected pitch/media/model/transfer/RNA and experimental plan remain unchanged. Private/draft/unsubmitted status and provisional nominees are unchanged.

## Lithium and administered astrocyte comparators, 23 September 2026

The ledger reaches 206 assessments: all 203 prior IDs remain, 202 prior objects are unchanged, LiCl is refined and three actually administered Figure 4 comparators are added. The new source supplies reporter-count/morphology positives, not neuronal function or BUBR1 restoration. BAPN lacks a reported concentration, pioglitazone retains the printed 1 mM exposure and the methyl-epicatechin preparation has conflicting names; automatic matching of the latter is withheld. Other computational nominees do not become administered agents.

Four new portable files contain source interpretation, a compact three-probe/eight-array summary and a standard-library replay with exact source hashes and native detection calls. One probe is absent throughout. Lower processed signals are not a validated differential-RNA, protein or MVA harm result. Raw papers, whole GEO datasets and figure images remain external. Root and independent reading scopes remain separate; no raw-array reprocessing or new outer-reading credit is claimed.

Only final-batch n390 changes; all original questions and other 46 objects remain. Eleven of 47 questions now have direct qualified information: six study-level answers within limits, five with unresolved components and 36 without direct answers. Original-queue totals remain 3,598 accepted / 5,759 combined outer / zero selected strings pending. The Track 2 abstract coverage is refreshed within 500 words. All FDA, SGC and NIH output tables remain unchanged. Track 1, its methods prefix, protected pitch/media/model/transfer/RNA and experimental plan remain unchanged. Private/draft/unsubmitted status and provisional nominees are unchanged.

## All-trans retinoic acid and culture-barrier function, 23 September 2026

The ledger reaches 207 assessments by adding one all-trans-retinoic-acid culture-barrier comparator, with every preceding 206 assessment object unchanged. A direct RA/no-RA comparison reports higher electrical resistance during differentiation. Transport properties belong to the optimized multicomponent model; BUB1B is a separate cross-model cell-cycle gene-set observation. Neither is converted into drug-induced BUBR1 restoration, disease repair or MVA efficacy. Original dose, timing and density reporting conflicts remain explicit. Tretinoin is a verified molecular synonym only; retinol, generic RA and cis isomers do not inherit the result.

Two new portable files contain original interpretation and compact provenance. Raw papers, figures, supplements and reading logs remain external. Source and independent reviewers' full scientific reading is separate from root's selected primary-page/figure/caption/result checks. No RNA reanalysis, biological experiment or outer-reading credit is added. Only final-batch n418 changes; all original questions and other 46 objects remain. Twelve questions have direct qualified information: six study-level answers within limits, six with unresolved source components and 35 without direct answers. Original-queue totals remain 3,598 accepted / 5,759 combined outer / zero selected strings pending.

All FDA, SGC and NIH output tables remain unchanged. The report's recent-comparison paragraph is consolidated around measured endpoints. Track 1, its methods prefix, protected pitch/media/model/transfer/RNA and experimental plan remain unchanged. Private/draft/unsubmitted status and provisional nominees are unchanged.

## PDLSC conditioned-medium source comparison, 23 September 2026

The 207-assessment ledger, catalog, aliases, decision joins and every FDA/SGC/NIH inventory file remain byte-identical. The new evidence is retained as a biological-preparation comparison: controlled osteogenic/mineralization improvements under LPS co-exposure and a context-dependent BUB1B RNA summary. The supplied normalized signal is 2.642 times LPS alone, while conditioned medium alone is 0.706 times control. No constituent or shared culture reagent receives a new drug identity or inherited efficacy. An active molecule, standardized active dose, BubR1 protein/function and MVA benefit are unestablished.

Four portable files contain original qualification, a compact targeted expression record and a standard-library replay. The latter verifies exact native-workbook identity, four-condition columns, the target/filtered gene rows and six supplied ratios; it adds no missing P/FDR, sample map or raw-sequencing analysis. The original workbook and raw primary assets remain external and unchanged. A full-workbook import attempt was intentionally stopped without output; only the completed targeted native-cell reading/replay receives credit. Source, independent and root reading scopes remain distinct.

Only final-batch n64 is updated; the other 46 record objects and all 47 original questions/flags remain. Thirteen questions have direct qualified information: six study-level answers within limits, seven with unresolved components and 34 without direct answers. Original-queue totals remain 3,598 accepted / 5,759 combined outer / zero selected strings pending; earlier 2,161 annotations were not independently re-audited. No new outer-reading or compound-count credit is claimed.

Track 1, its methods prefix, protected pitch/media/model/transfer/RNA and experimental plan remain unchanged. Private/draft/unsubmitted status and provisional nominees are unchanged.

## CENP-F intermediates and KT-NB source update, 23 September 2026

Five material-specific assessments extend the ledger from 207 to 212, preserving every prior object. Four interventions share one CENP-F paper; the fifth is an exact-source Kalanchoe n-butanol fraction. E,E-farnesol partly restores CENP-F localization during zoledronic-acid co-exposure, while higher-dose geranylgeraniol reversal and the clodronate null are unshown reports. KT-NB causes cancer-cell growth suppression and mitotic defects, with affected WI-38 controls. Neither study supplies measured BUBR1 restoration or useful normal-cell/MVA rescue. Preparation, geometry, timing and figure-reporting limits remain explicit.

Four portable qualifications contain original interpretation and source metadata. Raw articles, figures, PDFs, reading logs and source workbooks remain external. Source, independent and root readers retain distinct actual scopes. The original 3,598 revalidation queue and 5,759 combined records are unchanged; no new outer-reading credit follows. Only final-batch n122/n340 change, leaving the other 45 objects and all 47 original questions/flags intact. Fifteen questions have direct qualified information, of which nine retain unresolved components; 32 lack a direct answer.

Decision joins and the three ingredient/probe inventory layers are replayed against the updated ledger; aliases and source catalogs remain unchanged. A name match confers no ingredient, formulation, availability or efficacy equivalence. All 109 PR1 IDs remain, with 100 unchanged objects and the same nine prior refinements. Protected Track 1, methods prefix, media, model/transfer/RNA and experimental-plan files remain unchanged. The repository remains private, both PRs drafts and unmerged, and Track 2 unsubmitted.

## CENP-E/Aurora published-source comparison, 23 September 2026

Four assessments for centrinone, MLN8054, TCS7010 and alisertib extend the ledger to 216. Existing GSK-923295, Aurora-balancing, Eg5-balancing and MG-132 entries gain separately attributed source context; the other 208 prior objects remain unchanged. The source supplies partial congression positives in non-transformed RPE-1 cells and a short-exposure alisertib null. Figure identity, pooled MLN8054 dose, genetic-assay schedule and cell-selection limitations remain explicit. No BUBR1 restoration, faithful daughter yield or MVA benefit is assigned.

Three portable files provide original qualification and compact numeric source summaries. Root and an independent reviewer separately read the targeted three workbooks and reproduced descriptive values, retaining source P strings and missing experiment/dish labels. No raw-image or inferential-statistical reanalysis is claimed. Original source articles, images, full workbooks and logs remain external. Source, independent and root actual reading scopes remain distinct.

Only final-batch n85 changes. The 15 exact-source qualifications remain (six study-level answers and nine partials), with one separately qualified published counterpart added. Thirty-one questions have neither kind of qualification; 32 still lack an exact-original-source answer. The preprint’s full text and exact version remain unverified. All original questions/flags, other 46 record objects and original-queue counts remain unchanged, with no new outer-reading or independent-study credit.

All 109 PR1 IDs remain: 97 objects unchanged and twelve qualified refinements. Decision joins and ingredient/probe inventories are replayed without changing catalog identities or granting preparation/availability/efficacy equivalence. Track 1, its methods prefix, current pitch/media, model/transfer/RNA and protected experimental plan remain unchanged. Both PRs remain private drafts and unmerged; Track 2 remains unsubmitted.

## Common-compound functional comparisons, 23 September 2026

Three distinct assessments bring the ledger to 219, preserving all 216 prior decision objects. Urolithin A contributes exploratory randomized human muscle-function positives with primary endpoints null in both trials. Putrescine dihydrochloride contributes mouse reproductive outcomes, while L-ascorbic acid reuses a previously qualified oocyte chromosome-protection source as a separate chemical identity. Neither preparation synonymy nor these other-model benefits establish BUBR1 restoration, MVA benefit or nominee superiority.

Seven portable files contain original source qualifications, structured evidence and one numerical human-outcome summary. Source-reviewer, independent and root reading scopes remain distinct. The putrescine 2013 source remains abstract-only and its 2015 native figures/supplements are unavailable; the ascorbate maturation denominator and urolithin adjusted clinical confidence intervals remain unresolved. Full papers, images and reading logs remain external.

All 47 original source-question records and the original reading queue remain byte-identical. No new original-queue reading, independent-study or global-completion credit follows. PR1 retention remains 109 IDs, with 97 unchanged objects and twelve previously qualified refinements. Catalog and ingredient joins are refreshed from the exact new aliases; their matches do not establish availability, preparation equivalence or efficacy.

Track 1, its methods prefix, current pitch/media, model/transfer/RNA and protected experimental plan remain unchanged. Both PRs remain private drafts and unmerged; Track 2 remains unsubmitted.

## Reproductive function and vitamin K2 molecular evidence, 23 September 2026

CoQ10, spermidine and mitoquinone bring the ledger to 222 assessments. Vitamin K2 is the only refined prior object; the other 218 are unchanged. Its original public figures and supporting caption supply Bub1b RNA and BUB1B-associated tissue fluorescence recovery accompanying behavioral protection. Exact preparation, dose, route, assay specificity and cellular/causal interpretation remain unresolved. The source main and Methods remain unread.

Nine portable files provide original qualifications, structured evidence and a compact spermidine source-data summary. Spermidine retains direct mouse chromosome-count and reproductive positives with a downward, non-FDR-significant Bub1b RNA contrast. Root inspected all populated cells of 16 targeted source sheets and the Bub1b row/header and reproduced descriptive means; no new P values, raw-sequencing analysis or full-main reading follows. CoQ10 retains mouse function and an oral human null with a nonpositive Bub1b RNA assay. Mitoquinone retains culture-development benefits and separate PRKN results without BUB1B or nuclear-chromosome evidence. Source, independent and root reading scopes are distinct.

All 47 historical source-question records and original reading queue remain byte-identical. All 109 PR1 IDs remain: 96 unchanged and thirteen qualified refinements, now including vitamin K2. No original-queue reading, independent-study, global-completion or MVA-efficacy credit is inferred. Name relations are regenerated without salt, class, product or availability equivalence.

Raw articles, original images, full workbooks and reading logs remain external. Track 1, its methods prefix, current media, model/transfer/RNA workflows and protected experimental plan are unchanged. Both PRs remain private drafts and unmerged; Track 2 remains unsubmitted.

## Existing-compound functional counterevidence, 23 September 2026

The ledger remains 222 assessments. Resveratrol, icariin and TG101209 are refined; all 219 other decoded decision objects are unchanged. Resveratrol preserves prior human power/torque positives while adding a shared Copenhagen cohort’s smaller aerobic gain and null knee-endurance augmentation. Icariin’s exact abstract supports motor/marker follow-up, while its BUB1b assay, tissue and exposure route remain unresolved. It receives priority assay clarification rather than an unsupported inherited low rank. TG101209 receives a distinct early post-challenge lung-protection role with two positive and four null respiratory contrasts, descriptive survival and no BUB1B measurement.

Six portable qualification files retain attributed reading and access scopes. Existing resveratrol qualification files now link the full-main counterevidence review; the original positive-source records remain intact. The thermal-assay review links the separate TG lung study without upgrading the original CETSA validation. A bounded independent primary/native check supported the resveratrol and TG interpretations; its icariin priority reassessment was not independent source validation. Raw mains, figures and reading logs remain external.

All 222 assessment IDs and all 109 PR1 IDs remain. PR1 retention is 94 decoded objects unchanged and fifteen qualified refinements. The original selection, reading counts and all 47 source questions are unchanged. FDA, SGC and NIH inventory builders were replayed against the refined ledger; all 14 output tables remain byte-identical. This adds no compound identity, original-queue reading, independent study, MVA efficacy or worldwide-completion claim.

Track 1, its methods prefix, current media, models, transfer/RNA workflows and experimental plan remain unchanged. The repo remains private, both PRs are draft and unmerged, and Track 2 remains unsubmitted.

## Pediatric neuromuscular compound comparisons, 23 September 2026

Three intervention assessments bring the ledger to 225: creatine monohydrate, source-labelled glutamine, and L-citrulline plus metformin hydrochloride. All 222 preceding decision objects and all 109 PR1 IDs remain unchanged; PR1 retention remains 94 equal objects and fifteen earlier refinements. Creatine and glutamine share one CINRG source, which confers no additional study or participant count. The combination does not modify or inherit monotherapy metformin/AICAR evidence.

Five portable files contain two source qualifications and a compact set of published outcome estimates. Creatine/glutamine remain abstract-and-registry-level, with unavailable scientific mains and native tables explicitly uncredited. An independent reader supported creatine’s substantive pediatric strength priority, preserving the six-month primary-strength null and functional nulls. The combination qualification covers the scientific main, tables, native figures and supplement with separate protocol/registry limits. Root checked the main/table evidence and selected native supplementary tables. Early nominal strength/distal-function positives, the final primary null, subgroup interaction, MRI slower worsening and 43 complete primary outcomes are retained. No new inferential statistics or participant-level reanalysis was performed.

All original queue and 47 source-question files remain unchanged. Inventory replay adds two NIH exact-name relations, preserving all source fields and the other twelve FDA/SGC/NIH output tables. Name relations do not prove stereochemistry, formulation, regulatory eligibility, availability or efficacy. Raw papers, native figures and reading logs remain external. No original-queue reading, independent-study, MVA-efficacy or global-completion claim follows.

Track 1, its methods prefix, current media, model/transfer/RNA workflows and experimental plan are unchanged. The repository remains private; both PRs are draft and unmerged; Track 2 remains unsubmitted.

## Common muscle-compound comparisons, 23 September 2026

Two assessments bring the ledger to 227: L-leucine and calcium-HMB. All 225 previous objects and all 109 PR1 IDs remain, with PR1 retention unchanged at 94 equal objects and fifteen earlier refinements. Five portable files contain two source qualifications and compact published HMB outcome estimates; no raw papers, native figures or participant data are added.

Leucine qualification retains nominal walking/expiratory-pressure preservation from a 42/50 complete-case trial alongside LACE’s imprecise overall primary null. Positive subgroup estimates retain the published nonsignificant interaction values; no reanalysis or equivalence inference is claimed. LACE’s source-labelled leucine is not assigned unverified stereochemical identity. Calcium-HMB qualification retains a secondary gait-speed benefit, the protocol-specified calcium-unmatched control, primary mass interaction limits, and one cohort across two analyses. Its result mains remain inaccessible. Independent and root reading scopes are explicitly separate from the source reviewers’ scopes.

Inventory replay adds one FDA L-leucine name link and four NIH name relations in three groups. All source fields and the other eleven FDA/SGC/NIH output tables remain unchanged. Two punctuation aliases do not create distinct compounds. All original queue and 47 source-question files remain unchanged. No new independent-study, original-queue, MVA-efficacy or worldwide-completion credit follows.

Track 1, its methods prefix, current media, model/transfer/RNA workflows and experimental plan are unchanged. The repository remains private, both PRs draft and unmerged, and Track 2 unsubmitted.
