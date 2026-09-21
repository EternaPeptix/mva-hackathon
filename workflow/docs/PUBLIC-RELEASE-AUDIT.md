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
