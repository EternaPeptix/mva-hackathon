# jeweled — current research reports

**Revised 21 September 2026.** Track 1 model 1 has a confirmed automated score; its final report-persistence check is pending as detailed below. Track 2 remains a review draft. No clinical sign-out has occurred. The repository remains private; the new Skylar pitch awaits human review before public upload.

**[Open the Track 2 review package](REVIEW.md)** for the current report, evidence,
new video, editable slides, captions and clean-state reviewer assessment.

## Track 1 submission record

On **21 September 2026 at 17:41 UTC**, account/team **jeweled** sent one Track 1
submission: **model1_jeweled_evidence-reconciled-v2.csv**. The official public
[leaderboard](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026)
confirmed **100.0 rank points and F-max 1.000** at 17:42:54 UTC. This is the
prediction score; qualitative methods judging is separate.

The [exact report upload](jeweled_track1_model1_upload.md) combines the reviewed
Track 1 report and all 12 Track 1 methods answers, preserves the AI disclosure,
adds the required acknowledgement, and replaces pre-submission status wording.
The [CSV](jeweled_evidence-reconciled-v2.csv) remains byte-identical, with one
BUB1B pair and subjective EPCR 0.8. The submitted repository URL pins
[revision cda0b0c12249b446adeff61cd1838f148ccb1768](https://github.com/EternaPeptix/mva-hackathon/tree/cda0b0c12249b446adeff61cd1838f148ccb1768).
The separate report/methods drafts below retain their pre-submission source state.

**Verification boundary:** both temporary portal uploads were read back and
matched the local files byte for byte (CSV 480 bytes; report 26,390 bytes).
The leaderboard verifies stored metadata and score. The original UI became
unreadable while scoring, so the final “Submission received” response and durable
report save have not yet been independently confirmed. The service stores metadata
before its separate report upload. **Do not resubmit merely to resolve this check.**
Only one Submit action was sent; the pre-submit counter was 0/6.

## Current interpretation

The current [Track 1 report](jeweled_track1_report_submission.md), [Track 2 report](jeweled_track2_report_submission.md) and [methods-form answers](jeweled_methods_form_answers.md) govern interpretation. Track 1 adopts **subjective exact-pair EPCR 0.8** for BUB1B Leu737Ter plus Asn1002Lys. This replaces earlier judgments without new segregation or functional data. Trans is favored under the recessive model but remains unobserved; the missense remains a VUS. The [causal interpretation appendix](appendices/causal-interpretation.md) explains the evidence policy and remaining alternatives.

Track 2 selects **everolimus as the lead drug hypothesis and niacin (nicotinic acid) as the primary comparator**, tested separately rather than as a combination. The two routes are mTORC1 normalization and restoration of usable BUBR1. Everolimus has the clearer drug-to-target implementation and indication-specific pediatric pharmacology; niacin supplies an orthogonal proximal hypothesis through human muscle NAD engagement and separate SIRT2/BUBR1 experiments. Arimoclomol is a restoration reserve; dasatinib is a senolytic reserve requiring harmful burden and selective susceptibility. Prioritizing proximal causal evidence would put niacin first, so the pair is more defensible than a universal drug order. No relevant drug rescue is demonstrated; genomic EPCR is not a treatment-response probability.

The [one-row prediction CSV](jeweled_evidence-reconciled-v2.csv) records the current decision; its [provenance](jeweled_evidence-reconciled-v2.csv.provenance.json) covers format and reference checks, not a comparison with private truth.

The completed [public mouse RNA analysis](../workflow/track2-rna/README.md) did not corroborate an upregulated mTORC1 Hallmark program. That result does not rebut the published substrate-phosphorylation measurements. Underreported rapamycin nonbenefit and a different-cell-state everolimus-associated BUBR1 decrease remain concerns to test. The common comparison requires absolute functional benefit, preserved renewal and no worsened fidelity; only a claimed restoration mechanism requires corresponding BUBR1/segregation rescue.

## Pitch and delivery

The [current pitch](jeweled_track2_pitch_reconciled.md) is a new eight-slide, three-minute Skylar recording for the everolimus/niacin comparison. It replaces the earlier Lauren arimoclomol/NAC version as the proposed submission video. The old private release is archival. The new [movie](media/jeweled_track2_skylar.mp4), [editable deck](media/jeweled_track2_slides.pptx) and [captions](media/jeweled_track2_skylar.srt) are included here for human review. No new public video upload or Track 2 submission is claimed.

The [route decisions](appendices/track2-route-decisions.md) compare the alternatives and specify the result that would change each decision. A fresh reviewer with no conversation context checked the report and pitch against the official rubric. The revised main report explains the proposal in plain language; detailed experimental and source qualifications remain in appendices.

A subsequent [clean-state rescore](reviews/track2-clean-state-20260921.md) assessed
the corrected artifacts without seeing earlier reviews or scores. It is an
unofficial assessment, separate from the initial reviewer's correction check.
The current text additionally incorporates the human-review correction on
Pirinen's treated-control strength gains and the slide's patient-age description;
the archived score predates those changes.

Claude Code and Codex used personal subscriptions with training disabled, as confirmed by the user. Cartesia Sonic 3.6/Skylar uses the Free plan, verified on 21 September; its training/retention settings remain unverified. The cover is an AI-generated conceptual illustration.

## Standalone report files

The reports above use relative links for browsing this repository. For a separate
Markdown upload, export the two reports, methods answers and unchanged prediction
CSV from a committed revision:

```sh
python3 workflow/scripts/export_submission_reports.py ../mva-submission-ready
```

Run from the repository root and choose a new directory outside the checkout.
The exporter reads committed `HEAD` by default (`--revision` selects another
commit), changes relative file links to GitHub links pinned to that revision,
and adds a source header. Report prose and external citations are preserved;
the CSV is copied byte for byte. Uncommitted edits are not exported. The command
does not upload files, change repository visibility or grant access to private
links. Human review of the new pitch and a public video URL are still required before submission.

## Methods and retained appendices

The [public benchmark method](../workflow/docs/PUBLIC-BENCHMARK-METHOD.md) states the HG002 chromosome-20 scope, measured coverage, confident-region evaluation and separate SNV/indel results. A fresh full Mac Studio replay reproduced the primary counts and audit in 134.3 minutes; [reproduction documentation](../REPRODUCIBILITY.md) distinguishes that executed source from later launcher safeguards and limited Linux checks. Historical case runs still require authorized data and their recorded prerequisites.

The [Track 1 case-evidence map](appendices/track1-case-evidence.md) presents actual
prioritizer ranks, predictor sensitivities, candidate fragment support, locus and
ten-gene counts, depth and phase. Its [selected numerical results](../results/track1-selected-case-evidence.json)
separate supplied-VCF rankings from later four-lane native-call review and link
the claims to authorized reconstruction steps. Restricted full case outputs are
not included.

The new [Track 2 evidence pipeline](appendices/track2-evidence-pipeline.md) governs the lead selection, alternatives and proposed experiments. The imported [genome-informed analysis](../workflow/track2-model/README.md) reconstructs nominated coding consequences and illustrates population accounting under explicit assumptions; it does not calibrate patient risk or rank drugs.

The following twelve appendices preserve dated method detail. Their original results or proposed experiments remain useful, but historical status and interpretation statements do not override the current reports and governing evidence pipeline.

| Topic | Appendix |
|---|---|
| Secondary-variant inventory | [Secondary findings](appendices/secondary-findings.md) |
| Secondary-panel coding coverage | [Depth review](appendices/secondary-depth.md) |
| Locations and transcript membership of coverage gaps | [Saved-depth gap map](appendices/secondary-depth-gaps.md) |
| Nuclear-mixture model and limitations | [Nuclear-mixture review](appendices/nuclear-mixture.md) |
| Exact-allele public frequency fields | [Population-frequency review](appendices/population-frequency.md) |
| Public SV benchmark misses | [SV residual diagnostics](appendices/public-sv-residuals.md) |
| BUBR1-to-drug evidence | [Mechanistic bridge](appendices/mechanistic-bridge.md) |
| NAC analytes and exposure | [NAC exposure](appendices/nac-exposure.md) |
| Functional endpoints and decisions | [Assay decision plan](appendices/assay-decision-plan.md) |
| Prospective experimental design | [Pilot protocol](appendices/pilot-protocol.md) |
| First-division endpoint accounting | [Founder-ledger specification](appendices/founder-ledger.md) |
| Tumor protection and treatment interactions | [Oncology counter-screen](appendices/oncology-counterscreen.md) |

See [provenance and evidence scope](PROVENANCE.md) for what is included, omitted and still required before publication.
