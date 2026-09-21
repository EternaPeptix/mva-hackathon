# jeweled — current research reports

**Revised 21 September 2026.** These are review drafts, not a competition submission or clinical sign-out. The repository and approved GitHub media copy remain private; public visibility and Vimeo delivery have not been completed.

## Current interpretation

The current [Track 1 report](jeweled_track1_report_submission.md), [Track 2 report](jeweled_track2_report_submission.md) and [methods-form answers](jeweled_methods_form_answers.md) govern interpretation. Track 1 adopts **subjective exact-pair EPCR 0.8** for BUB1B Leu737Ter plus Asn1002Lys. This replaces earlier judgments without new segregation or functional data. Trans is favored under the recessive model but remains unobserved; the missense remains a VUS. The [causal interpretation appendix](appendices/causal-interpretation.md) explains the evidence policy and remaining alternatives.

Track 2 selects **everolimus as the lead drug hypothesis and niacin (nicotinic acid) as the primary comparator**, tested separately rather than as a combination. The two routes are mTORC1 normalization and restoration of usable BUBR1. Everolimus has the clearer drug-to-target implementation and indication-specific pediatric pharmacology; niacin supplies an orthogonal proximal hypothesis through human muscle NAD engagement and separate SIRT2/BUBR1 experiments. Arimoclomol is a restoration reserve; dasatinib is a senolytic reserve requiring harmful burden and selective susceptibility. Prioritizing proximal causal evidence would put niacin first, so the pair is more defensible than a universal drug order. No relevant drug rescue is demonstrated; genomic EPCR is not a treatment-response probability.

The [one-row prediction CSV](jeweled_evidence-reconciled-v2.csv) records the current decision; its [provenance](jeweled_evidence-reconciled-v2.csv.provenance.json) covers format and reference checks, not a comparison with private truth.

The completed [public mouse RNA analysis](../workflow/track2-rna/README.md) did not corroborate an upregulated mTORC1 Hallmark program. That result does not rebut the published substrate-phosphorylation measurements. Underreported rapamycin nonbenefit and a different-cell-state everolimus-associated BUBR1 decrease remain concerns to test. The common comparison requires absolute functional benefit, preserved renewal and no worsened fidelity; only a claimed restoration mechanism requires corresponding BUBR1/segregation rescue.

## Pitch and delivery

The [pitch document](jeweled_track2_pitch_reconciled.md) contains a proposed text-only revision explaining the everolimus/niacin comparison and decisive experiments, followed by the unchanged approved narration. The approved 180-second Lauren / Cartesia Sonic 3.6 movie and its 24 caption segments use the previous 293-word arimoclomol/NAC script and remain [private GitHub release attachments](https://github.com/EternaPeptix/mva-hackathon/releases/tag/track2-pitch-lauren-20260920); access requires repository permission. That recording does not present the current everolimus/niacin comparison. The proposed revision is not recorded or approved. Media work and Vimeo are on hold at the user's request; competition submission remains pending.

Claude Code and Codex used personal subscriptions with training on the user's content disabled, as confirmed by the user. Cartesia plan/tier and training/retention settings remain unverified and are disclosed separately.

## Methods and retained appendices

The [public benchmark method](../workflow/docs/PUBLIC-BENCHMARK-METHOD.md) states the HG002 chromosome-20 scope, measured coverage, confident-region evaluation and separate SNV/indel results. A fresh full Mac Studio replay reproduced the primary counts and audit in 134.3 minutes; [reproduction documentation](../REPRODUCIBILITY.md) distinguishes that executed source from later launcher safeguards and limited Linux checks. Historical case runs still require authorized data and their recorded prerequisites.

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
