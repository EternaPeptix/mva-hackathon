# Report provenance and scope

**Review revision: 21 September 2026.** This directory brings the revised research reports, methods answers, prediction CSV, proposed pitch revision and unchanged approved narration into one browsable checkout. It is not a new competition submission, a public release or a new scientific run.

## Evidence scope

The current [Track 1 report](jeweled_track1_report_submission.md), [Track 2 report](jeweled_track2_report_submission.md) and [methods form](jeweled_methods_form_answers.md), together with the [causal interpretation appendix](appendices/causal-interpretation.md), state the adopted interpretation. The subjective exact-pair EPCR is **0.8**. Its revision from 0.35 reflects reappraisal of the supplied clinical context and completed evidence, without new segregation or functional data. Unknown phase is retained; N1002K remains a VUS. No treatment response is claimed.

The report and script texts are synchronized with the current internal research documents; local navigation differs, including an explicit label for the omitted detailed Track 2 revision. The methods form retains all 22 original question headings. The pitch document preserves the approved six-section narration, 24 caption segments and 293 words, and adds a clearly separated text-only proposal after review identified missing explanation of the abundance test and drug ranking. This proposal is not recorded or approved and does not add or alter the media.

The prediction CSV and its provenance are byte-identical to the revised internal export. The CSV SHA-256 is `40232aa87da1b498d44e73d82caf485aa390c7d02547bf4f817c557495cf1355`. Its one row is a manually reviewed research prediction, with phase unknown; the recorded builder checks establish format and reference compatibility, not clinical validation or private-truth accuracy.

The twelve dated appendices were selected from the earlier curated review package. Their scientific bodies are preserved, with a notice explaining their historical status. Links labelled “underlying file not included” lead here deliberately: the named receipt, ledger or output is **not** included or independently reproducible from this directory. These links document a scope limit, not a substitute download or an assertion that the omitted file exists in this checkout.

Excluded material includes raw case reads, read names, alignments, pileups, raw/genome-scale case VCFs, full ledgers, large interval arrays, execution logs, machine-specific working paths, credentials and internal conversation/Git history. The included reports still disclose selected candidate-level findings and must be handled consistently with the organizers' data-use terms. Curating a report does not itself authorize public disclosure.

## Archival material

The retained appendices describe their original dates and analysis scopes. Statements such as “unchanged EPCR,” “no submission change” or “frozen package” refer to those original analyses; the current reports and [current index](README.md) take precedence. Older Track 2 revision 11, earlier probability judgments, detailed working reports, approval receipts and prior package manifests are not copied here. Historical clinical or mechanistic assertions are not silently reinstated by an appendix link.

[REPRODUCIBILITY.md](../REPRODUCIBILITY.md) and the workflow sources are an archive of methods and their declared qualifications. Some exact historical launch plans depend on native paths, separately provisioned public resources or authorized case inputs. Those plans are not fresh-clone execution recipes. The [public HG002 benchmark method](../workflow/docs/PUBLIC-BENCHMARK-METHOD.md) separates measured public RTG vcfeval results from representation sensitivities and does not claim an additional hap.py comparison. Source/control evidence and the selected scientific summaries have different scopes.

## Media and AI disclosure

The [private Lauren release](https://github.com/EternaPeptix/mva-hackathon/releases/tag/track2-pitch-lauren-20260920) holds the approved movie and captions; no binaries are copied into this report directory. Its narration does not quote EPCR and presents drug effects as hypotheses, but underexplains the current report's lead mechanism and ranking. The proposed text revision addresses those presentation gaps. The old movie does not contain that text; no new narration or rendering has occurred. Vimeo is on hold at the user's request, and competition submission remains pending. Private GitHub attachment delivery is not public publication.

Claude Code and Codex were used on personal subscriptions with training on the user's content disabled, according to the user's confirmation. Exact historical model/plan names are unverified. Candidate-level evidence and report text were reviewed with hosted assistants. Cartesia Sonic 3.6 supplied synthetic stock-Lauren narration; its plan/tier and training/retention settings are separately unverified. The Claude/Codex confirmation does not establish Cartesia settings.

## Acknowledgement

The [organizer rules](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/rules.py) prescribe this acknowledgement:

> This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.

## Dataset provenance

The case context and data provenance are described by the [organizer overview](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/about.py), the [challenge dataset resource](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data) and the official Synapse project [“Rare Disease, Real Kid. The 2026 MVA Hackathon” (syn76251147)](https://www.synapse.org/Synapse:syn76251147/wiki/642892). The rules require the dataset reference supplied on that Synapse page at publication time. The current public page supplied no DOI or formatted dataset citation; recheck it before publication rather than inventing one. These provenance links do not imply that the case data are public.
