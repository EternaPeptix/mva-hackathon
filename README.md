# MVA research workflows

This is the private reproducibility repository for team **jeweled**'s MVA hackathon
research. It contains the analysis code, configurations, recorded Mac Studio
environments, and selected public/synthetic validation evidence. The clinical
analysis and Track 1/Track 2 submissions remain unfinished.

The working hypothesis and eventual drug-repurposing proposal require careful
variant, phase, functional and exposure review. No diagnosis or treatment benefit
is established by this repository. The public benchmark is a workflow check, not
clinical sensitivity for the child.

## Start here

- [Reproduction and data requirements](REPRODUCIBILITY.md)
- [Apple Silicon environment capture](workflow/environments/20260919/runtime.json)
- [Frozen small-variant comparison](workflow/configs/public-small-variant-comparison.json)
- [Original public benchmark result](results/germline-public-comparison-20260919.json)
- [Exploratory multiallelic comparison](results/germline-retained-multiallelic-comparison-20260919.json)
- [Exact truth-identity and scope reconciliation](results/germline-representation-reconciliation-20260919.json)
- [Native CRAM compatibility control](results/gatk-cram-compatibility-native-control-20260920.json)

The workflow covers reference preparation, alignment, native GATK/bcftools calling,
phenotype-assisted and phenotype-independent prioritisation, read evidence, phasing,
copy number, homozygosity, structural discovery, public controls, and report/pitch
production. Some components are completed controls or prepared follow-ons; presence
of a script does not mean it has produced a patient result.

The tested platform is **Apple Silicon macOS on a Mac Studio**. Conda specifications
are for `osx-arm64`; Linux, x86 and GPU-only paths are not substitutes for this setup.
The older Docker/three-caller scripts remain historical and are not the active method.
Use the native entry points described in the reproduction document.

## Measured public checks

On public HG002 chromosome 20, over the identical 82,589 GIAB v4.2.1 truth records:

| Caller | Original split-call recall / precision | Exploratory retained-genotype recall / precision |
|---|---:|---:|
| GATK | 99.1040% / 98.5365% | 99.3680% / 99.0537% |
| bcftools | 98.2867% / 97.8593% | 98.5289% / 98.3415% |

The retained-genotype arm was motivated by observed discordances. It uses the same
comparison rules and does not overwrite the original result. Splitting complex
genotypes can introduce conflicting reference assertions or obscure a wrong second
allele. Both gains and losses are counted; changed call-record denominators are not
one-to-one counts of corrected biological errors. Neither result establishes phase
or causal probability for the child.

Twelve original native comparator controls, six audit controls, six multiallelic
mechanism controls, five metadata checks and three interval-geometry checks passed.
Separately, seven reconstructed environments passed 22 selected existing native
controls. These are scoped checks, not a claim that every workflow has been reproduced
end to end from a fresh public download.

The pinned GATK/HTSJDK reader rejects CRAM 3.1. A separate recovery path creates
CRAM 3.0 and directly compares every decoded record, including all FLAG bits and
auxiliary data, before calling. Its native 195-record synthetic control recovers
the known heterozygous SNV and deletion with both callers, and rejects deliberate
quality, tag, duplicate-flag and record-order changes. The direct check covers a
flag-preservation gap reproduced in the pinned samtools checksum implementation.

## Data and release status

Raw reads, alignments, genotype-scale tables, patient-derived result files, tool
archives and downloaded references are absent from this source snapshot. Configurations
include permitted HPO terms and a bounded list of named candidate sites needed to
reproduce the targeted analyses. Included empirical result files use public HG002
or fictional controls; runtime records describe software reconstruction.

Access to the challenge dataset must be obtained separately through its authorised
process. The organisers' [data-handling clarification](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/2)
distinguishes retained findings/code from genomic data and its reformats. The final
reports, reviewed candidate CSV and pitch URL will be added after their remaining
analysis and review. This repository stays private during the current research phase;
the required later public release remains a separate step.

`SOURCE-SNAPSHOT.json` records source and release hashes. It also records a
comment-only substitution of an unnecessary acquisition-specific read header/barcode
example in the exported alignment script. No executable lines change. The internal
Git history and private analysis outputs were not copied. Historical source hashes
in execution receipts refer to the original analysis files, not an assertion that
every later documentation or comment revision has identical bytes.

**AI disclosure:** Anthropic Claude Code and OpenAI Codex, personal subscriptions,
training on the user's content disabled, user-confirmed 19 September 2026. Exact
historical model/plan names remain unverified. This records the user's configuration;
it is not an independent audit of all historical provider retention or feedback settings.

Code is covered by the [MIT license](LICENSE); written reports and documentation use
CC BY 4.0 under its scope note. Third-party software and reference data retain their
own licenses and access requirements.
