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
- [Python builds from published sources](results/runtime-published-source-build-20260920.json)
- [Fresh runtime validation](results/runtime-published-source-validation-20260920.json)
- [Public structural comparison reproduced with the new runtime](results/runtime-published-source-public-sv-20260920.json)
- [Pinned public CNV context resources](results/cnv-context-public-sources-20260920.json)
- [CNV context representation controls](results/cnv-context-method-controls-20260920.json)
- [Structural/CNV selection controls](results/structural-cnv-context-controls-20260920.json)
- [Literal tandem-junction and transcript-geometry controls](results/tandem-junction-method-controls-20260920.json)
- [Nearby-pair phase, mapping and mate-handling controls](results/paralog-pair-method-controls-20260920.json)
- [Prospective assay-capacity calculation](results/assay-capacity-planning-20260920.json) — planning assumptions, not an experiment or power estimate
- [SNV cluster and fragment controls](results/snv-cluster-method-controls-20260920.json)
- [Public sequence comparison controls](results/public-clone-method-controls-20260920.json)
- [Target-base orientation controls](results/cluster-target-projection-controls-20260920.json)
- [Four-lane indel sequence/genotype controls](results/four-lane-indel-controls-20260920.json)
- [Indel window and orientation diagnostics](results/four-lane-indel-context-controls-20260920.json)
- [Fixed-region depth and missingness controls](results/chromosome-depth-method-controls-20260920.json)
- [Regional-ploidy interval controls](results/ploidy-region-method-controls-20260920-v2.json)
- [Diploid/haploid calling on fictional sequence](results/ploidy-calling-native-controls-20260920.json)
- [Ploidy-boundary genotype and ownership controls](results/ploidy-boundary-method-controls-20260920.json)
- [Native boundary reproduction and fixed-ploidy recovery](results/ploidy-boundary-native-controls-20260920.json)
- [Pinned GATK regional-ploidy semantics](results/ploidy-boundary-public-source-review-20260920.json)
- [Versioned mitochondrial references and indexes](results/mitochondrial-public-resource-verification-20260920.json)
- [Mitochondrial geometry and read-preservation controls](results/mitochondrial-method-controls-20260920.json)
- [Circular-read native calling control](results/mitochondrial-native-control-20260920.json)
- [Original multi-allele context during mitochondrial comparison](results/mitochondrial-review-method-controls-20260920.json)
- [Mitochondrial read-context and coordinate controls](results/mitochondrial-context-method-controls-20260920.json)
- [Pinned mitochondrial mixture-screen runtimes and sources](results/mutserve-runtime-source-review-20260920.json)
- [Native caller/importer contracts on fictional reads](results/mutserve-contract-controls-20260920.json)
- [Mixture-screen selection and allele-rank controls](results/mitochondrial-mixture-method-controls-20260920.json)

The workflow covers reference preparation, alignment, native GATK/bcftools calling,
phenotype-assisted and phenotype-independent prioritisation, read evidence, phasing,
copy number, homozygosity, structural discovery, public controls, and report/pitch
production. Some components are completed controls or prepared follow-ons; presence
of a script does not mean it has produced a patient result.

The revised local pitch tools include segment-by-segment decoding and caption/audio
integrity checks. The [media reproduction recipe](workflow/docs/pitch-v2-reproduction.md)
documents the tested Mac Studio environment, storyboard contract and limits. The
case-specific storyboard, clinical reports and media remain outside this source snapshot.

The [changed-call read-review recipe](workflow/docs/changed-call-read-reproduction.md)
adds original ploidy/allele-preserving context, competing-alternate fragment counts
and exact-window sensitivities. Three fictional native-pysam controls qualify this
path; selected case inputs and results remain separately controlled.

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

The source-build route now obtains 40 pinned publisher wheels and builds the other
five packages from pinned published source archives without reading the old cached
wheels. Seven further fresh environments pass the same 22 controls; six additional
native-extension cases match the old runtime. Repeating both public structural
comparisons reproduces every compared allele, genotype and decision, including
69/802 matched DELLY truth alleles and 171/802 SvABA truth alleles. Those repetitions
reuse the original caller outputs. They do not close the structural sensitivity gap
or establish clinical accuracy. See the installation recipe and limits in
[REPRODUCIBILITY.md](REPRODUCIBILITY.md).

The pinned GATK/HTSJDK reader rejects CRAM 3.1. A separate recovery path creates
CRAM 3.0 and directly compares every decoded record, including all FLAG bits and
auxiliary data, before calling. Its native 195-record synthetic control recovers
the known heterozygous SNV and deletion with both callers, and rejects deliberate
quality, tag, duplicate-flag and record-order changes. The direct check covers a
flag-preservation gap reproduced in the pinned samtools checksum implementation.

Copy-number context code now compares nominal intervals with the complete gnomAD
v4.1 SV sites release and dated ClinGen dosage tables, and retains every common-panel
marker including missing/excluded evidence. Nine focused representation controls
pass on the Studio. Population allele counts and multiallelic sample copy-state
distributions remain distinct; overlap alone does not classify a case allele.
The corresponding patient-derived result tables and figures are excluded.

A fixed-region depth workflow preserves explicit zeros and separate reference-span
and concrete-base denominators. Four native controls cover selection, missing
positions, overlapping mates and filters. The regional-ploidy runner has three
interval controls; both GATK and bcftools recover a diploid SNV, haploid SNV and
haploid insertion on fictional sequence. The initial index-collision source is
preserved separately. These controls do not validate a karyotype, mosaicism
sensitivity, mitochondrial interpretation or clinical negative result. Patient
depth measurements and case output files are outside this export.

A boundary-specific native probe reproduces GATK 4.6.1.0's documented use of the
maximum ploidy overlapping an active region. It plants 12 SNVs across three
transitions, including exact interval endpoints. Six regional genotype-shape
exceptions recur; independent constant-ploidy calls and the second caller recover
all 24 expected caller/allele genotypes. Four focused tests check original
multiallelic genotype/likelihood preservation, invalid and partial evidence, and
explicit handling of boundary-spanning records. The fixed-ploidy runner preserves
the original output and reuses a checksum-identified completed bcftools VCF.

The mitochondrial adaptation verifies the literal rCRS sequence and 8,000-base
rotation, preserves primary read sequence/quality, and keeps original multi-allele
context through calling, liftover and comparison. Nine focused checks and four
planted alleles in 30,000 synthetic circular read pairs pass. Initial filters assume
zero contamination; the method has not validated that assumption, nuclear-copy
discrimination or very low heteroplasmy. Case calls, target-level read counts and
interpretations remain outside the source export.

A further context helper preserves CIGAR-aware mismatch and read-cycle evidence,
BAQ missingness, reciprocal mates and reported full-reference remapping alternatives.
Seven focused controls pass. Its symmetric read-context filters are exploratory;
they do not establish a heteroplasmy detection limit or exhaustive NUMT exclusion.
Detailed read identities, sequences and case-derived summaries are not exported.

## Data and release status

Raw reads, alignments, genotype-scale tables, patient-derived result files, tool
archives and downloaded references are absent from this source snapshot. Configurations
include permitted HPO terms and a bounded list of named candidate sites needed to
reproduce the targeted analyses. Included empirical result files use public HG002
or fictional controls; runtime records describe software reconstruction. The CNV
source receipt contains public download identities and metadata, not case matches.
The structural/CNV control receipt covers fictional genotype/filter/query-boundary
selection cases; the corresponding patient call comparisons remain outside this export.

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

## Autosomal discovery handoff

The [completed-shard handoff audit](workflow/docs/discovery-handoff-reproduction.md)
checks native allele/genotype preservation, the existing screen's quality-field
compatibility, and complete-record normalization replay before downstream ranking.
Five [native fictional controls](workflow/source-snapshots/discovery-handoff-controls-20260920/controls.json)
pass, including wrong-index and changed-input refusal. Pending chromosomes remain
explicit. This adds a method check; it does not establish a complete genomic analysis,
new causal candidate, quality policy or clinical sensitivity. Case audit outputs are
not included in this source snapshot.

## Native coding inventory and descriptive models

The [native coding review](workflow/docs/NATIVE-CODING-REVIEW.md) links each
protein-altering allele to its complete original caller record, including
multiallelic genotype and likelihoods. Missing or low quality remains visible;
no PASS label or GQ is imputed. Caller-separated recessive combinations retain
frequency uncertainty and source-scoped phase.

Eighteen [native fictional checks](workflow/source-snapshots/native-coding-controls-20260920/controls.json)
pass, including exact sequence-equivalent normalization, complete annotation
payload preservation, original multiallelic recovery, input-tamper refusal and
explicit partial scope. These checks do not establish biological sensitivity or
a complete candidate ranking. Case inventories remain outside this source export.

## Joint sequence hypotheses and fragment evidence

The [joint-haplotype review](workflow/docs/JOINT-HAPLOTYPE-REVIEW.md) enumerates
compatible sequence hypotheses within each caller, including retained uncalled
ALTs. Identical sequences preserve all source lineages while sharing one
observed fragment count. Original read groups, unmatched sequences, mate
conflicts and raw-quality sensitivity remain explicit; constructed sequences
do not become genotype or parental-phase calls.

Twenty [fictional helper and native integration controls](workflow/source-snapshots/joint-haplotype-controls-20260920/controls.json)
passed on Mac Studio 512S2. They exercise literal sequence equivalence, native
CRAM extraction, original caller context and constructed-sequence annotation.
They do not establish clinical sensitivity, causality or treatment benefit.

The private source configuration includes selected targets, sample labels and
input-identity hashes needed to reproduce the exact method. Individual read
records, read names, case evidence/result payloads, clinical reports and media
are excluded from this source export.
