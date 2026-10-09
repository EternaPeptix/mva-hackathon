# Reproduction scope and order

Start with the public HG002 benchmark below to exercise download, calling and comparison without any private case inputs. Historical case-workflow restoration follows separately.

## Public benchmark quick start

From a fresh clone, with Bash and system Python 3.9+ installed:

```bash
bash workflow/public-benchmark/run.sh ../hg002-public-replay
```

The new [public replay entrypoint](workflow/public-benchmark/README.md) downloads
checksum-pinned public HG002 chr20 reads, GRCh38 reference and GIAB v4.2.1 truth;
installs an isolated, explicitly locked runtime; runs the existing twelve RTG
controls; prepares the reference; runs the unchanged GATK/bcftools caller over
the full public chromosome; then runs the genotype-aware comparison and
accounting audit. It preserves full chromosome context and
official confident-region evaluation. This is a prealigned caller-component
benchmark, not FASTQ alignment or clinical validation. It has no private-input
requirement and does not copy historical performance numbers into a new result.

Choose a new output tree outside the clone. Public data total 2.15 GB, with about
another 1 GB for runtime archives; plan 15–25 GiB working storage and at least
40 GiB free. The full Mac Studio replay took 134.3 minutes, including downloads,
setup, calling, comparison and audit; another host's runtime is unmeasured. The default two-worker run has
a six-hour, 32-GiB sampled group-memory and 30-GiB output guard. Only its own
process group is controlled. Linux x86-64/aarch64 and macOS arm64 have separate
package locks; package resolution alone does not establish execution on Linux.
The [replay guide](workflow/public-benchmark/README.md) gives exact scope,
prerequisites, limits, output locations and a clearly labelled controls-only mode.

On 21 September 2026, the full macOS arm64 replay completed with exit 0 after downloads, setup, twelve existing RTG controls, both callers, comparisons and accounting audit. It exactly reproduced the original primary split-call counts and audit partitions over 82,589 scored truth records. Total wall time was 8,057.14 seconds (134.3 minutes), with peak sampled process-group memory 5.31 GiB and output files 12.65 GiB. These are observations from one run, not resource guarantees.

Linux x86-64 under QEMU completed the pinned environment, all twelve controls, and actual GATK/bcftools calling plus both RTG comparisons on a predefined public 100-kb interval. Full Linux chr20 replay remains unverified; Linux aarch64 has a resolved lock only.

The final launcher adds shell-safe Java-path quoting and inherited-JVM-option sanitation after the Mac run started. Those safeguards received separate syntax/path checks; the complete run proves the recorded launched source, not full execution of the final hardened launcher. Scientific helper bodies and settings are unchanged. The [platform record](results/public-benchmark-replay-20260921.json) retains the executed/final source identities and terminal evidence.

## Historical authorized-case workflow (Apple Silicon)

This snapshot supplies the workflow source and configuration files. Reproduction
requires the separately acquired inputs and software below. Exact executed commands,
input identities and output counts are retained in the selected JSON receipts.
Paths in historical receipts point to the original Studio and must be mapped to a
new project root; they are provenance, not portable shell commands to run blindly.

## 1. Prepare a separate Apple Silicon environment

Use a new project root with the source checkout in `ROOT/repo`, references in
`ROOT/ref`, tool installations in `ROOT/tools`, and outputs outside the checkout.
The tested environment was macOS 26.4.1 arm64, Python 3.12.14, samtools/bcftools 1.24,
GATK 4.6.1.0 with Java 21.0.10, and BWA 0.7.19.

The [environment inventory](workflow/environments/20260919/runtime.json) contains four
Conda archive specifications and three Python environments. Install Micromamba 2.9.0
from its official distribution, then use the exact `*-osx-arm64.explicit.txt` files:

```sh
MVA_ROOT=/path/to/new/mva-project
MVA_REPO="$MVA_ROOT/repo"
MVA_MM=/path/to/micromamba
for MVA_ENV in bio delly-2.6.0 mafft-7.526 svaba-1.2.0; do
  "$MVA_MM" create --no-rc --no-env --yes --always-copy --no-pyc \
    --safety-checks enabled --extra-safety-checks \
    --root-prefix "$MVA_ROOT/tools/mamba" \
    --prefix "$MVA_ROOT/tools/$MVA_ENV" \
    --file "$MVA_REPO/workflow/environments/20260919/$MVA_ENV-osx-arm64.explicit.txt"
done
```

The recorded native reconstruction used verified cached archives in offline mode;
the example permits fetching the same public archives. Install the separately pinned
GATK JAR and standalone Java tools from their official sources as required. The
inventory records identities; it does not include their binaries.

CNVpytor 1.3.2 and WhatsHap 2.8 use virtual environments inheriting the `bio` Python
base; Truvari 5.5.0 uses an isolated environment. A tested recipe now obtains 40 pinned
publisher wheels and builds the other five packages from pinned published source
archives. It does **not** require the old cached wheels or a distributed binary bundle.
See the [source identities](results/runtime-wheel-sources-20260919.json),
[build-tool pins](workflow/environments/20260920/source-build-tools.json),
[completed build](results/runtime-published-source-build-20260920.json) and
[fresh verification](results/runtime-published-source-validation-20260920.json).

After installing the pinned `bio` environment and an Apple C/C++ toolchain, build into
a new directory. The measured build used CPython 3.12.14, Apple Clang 21.0.0 and macOS
SDK 26.4. The interpreter is locally pinned after the Conda installation because
relocation can change its bytes; that local hash is an integrity guard, not an
independent publisher attestation.

```sh
MVA_PY="$MVA_ROOT/tools/bio/bin/python"
MVA_BUILD="$MVA_ROOT/results/published-python-build"
MVA_PY_SHA=$("$MVA_PY" -c 'import hashlib,sys; print(hashlib.sha256(open(sys.executable,"rb").read()).hexdigest())')
"$MVA_PY" "$MVA_REPO/workflow/scripts/rebuild_published_python_wheels.py" \
  --source-receipt "$MVA_REPO/results/runtime-wheel-sources-20260919.json" \
  --build-tools "$MVA_REPO/workflow/environments/20260920/source-build-tools.json" \
  --runtime-export "$MVA_REPO/workflow/environments/20260919/runtime.json" \
  --python "$MVA_PY" --output "$MVA_BUILD" \
  --expected-source-receipt-sha256 bd7bfac5fb26c2358cc77e662cfbb2e462a832700d2a61d379f303a991ee7531 \
  --expected-build-tools-sha256 38d9e65a18d8acb6625c5e2f4ffc2ae911c23f9f01f471460f1cd0867d978369 \
  --expected-runtime-export-sha256 b508fe997179a9e8a82042f96ef0b479ed229b5cc3af104daa40572889ecffde \
  --expected-python-sha256 "$MVA_PY_SHA"
```

Require a complete `status.json` with five source builds and 45 wheels. Original
archives are preserved; prebuilt native files are omitted only from fresh source
copies. Pywfa's upstream wheel tag incorrectly suggests pure Python despite native
content, so this route deliberately requires macOS arm64 and CPython 3.12. Compiler
changes can change output hashes. Use the **newly generated** requirements, rather
than the older cached-wheel hashes in `workflow/environments/20260919/`:

```sh
for MVA_ENV in cnvpytor-1.3.2 whatshap-2.8; do
  "$MVA_PY" -m venv --copies --system-site-packages "$MVA_ROOT/tools/$MVA_ENV"
done
"$MVA_PY" -m venv "$MVA_ROOT/tools/truvari-5.5.0"
for MVA_ENV in cnvpytor-1.3.2 whatshap-2.8 truvari-5.5.0; do
  "$MVA_ROOT/tools/$MVA_ENV/bin/python" -m pip --isolated install \
    --no-index --find-links "$MVA_BUILD/wheels" --no-deps --require-hashes \
    --ignore-installed --no-compile -r "$MVA_BUILD/$MVA_ENV-local.wheels.txt"
  "$MVA_ROOT/tools/$MVA_ENV/bin/python" -m pip --isolated check
done
```

Use previously absent environment paths. The automated validation helper
`verify_python_rebuild.py` performs this installation, checks package identities and
inheritance, and executes the existing controls when given a completed
`verify_runtime_rebuild.py` result. For new source wheels, it requires
`--published-source-build --expected-wheelhouse-sha256` with the actual new build
receipt hash. The validation receipt records the full executed invocation. The
Conda verification helper uses a verified archive cache; its installation counterpart
above can populate that cache from the exact public URLs. It also needs the separately
pinned GATK JAR at the recorded relative location.

Both the earlier preserved-wheel reconstruction and this new source-build route
passed 22 selected existing controls in seven fresh environments. Six additional
native-extension cases match the old runtime. Two repeats on retained public HG002
structural caller outputs reproduce the original initial/refined counts and every
compared allele/genotype/decision record; see the
[public comparator receipt](results/runtime-published-source-public-sv-20260920.json).
These repeats did not rerun alignment or DELLY/SvABA discovery. The five new wheels
are not byte-identical to the historical builds, and full independent reproduction
of every analysis is not claimed. Exomiser/LIRICAL dependencies were inventoried,
not freshly installed in either exercise.

## 2. Reproduce the public small-variant comparison

### Historical machine-bound entrypoints

The instructions and receipts below retain the original native execution. Their
Mac executable hashes and private-layout runtime locations should not be used as
a claim of fresh-clone portability; use the new adapter above for that purpose.

The input is public HG002 NovaSeq PCR-free chr20 reads with their exact GRCh38 no-alt
reference; small-variant truth is GIAB v4.2.1 with its official confident-region BED.
The [input receipt](results/germline-public-inputs-20260919.json),
[reference preparation](results/germline-public-reference-20260919-v2.json), and
[completed native caller receipt](results/germline-public-pilot-complete-20260919.json)
record URLs, paths, hashes and scope. Primary small-variant truth-file hashes are
local integrity pins; the publisher's supplied checksum list did not cover those
three primary files. The source BAM has separate publisher-integrity evidence.

Use new output directories, preserve full chromosome context, and follow this order:

1. Prepare/download the pinned public alignment and reference using
   `prepare_public_sv_benchmark.py` and `prepare_public_assembly_reference.py`.
   These inputs are reused across the separate small-variant and structural studies.
2. Install the pinned RTG Tools 3.13 runtime and configuration described by the
   [runtime](results/rtg-runtime-20260919.json) and
   [installation](results/rtg-install-20260919.json) records.
3. Run `prepare_public_small_variant_benchmark.py` with the frozen configuration;
   verify all chr20 reference bases and the truth/BED identities.
4. Run `probe_rtg_comparator.py`; require all twelve native controls to pass.
5. Run `run_public_small_variant_pilot.py` and require completed native GATK/bcftools
   callsets. The original controller is not needed to reproduce this completed study.
6. Run `compare_public_small_variants.py` with `--root`, `--output`, `--config`,
   `--prepared`, `--controls` and `--pilot` pointing to those complete artifacts.
7. Run `audit_public_small_variants.py` for independent overlap/count diagnostics.

Each command's full argument list is recorded in the linked receipts; the executable
entry points also provide `--help`. The official BED is passed to RTG via evaluation
regions; do not preclip the VCFs. The 85,951 prepared truth records become 82,589 scored
records. Accuracy denominators are actual native output identities, not the prepared
record count. This is a single public chromosome within confident regions.

The explicit post-observation sensitivity is separately reproducible with
`probe_multiallelic_representation.py`, `prepare_public_multiallelic_sensitivity.py`,
the unchanged comparator, `summarize_public_representation.py`, and
`review_public_representation_losses.py`. It retains the original raw multiallelic
genotypes. GATK gains 219 truth matches and loses one; bcftools gains 204 and loses
four. [All losses](results/germline-representation-loss-review-20260919.json) remain
available alongside the gains. The frozen primary result is not replaced.

## 3. Reproduce the authorised case analysis

Start with the [case-evidence and reconstruction map](reports/appendices/track1-case-evidence.md).
It connects the reported ranks, predictor sensitivities, selected allele support,
locus/ten-gene census, coding depth and phase to the corresponding inputs and
entry points. The supplied-VCF ranking runs and final four-lane native-call
reviews are separate analyses. The map states which original outputs remain
restricted or unavailable in this checkout; it is not a fresh-clone case replay.

Obtain the challenge data through its access process and retain it outside the Git
checkout. The internal analysis used all eight source FASTQs from four lanes, with
publisher-checksum verification. Do not create an inferred clinical dataset or treat
the public HG002 reads as patient data.

- `01_align.sh` performs reference-aware BWA alignment and duplicate marking,
  preserves separate lane artifacts and validates the completed CRAM.
- `call_genome.py` runs native diploid GATK and bcftools discovery by whole-reference
  contig shards. `resume_reads.py` and `resume_calling.py` are dependency-aware
  controllers for the original setup; inspect their configuration before reuse.
- GATK 4.6.1.0 / HTSJDK 4.1.3 cannot consume CRAM 3.1. For that input format,
  `prepare_gatk_compatible_cram.py` creates a separate CRAM 3.0 copy, preserves
  the original, verifies every decoded record, checks headers/index/reference,
  and runs an indexed GATK read probe. Supply the actual source/reference SHA-256
  values, expected total alignment-record count, sample and a probe region using
  its CLI options. Optional `--call-output` requires the expected unchanged caller
  SHA-256 and a new directory. Do not reuse the original failed output or restart
  its old waiting controller. The helper compiles `verify_alignment_records.cpp`
  with the installed Apple C++ compiler and pinned HTSlib headers/library.
- `prioritize.py` and the Exomiser configurations perform the phenotype-assisted
  analyses. `run_lirical.py` and related summary/resource helpers provide the separate
  disease comparison. `rare_biallelic_baseline.py` and transcript helpers implement
  the phenotype-independent review. Exact resource versions and paths must be supplied.
- Read, mate and WhatsHap helpers retain source genotype/phase context. Full original
  allele lists, GT, AD and PL are needed alongside split allele rows; a matching
  component must not hide an incompatible second allele or imply parental phase.
- `audit_final_alignment_targets.py` verifies the completed four-lane receipt,
  extracts bounded SNP windows without changing whole-alignment duplicate marks,
  and uses the existing explicit-reference and fragment evidence helpers. The
  separate four-lane phasing configuration declares its exact regional scope.
- The CNVpytor, structural and ROH workflows require final alignment/QC, appropriate
  masks, ploidy/callability review and their own interpretation. An ordinary diploid
  callset does not establish mosaic, mitochondrial, sex-chromosome or structural coverage.
- `review_four_lane_roh.py` reconciles the completed panel/ROH outputs, original
  genotypes, exclusions, segments and model summaries without changing the analysis.
  Its descriptive length bins are not clinical cutoffs or evidence of parental origin.
- `review_four_lane_copy_number.py` checks all native CNV fields and stored read
  counts, compares reference-gap metrics to the FASTA, and annotates every call
  against the pinned Ensembl GFF3. Its separately configured post-run review subsets
  retain all rejected records and do not classify dosage pathogenicity.
- `fetch_cnv_context_resources.py --output <new-source-directory>` obtains the
  generation-pinned gnomAD v4.1 genome-SV VCF/index and current ClinGen GRCh38 gene/
  region TSVs. The recorded public-source receipt pins the dated files used here;
  ClinGen daily URLs may later change, so preserve and reuse the verified snapshot
  for exact dated reproduction. gnomAD publisher MD5, generation, size, gzip
  integrity and local SHA-256 are checked; CRC32C is recorded without rechecking.
- `review_cnv_context.py` takes `--cnv-review`, `--roh-analysis`, `--resources`,
  `--config workflow/configs/cnv-context.json`, a new `--output` directory, and
  the corresponding `--expected-cnv-review-sha256`, `--expected-roh-summary-sha256`
  and `--expected-resources-sha256`. It verifies shared alignment/reference identities,
  retains all raw CNV records and panel markers, and scans the full population VCF.
  Nominal gnomAD intervals follow the official BED convention (VCF POS..END inclusive).
  Simple PASS overlap, multiallelic copy states, filtering and dosage curation remain
  distinct; missing matches are not assigned a zero population frequency.
- `summarize_cnv_context.py --context <completed-context> --expected-summary-sha256
  <hash> --output <new-json>` checks the upstream artifacts and all overlapping
  multiallelic sample count/frequency denominators. `plot_cnv_context.py` can plot
  explicitly selected calls using the same context plus `--analysis`, `--call ID=label`
  and a new output directory. The plot uses the recorded CNVpytor environment with
  NumPy/Matplotlib. Patient-derived context tables, plots and point data stay outside
  this source repository. The public-source and nine-control receipts are included.
- If DELLY fails on the original compressed-read decoder path, the existing
  fully record-verified CRAM 3.0 copy can be used through
  `resume_structural_from_verified_copy.py`. This pins the completed census,
  conversion proof, native synthetic format control and unchanged policy, and
  preserves the original failed run. `probe_structural_cram30.py` supplies the
  bounded synthetic format/caller check; its imported fixture parser requires NumPy
  (the recorded control uses the CNVpytor Python environment). See the
  [completed control](results/structural-cram30-control-20260920.json).
- `review_four_lane_structural.py --analysis <completed-structural-directory>
  --cnv-review <completed-CNV-review> --gff <pinned-GFF> --output <new-directory>`
  takes `--expected-summary-sha256`, `--expected-cnv-review-sha256` and
  `--expected-gff-sha256`. It verifies all native outputs and reuses the original
  decoder/serializer to reconcile the complete call census byte for byte. It
  compares all primary CNV records with same-type nominal structural spans,
  retains reference/missing/filtered/sex/query-boundary exclusions, and describes
  selected-gene endpoints against the CDS union across annotated transcripts.
  This is not independent decoder validation, complex-allele reconstruction,
  transcript consequence prediction or clinical classification.
- `review_tandem_junction.py` performs a bounded follow-up of a selected precise
  `5to3` DUP without an inserted sequence. It takes `--analysis`, `--record-id`,
  `--expected-summary-sha256`, `--gff`, `--expected-gff-sha256`, `--config`,
  `--samtools` and a new `--output` outside this repository. The configuration is
  `workflow/configs/tandem-junction-review.json`. It rechecks the complete input
  alignment/reference/index identities, retains every exact local consensus
  decomposition, and stops for manual review if the decomposition is ambiguous.
  Unfiltered regional extraction preserves the original duplicate marks and all
  alignment records. Exact sequence, observed split segments and outward pairs
  are reported as overlapping fragment sets. Transcript geometry is conditional
  on a simple tandem arrangement, not proof of expression or gene disruption.
  Sequence/read tables and patient results remain outside this source export.
- `review_paralog_candidate_pair.py --root <ROOT> --config <JSON> --output <new-directory>`
  checks a nearby candidate pair against the completed raw small-variant shard and
  selected structural calls. The example configuration is
  `workflow/configs/fancd2-paralog-review.json`; its identities are case-specific and
  must refer to the intended completed inputs. The script rechecks input hashes,
  retains raw regional calls/reads, and counts RG/QNAME fragments at several mapping
  thresholds. It remaps original primary reads, actual mates and reference windows
  individually with the pinned full-reference BWA index, retaining reported alternatives.
  `summarize_paralog_mates.py --analysis <completed-directory>
  --expected-summary-sha256 <SHA256> --output <new-directory>` describes companion
  support from those saved mappings, including reciprocal coordinates and clipping.
  This sensitivity removes pairing information; weak single-end mapping does not
  establish an artifact, and a favored mapping does not demonstrate unique origin.
  Read and mate observations are not independent molecules. Case results and sequences
  remain outside the source export.

At this source snapshot the four-lane alignment and format verification are complete;
small-variant whole-reference calling and broader interpretation remain unfinished.
The structural run and initial CNV/structural reconciliation are complete. The patient-derived ranking/read summaries,
draft submissions, final reports and video are not included here. They remain separate
work products awaiting their stated review; this source repository does not complete
those submission requirements.

## 4. Check components and preserve evidence

The mitochondrial workflow is a native adaptation of the pinned Broad regular/
shifted-reference method, not a completed contamination or clinical-validation
pipeline. Acquire the checksum-pinned public FASTAs, chain, intervals and mask:

```sh
"$MVA_ROOT/tools/bio/bin/python" "$MVA_REPO/workflow/scripts/prepare_mitochondrial_resources.py" \
  --root "$MVA_ROOT" --sources "$MVA_REPO/workflow/configs/mitochondrial-public-resources.json" \
  --output "$MVA_ROOT/ref/mitochondrial-20260920"
"$MVA_ROOT/tools/bio/bin/python" -m unittest discover \
  -s "$MVA_REPO/workflow/tests" -p test_mitochondrial_discovery.py -v
MVA_BCFTOOLS="$MVA_ROOT/tools/bio/bin/bcftools" "$MVA_ROOT/tools/bio/bin/python" \
  -m unittest discover -s "$MVA_REPO/workflow/tests" -p test_mitochondrial_review.py -v
"$MVA_ROOT/tools/bio/bin/python" "$MVA_REPO/workflow/scripts/probe_mitochondrial_discovery.py" \
  --root "$MVA_ROOT" --resources "$MVA_ROOT/ref/mitochondrial-20260920" \
  --output "$MVA_ROOT/results/mitochondrial-control-reproduction"
```

Use new output directories. Require seven method checks, two native normalization
checks and all four planted alleles to pass. Preserved pre-repair sources record
read-group escaping and unmapped-BAM dictionary handling failures; use the current
scripts. Synthetic mixtures test fractions 0.10, 0.20, 0.25 and 1.00, not a validated
1% detection limit, sequencing-error model or NUMT discrimination.

`run_mitochondrial_discovery.py` requires `--root`, `--reference`, `--alignment`,
`--resources`, `--config` and `--output`. Use the verified compatible CRAM and
`workflow/configs/mitochondrial-discovery.json`. It verifies the source/reference
hashes and mitochondrial sequence, admits the upstream read subset, compares all
primary forward sequences/qualities and preserves separate mapping/duplicate decisions.
The calling partition is regular positions 576–16024 plus shifted positions
8025–9144, covering each original position once. Preserve lifted rejects even when
empty, both raw VCFs, original multi-allele context and full depth missingness.

Initial output filenames explicitly state that filtering assumes no contamination.
The upstream contamination component and independent nuclear-copy checks remain
unimplemented. `review_mitochondrial_calls.py` compares a completed `--run` with
`--vendor-vcf` and `--clinvar-vcf`, requiring their pinned identities, plus `--root`
and a new `--output`. It records the original ALT index and full genotype/filter
context; a split row never replaces the original genotype. `audit_mitochondrial_targets.py`
uses `--root`, `--run`, `--config` and `--output` for the bounded original/remapped
read check. Its configuration is selected after call review, not an independent
performance benchmark. Unresolved clinical annotations, low-quality calls or a
negative contamination screen must not be promoted to biological truth.

`review_mitochondrial_read_context.py` accepts `--root`, `--config` and a new
external `--output` directory. The context configuration pins a completed case
receipt and its prior target audit. It checks original primary reads and reciprocal
mates, reconciles earlier BAQ fragment counts, and reports whole-reference
single-end BWA alternatives with explicit target-base projection. Detailed read
identities, sequences and alignments stay outside this source snapshot. Run the
seven geometry, missingness, fragment and symmetric-filter controls first:

```sh
"$MVA_ROOT/tools/bio/bin/python" -m unittest discover \
  -s "$MVA_REPO/workflow/tests" -p test_mitochondrial_read_context.py -v
```

Read-end exclusion, base-quality thresholds and absence of additional aligned
edits are post hoc sensitivity analyses. They are not calibrated heteroplasmy
rules. Reference differences may be genuine linked alleles, and BWA's reported
alternatives do not exhaustively exclude NUMTs or contamination.

The chromosome-depth survey freezes intervals from reference coordinates before
read inspection, requires explicit depth for every selected position and preserves
zero-depth windows. Reference eligibility is based only on concrete-base fraction
and GC composition. Reproduce it with the completed four-lane state:

```sh
"$MVA_ROOT/tools/bio/bin/python" "$MVA_REPO/workflow/scripts/survey_chromosome_depth.py" \
  --root "$MVA_ROOT" --alignment-state "$MVA_ROOT/state/raw-reads.json" \
  --config "$MVA_REPO/workflow/configs/chromosome-depth-survey.json" \
  --output "$MVA_ROOT/results/chromosome-depth-survey-reproduction"
```

The separate `call_ploidy_regions.py` runner needs `--root`, `--reference`,
`--alignment`, `--config` and a new `--output` directory. Its frozen model is
`workflow/configs/sex-chromosome-ploidy-model.json`; it verifies exact reference and
alignment hashes, read groups, reference lengths and the complete Y PAR N masks.
The model assigns diploid X PAR and haploid other X/Y sequence without changing
subject metadata. Use the verified compatible CRAM, not an unvalidated format
conversion. Preserve raw multi-allele VCFs and later reconcile all boundary-crossing
records and caller discordances. Variant-only absence is not reference evidence.
GATK uses its native default filters; bcftools pileup uses MQ20/BQ20, excluded flags
0xF04 and a 10,000 depth cap. These are raw calls without BQSR or clinical filtering.

GATK 4.6.1.0 explicitly selects the greatest ploidy among intervals overlapping an
active region. A variant near a boundary can therefore have a genotype shape that
differs from the model at its start. The strict regional runner rejects that
condition after retaining both original caller outputs. This behavior is documented
in the [pinned HaplotypeCaller source](https://github.com/broadinstitute/gatk/blob/4.6.1.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/HaplotypeCaller.java).

`call_fixed_ploidy_regions.py` uses the same five required path arguments and adds
`--bcftools-vcf <COMPLETED_ORIGINAL_VCF>` and `--bcftools-sha256 <EXPECTED_SHA256>`.
It independently calls the complete selected intervals once per constant ploidy,
without regional overrides, and reuses the identified bcftools file unchanged.
Per-ploidy raw VCFs retain every emitted record. Combined start ownership is assigned
before normalization; outside-owner records are logged, and REF spans crossing a
boundary remain explicit review items. Original ALT order, GT, AD and PL are
preserved. No genotype is coerced to fit a model.

Run the boundary-specific controls on the pinned runtime:

```sh
"$MVA_ROOT/tools/bio/bin/python" -m unittest discover \
  -s "$MVA_REPO/workflow/tests" -p test_fixed_ploidy_regions.py -v
"$MVA_ROOT/tools/bio/bin/python" "$MVA_REPO/workflow/scripts/probe_ploidy_boundaries.py" \
  --root "$MVA_ROOT" --output "$MVA_ROOT/results/ploidy-boundary-control-reproduction"
```

The native probe uses a seeded 12 kb fictional genome, 600 reads and 12 SNVs at
three ploidy transitions, including each interval's last/first bases. It requires
the original regional runner to finish both callers and then fail its strict
boundary check. Separate fixed-ploidy calls must recover all 12 expected genotypes,
as must the reused bcftools output: 24 caller/allele checks. The four focused tests
also exercise invalid/partial genotypes, missing likelihoods, original multiallelic
index preservation and retention of a constructed boundary-spanning deletion.
These controls do not establish sensitivity in noisy or clinical data.

Run the four depth and three interval controls with the pinned native runtime:

```sh
MVA_SAMTOOLS="$MVA_ROOT/tools/bio/bin/samtools" "$MVA_ROOT/tools/bio/bin/python" \
  -m unittest discover -s "$MVA_REPO/workflow/tests" -p test_chromosome_depth_survey.py -v
"$MVA_ROOT/tools/bio/bin/python" -m unittest discover \
  -s "$MVA_REPO/workflow/tests" -p test_ploidy_regions.py -v
"$MVA_ROOT/tools/bio/bin/python" "$MVA_REPO/workflow/scripts/probe_ploidy_regions.py" \
  --root "$MVA_ROOT" --output "$MVA_ROOT/results/ploidy-native-control-reproduction"
```

The final command creates fictional sequence/reads and exercises both actual callers
on three alleles spanning diploid and haploid intervals. Require all six caller/allele
genotype comparisons to pass. The preserved pre-repair source documents an index
collision fixed before the passing run; it is not the current runner. These controls
provide software-contract evidence, not a clinical sensitivity estimate. The survey
is not gene-level callability, and mitochondrial circular-boundary, NUMT, heteroplasmy
and contamination analyses remain separate work.

The two selected baseline indels use `audit_final_alignment_indels.py --root <ROOT>
--alignment-state <ROOT>/state/raw-reads.json --config
workflow/configs/baseline-four-lane-indels.json --output <NEW_DIRECTORY>`. This preserves
source records and full multi-allele caller genotypes, then reports exact sequence windows
with three flank sizes and three MQ thresholds. `diagnose_indel_windows.py --run
<COMPLETED_DIRECTORY> --reference <ROOT>/ref/realignment/GRCh38_nochr.fa
--output <NEW_JSON>` explains missing boundaries and low-quality windows without changing
filters. Run the four `test_four_lane_indels.py` and two `test_indel_window_diagnostics.py`
controls in the pinned native Python environment. These are bounded software controls;
zero qualifying windows or a zero-depth reference call is not an indel exclusion.

The baseline SNV follow-up uses the existing `audit_final_alignment_targets.py` with
`workflow/configs/baseline-four-lane-snvs.json`, a completed four-lane `state/raw-reads.json`,
and a new external output directory. `review_snv_cluster_mapping.py --root <ROOT> --config
workflow/configs/serpina1-cluster-review.json --output <NEW_DIRECTORY>` then retains selected
aligned-base configurations, reciprocal mates and full-reference BWA alternatives.

For the independent public-sequence comparisons, obtain the EMBL records
[AL117259.6](https://www.ebi.ac.uk/ena/browser/api/embl/AL117259?download=false) and
[JX680599.1](https://www.ebi.ac.uk/ena/browser/api/embl/JX680599?download=false), retain their
archive metadata and verify the complete-file SHA-256 values specified in
`serpina1-public-clone-comparison.json` and `serpina1-public-transcript-comparison.json`.
Store them at each configuration's `public_embl` path under the analysis root. The second
record is a spliced mRNA/CDS sequence; its alignments and spans cannot be interpreted as
genomic insert lengths, expression or a physical haplotype.

Run `compare_cluster_public_sequence.py --root <ROOT> --config <COMPARISON_CONFIG>
--output <NEW_DIRECTORY>` separately for each public record. Then run
`project_cluster_target_bases.py --root <ROOT> --comparison <TRANSCRIPT_COMPARISON_DIRECTORY>
--output <NEW_JSON>` to check that the original target bases themselves participate in the
matches, with both read orientations accounted for. Small-reference mapping quality is
not genome-wide confidence; missing or unmapped score comparisons remain unavailable.

Example configurations pin the original completed receipts. A fresh reproduction creates
new timestamped receipts: first validate its inputs, commands and complete outputs, retain
the original and new receipt identities, and set the dependent configuration's receipt
pin to that validated local result. Do not bypass a mismatch or silently reuse old results.
Five cluster, four public-comparison and two projection controls ran on 512S2 with source
identities recorded in the selected receipts. Case reads and findings are outside this export.

`plan_assay_capacity.py --config workflow/configs/assay-pilot-capacity.json
--output <new-output.json>` uses standard-library arithmetic to calculate a prospective
well/imaging/storage scenario. The example retains independent culture preparations when
counting minimum instrument sessions and includes both temporal endpoints. It assumes no
compression savings or image-processing overhead; the result is not measured laboratory
throughput, statistical power, a drug exposure range or an observed experiment. The clinical
pilot protocol remains a separate reviewable work product outside this source export.

The `workflow/tests` directory contains unit and native integration checks. Some
require external runtimes or previously completed public-control artifacts. A skipped
test is not a passing native test. Run the tests relevant to the installed component,
inspect skips and native outputs, and keep the executed source/configuration hashes.

The nine self-contained CNV context controls need only Python's standard library:

```sh
python3 -m unittest discover -s "$MVA_REPO/workflow/tests" -p test_cnv_context.py -v
```

The two structural/CNV selection controls also use the standard library:

```sh
python3 -m unittest discover -s "$MVA_REPO/workflow/tests" -p test_structural_cnv_context.py -v
```

The tandem-junction suite has five standard-library controls and one native BAM
control requiring pysam. All six ran on the Studio; a missing-pysam skip is not
a passing native control:

```sh
python3 -m unittest discover -s "$MVA_REPO/workflow/tests" -p test_tandem_junction.py -v
```

The six nearby-pair/mate controls require pysam and include mapping-quality 255,
fragment conflicts, tied/absent alternative mappings, reciprocal mates and clipping:

```sh
"$MVA_ROOT/tools/bio/bin/python" -m unittest discover \
  -s "$MVA_REPO/workflow/tests" -p test_paralog_candidate_pair.py -v
```

They test representation and missingness rules, including checksum corruption and
different allele/sample denominators. They are not a clinical performance benchmark.

The self-contained synthetic format control needs the prepared native runtime but
no challenge data. From any reviewed source copy:

```sh
MVA_PROJECT_ROOT="$MVA_ROOT" "$MVA_ROOT/tools/bio/bin/python" -m unittest discover \
  -s "$MVA_REPO/workflow/tests" -p test_gatk_cram_compatibility.py -v
```

Require the native test to run and pass. Its default temporary fixtures contain
only generated sequence. It exercises the exported recovery code and both real
callers, and is not a whole-genome clinical benchmark.

The source snapshot copy is byte-checked. Its only alignment-script change is a
documented comment replacement removing an unnecessary acquisition-specific header
and barcode example; executable lines are identical. Python compilation and shell
syntax checks validate the copy's syntax, not end-to-end scientific correctness.
Controlled sequencing data, indexes and genotype-scale intermediates are not inputs
to this Git repository and must not be committed to make tests convenient.

## Mitochondrial mixture screen

The [runtime pins](workflow/configs/mitochondrial-contamination-runtime.json)
identify the official haplocheck archive, bundled mutserve v2.0.0-rc12,
haplocheck 1.3.3 and rCRS reference. Publisher-provided SHA-256 metadata was absent;
the recorded hashes identify the downloaded bytes. Runtime binaries are acquired
separately. The [exact native control sources and setup](workflow/source-snapshots/mitochondrial-contamination-20260920/README.md)
document the fixed study directories and reproduction order. The qualification
receipt and its unchanged source/control files are required before a case run.

After acquiring and verifying the native runtimes, completing the controls, and
setting `MVA_REFERENCE` and `MVA_ALIGNMENT` to the exact inputs in the reviewed
configuration, use a **new output directory with no periods in its entire path**:

```sh
"$MVA_ROOT/tools/bio/bin/python" "$MVA_REPO/workflow/scripts/run_mitochondrial_contamination.py" \
  --root "$MVA_ROOT" --reference "$MVA_REFERENCE" --alignment "$MVA_ALIGNMENT" \
  --config "$MVA_REPO/workflow/configs/mitochondrial-discovery.json" \
  --qualification "$MVA_ROOT/results/mutserve-contract-controls-20260920/qualification-verification.json" \
  --parser-probe "$MVA_REPO/workflow/source-snapshots/mitochondrial-contamination-20260920/HaplocheckParserProbe.java" \
  --output "$MVA_ROOT/results/mitochondrial-contamination-NEW"
```

The wrapper preserves original whole-genome mappings, excludes unsupported read
flags and unknown MQ/AS, and invokes native MQ30/BQ20/AS30 filtering. Overlapping
mates remain separate observations. It reconciles the full caller text table with
the VCF, retains qualified biallelic SNP records without GT/AF recoding, verifies
their interpretation through the actual haplocheck importer, and checks fresh
text/raw/HTML reports in hotspot-inclusive and hotspot-excluded arms. A three-allele
control shows that the caller's VCF can drop one allele; an additional native
rank-switch control shows that changing the rank of another minor allele can
suppress an otherwise supported allele. These are explicit limitations, not
clinical validation. NO/ND cannot be reported as zero or nuclear contamination.

[Native contract results](results/mutserve-contract-controls-20260920.json),
[source-linked controls](results/mitochondrial-mixture-method-controls-20260920.json).

## Local narrated pitch production

The [second pitch-template recipe](workflow/docs/pitch-v2-reproduction.md) describes
the four native/Python production and inspection programs. It records the tested
macOS 26.5.1 environment separately from the genomic runtime above. Supply a reviewed,
authorized storyboard outside this source snapshot and produce a new output tree.
Media checks cover measured duration, captions, decoded speech windows and every
planned static-segment midpoint; they do not establish scientific validity, complete
pronunciation review or upload status. Earlier pitch source remains preserved.

## Read evidence for changed original calls

The [changed-call recipe](workflow/docs/changed-call-read-reproduction.md) reuses the
qualified BAQ/fragment and exact-window helpers without a new genotyper run. It
retains native haploid/multiallelic caller context and rejects unknown mapping
quality before using the SNP helper. [Fictional controls](workflow/source-snapshots/changed-call-read-controls-20260920/controls.json)
and the [observed runtime](workflow/source-snapshots/changed-call-read-controls-20260920/runtime.json)
record the executed qualification. Overlapping hypotheses, windows and quality
thresholds are not independent samples and must not be added. Read support and
annotation proximity do not establish locus assignment, pathogenicity or causality.

## Completed-autosome handoff check

Follow [the handoff reproduction instructions](workflow/docs/discovery-handoff-reproduction.md)
using the verified native environment and a new output directory. The provided
calling-definition hash and original reference must match. The audit snapshots only
completed autosomes and requires every chromosome when `--require-complete` is set.
It exercises indexes through explicit contig queries, replays split normalization,
and retains links to newly fingerprinted original multiallelic files. Counts under
the older quality screen diagnose compatibility; they do not adopt new filtering
thresholds or infer missing GQ. No native caller or completed output is modified.

Five [native fictional controls](workflow/source-snapshots/discovery-handoff-controls-20260920/controls.json)
passed on Python 3.12.14, pysam 0.24.1 and bcftools/htslib 1.24 on a Mac Studio.
They test preservation, missingness, raw-genotype tampering, wrong indexes, incomplete
scope and empty valid input. These controls concern handoff consistency, not
biological sensitivity. The full case audit and case callsets remain outside this
source-only export.

## Native coding review after the handoff

Use the [coding review recipe](workflow/docs/NATIVE-CODING-REVIEW.md) with a
checksum-verified completed handoff, the same reference/index, a pinned GFF, the
review configuration and a new caller-specific output directory. A partial
autosomal snapshot requires the explicit flag and remains incomplete. The
workflow creates derived ordinal-tagged and normalized records while retaining
full original GT/ALT/AD/PL context; annotation may add consequences but may not
drop records or change non-consequence fields.

Quality fields are descriptive. Population counts are joined by exact normalized
allele; missing observations remain unresolved. Model generation keeps callers
and source phase blocks separate, quarantines duplicate keys, and retains cis or
overlap combinations for representation review. These outputs are descriptive
combinations, not clinical findings or EPCR estimates.

Eighteen [fictional helper and native integration checks](workflow/source-snapshots/native-coding-controls-20260920/controls.json)
passed on 512S2 with Python 3.12.14 and pysam 0.24.1. Their frequency inputs are
fictional exact-response files; they do not validate the public frequency database
or measure patient accuracy. The case execution uses the pinned local public
database and requires separate output reconciliation.

## Joint sequence and fragment review

Follow the [joint-haplotype recipe](workflow/docs/JOINT-HAPLOTYPE-REVIEW.md) with
the exact source/configuration pins, authorized external inputs and a new output
directory. The workflow preserves native caller context, constructs each
caller's compatible local sequence hypotheses separately and retains rejected
overlapping edit combinations. Sequence-identical hypotheses retain their
individual source labels without multiplying their observed fragment support.

Reproduce the [20 fictional controls](workflow/source-snapshots/joint-haplotype-controls-20260920/controls.json)
on the native Studio environment before analysing authorized reads:

```sh
MVA_PROJECT_ROOT="$MVA_ROOT" "$MVA_ROOT/tools/bio/bin/python" -m unittest discover \
  -s "$MVA_REPO/workflow/tests" -p 'test_joint_haplotype*.py' -v
```

Require all tests to run, with no missing-runtime skips. The recorded run used
Python 3.12.14 and pysam 0.24.1. It exercises native CRAM/raw-call extraction and
local annotation alongside the sequence/fragment helper controls. Fictional
inputs establish the exercised software contract, not biological sensitivity.

Exact read-window comparisons use raw base qualities, not BAQ. Each threshold
and window retains its own observations, read-group/strand axes, fragment
identities, unmatched sequences and conflicts. Nested windows, quality views
and alternative source lineages reuse evidence and must not be added as
independent observations. Constructed local consequences do not establish a
patient genotype, parental phase, locus assignment or functional rescue.

The private configuration preserves target/sample/input-hash selection pins.
Acquire those authorized external inputs separately; do not commit individual
read records, read names, case evidence/result payloads, clinical reports or
media to reproduce this source-only snapshot.

## Full-query fitting and whole-reference remapping

Use the [clipped-read recipe](workflow/docs/CLIPPED-READ-REVIEW.md) with the
pinned native Studio runtime and authorized external inputs. The private
configuration retains the selected sample, intervals, input identities and
resource limits. The runner accepts `--root`, `--config` and a new `--output`
directory. Original alignment/seed indexes, the reference FAI, both complete
raw caller files and the BWA index receipt and files must match their pins.
Original data, indexes and prior outputs are never rewritten.

The [qualification receipt](workflow/source-snapshots/clipped-read-controls-20260920/controls.json)
records 12 tests passed with no skips at `2026-09-20T15:21:44.485903+00:00` on
Mac Studio 512S2: seven helper controls and five native integration controls,
using Python 3.12.14 and pysam 0.24.1. It pins five workflow scripts, the
configuration and two test files. Reproduce the complete suite before using a
revision; reject skipped, empty or incomplete execution:

```sh
MVA_PROJECT_ROOT="$MVA_ROOT" MVA_TEST_DIRECTORY="$MVA_REPO/workflow/tests"   "$MVA_ROOT/tools/bio/bin/python" - <<'PYCONTROLS'
import os
import unittest
suite = unittest.defaultTestLoader.discover(
    os.environ['MVA_TEST_DIRECTORY'], pattern='test_clipped_read*.py')
result = unittest.TextTestRunner(verbosity=2).run(suite)
if not result.wasSuccessful() or result.skipped or result.testsRun != 12:
    raise SystemExit(1)
PYCONTROLS
```

The native controls use generated sequence and reads, real BAM/CRAM/VCF/BCF
files, and actual BWA indexing and remapping. They do not run native genotypers
or measure sensitivity on clinical data. The retained receipt identifies the
qualified bytes; source, configuration or runtime changes need their own
qualification evidence.

Local comparison consumes every available query base with unit edit costs.
It fits eligible primary seed ends only; recovered mate-only ends receive
original-context review and genome remapping. Missing qualities permit an
explicitly unweighted local fit when the sequence is concrete, but exclude
FASTQ export. Multiple optimal placements and capped witnesses must remain
explicit. Target offsets refer to constructed sequence; indel-containing
panels require edit-coordinate mapping before genomic interpretation.

The single and paired remapping arms preserve every reported BWA record and
check full exported sequence/quality payloads where present. Read ends,
identifier-group fragments, source lineages and reused remap arms are distinct
denominators. The workflow does not assign allele-support fragment counts.
Each actual run requires separate status/summary, fresh hash, end/export,
fit/CIGAR and fragment reconciliation. Keep those private read-level outputs
outside this source snapshot; fictional qualification alone does not establish
case correctness, locus assignment, parental phase or clinical validity.

## Prospective founder-cell ledger

Follow the [founder-ledger schema and execution recipe](workflow/docs/FOUNDER-LEDGER.md).
Provide four strict JSON inputs: the frozen design, planned-well registry,
starting-cell roster and adjudicated outcomes. Supply a separately recorded SHA-256
for every input; the design must also pin the roster. Use a new output directory.
The runner requires only the Python standard library and records its own source
identity along with every input identity.

The [qualification receipt](workflow/source-snapshots/founder-ledger-controls-20260920/controls.json)
pins the source and test files used for fictional controls on the Mac Studio.
Reproduce that complete suite and compare its test count to the receipt before
using the recorded revision:

```sh
python3 -m unittest discover -s "$MVA_REPO/workflow/tests" -p 'test_founder_ledger.py' -v
```

Keep zero-denominator fractions unavailable, retain unscorable anaphases separately
from tracking loss, and keep failed/unrun wells in the planned registry. A known
starting census from an interrupted acquisition still contributes its adjudicated
outcomes. An unknown starting census does not become a valid empty well. Cells
already mitotic at time zero retain their left-censored entry rather than an
invented entry time. First-event classifications cannot be overwritten by later
secondary outcomes supplied to the ledger.

The output validates the supplied adjudication contract. It does not reconstruct
images, detect an earlier event omitted from every input, prove independent
biological replicates, infer drug efficacy or verify that prospective freezing
actually occurred. Experimental observations and case-specific outputs belong
outside this source-only repository.


## Coding-region depth method

The source snapshot now includes the fixed-gene CDS-depth method, its exact
configuration and ten native fictional controls. The method retains explicit
zero coverage and reference denominators, refuses unavailable mapping/base
qualities, and reports three declared quality views. Depth alone does not call
genotypes or establish clinical callability. The configuration contains the
authorised sample/input identities and selected gene list; alignment data and
case results are excluded from this source export.

See the [reproduction guide](workflow/docs/CODING-DEPTH-REPRODUCTION.md) and
[qualified fictional-control record](workflow/source-snapshots/coding-depth-controls-20260920/controls.json).


## Native overlay and secondary-inventory additions

Use the [overlay guide](workflow/docs/ADAMTS1-NORMALIZED-OVERLAY.md) and the
[secondary reproduction guide](workflow/docs/NATIVE-SECONDARY-REPRODUCTION.md).
The [secondary contract](workflow/docs/NATIVE-SECONDARY-CONTRACT.md) remains an
unchanged historical preparation document; the reproduction guide records the
final source layout and prerequisites. These additions contain no actual case
observations or clinical reporting decisions.

Four original native fictional-control receipts are retained unchanged. Their
source paths map to this release as follows:

| Original receipt | Controls | Receipt source-name mapping |
|---|---:|---|
| [Overlay](workflow/source-snapshots/adamts1-normalized-overlay-20260920/controls.json) | 21 | All five `source_sha256` names resolve below `workflow/source-snapshots/adamts1-normalized-overlay-20260920/`. |
| [Independent overlay checker](workflow/source-snapshots/adamts1-normalized-overlay-20260920/independent-verification/controls.json) | 23 | Both source names resolve beside that receipt. |
| [Secondary producer](workflow/source-snapshots/native-secondary-controls-20260920/controls.json) | 35 | The selected historical names map under `workflow/`: `docs/NATIVE-SECONDARY-CONTRACT.md`, `scripts/qualify_native_secondary.py`, `scripts/review_native_secondary.py`, both `scripts/vendor/` dependencies and `tests/test_native_secondary.py`. |
| [Independent secondary checker](workflow/source-snapshots/native-secondary-independent-20260920/controls.json) | 36 | All three `source_pins` names resolve below its adjacent `source/` directory. |

The secondary producer's original receipt also pins seven historical files not
exported here: three earlier configuration drafts, `docs/README.md`,
`scripts/launch_secondary.py`, `scripts/preflight_public_secondary.py` and
`tests/test_launch_secondary.py`. No full historical workflow snapshot is claimed.
The approved version 2 configuration is included as separately reviewed exact
execution configuration; it is not relabeled as one of those historical draft
configuration bytes. Reproducing an independent qualification requires the entire
adjacent producer workflow snapshot specified by the guide, or a newly qualified
relocated snapshot and its new receipt. Copying a historical receipt alone does
not satisfy this binding.

The original 21/23/35/36 control results are reused because their exported source
bytes are unchanged. They were not rerun for this source export. Configuration
sample identifiers, targets, hashes and source-selection counts are intentionally
included under the private-source policy; original alignments, variant payloads,
case observations, clinical reports and media remain separate authorized inputs.
Neither fictional controls nor source publication establishes phase, locus origin,
biological function, complete panel callability, reportability or causality.


## v26 source and native-layout evidence

The approved 62-file source/configuration/fictional-control selection is preserved without rewriting qualified methods. Nine original suites retain 6, 16, 21, 40, 8, 25, 25, 34 and 18 controls respectively, with all 41 source-field entries mapped. These are separate historical tests, not a pooled validation score. The 18 population controls are Java self-tests; no separate guard-test suite is implied. No scientific controls were rerun for this release preparation.

Read [the dependency and replay boundaries](workflow/docs/V26-SOURCE-REPRODUCTION.md), [exact restoration paths and hashes](workflow/docs/V26-NATIVE-RESTORATION.json), and [receipt-to-source mappings](workflow/docs/V26-QUALIFICATION-MAP.json). The four unchanged depth sources and original qualification already in v25 are explicitly mapped to the historical native layout. Approved configurations remain distinct from qualified draft configurations; both are retained where the launcher checks their authorization delta.

Complete source-field mapping does not include every receipt-referenced runtime, generated fixture, log, public resource, private input or case result. In particular the population wrapper still needs omitted case-selected query/prior-response files even in its historical qualification mode, while the Java self-test needs its recorded runtime only. Native path/hash/output guards remain intact. No fresh-clone patient execution, biological sensitivity, clinical callability or sample-purity claim follows from this source package.


## v27 source and replay boundary

See the [v27 reproduction notes](workflow/docs/V27-SOURCE-REPRODUCTION.md) and [source map](workflow/docs/V27-SOURCE-MAP.json). Exact historical source/receipt bindings are preserved separately from generic release adaptations and their local checks. Saved-output parity is not new biological evidence. Replaying the historical research still requires authorized private inputs, pinned public resources, runtime identities and preserved native layouts. Speech synthesis additionally requires a user-supplied provider credential, permitted credit use and a narration manifest; no speech requests are made by the archive preparation or verification.


## Public structural-variant residual annotation sources

Follow the [source snapshot instructions](workflow/source-snapshots/public-sv-residual-review-20260920/README.md) for exact source mappings, separate control scopes, the reused guard-helper dependency and the twelve saved public input requirements. Original source bytes and historical plans are preserved. The plans bind a particular completed run; reproduce deliberately with a new owned output directory and a reviewed plan binding the fresh paths and hashes. Do not edit historical plans in place or reuse a completed output directory.

This release verifies exported source identity and documentation. The 29 producer and nine independent helper checks are retained historical native receipts, not a combined validation score. The latter do not qualify the full checker end to end. The actual diagnostic and separate full saved-output readback completed on Mac Studio 512S2 before this source release; no caller, patient analysis, speech request or benchmark is run by preparing or checking this archive.
