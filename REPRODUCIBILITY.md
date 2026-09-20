# Reproduction scope and order

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
base; Truvari 5.5.0 uses an isolated environment. The `*-local.wheels.txt` files require
the exact tested wheel hashes. Forty wheels matched PyPI publisher artifacts; five
were local builds whose source-archive identities matched PyPI. See the
[wheel-source record](results/runtime-wheel-sources-20260919.json).

The exact 45-wheel recovery bundle remains on owned storage and is **not distributed
in this initial source snapshot**. Its SHA-256 is
`2cf09339460603fee2149461d21f6f0ecc9819360c385afcb7c27db2165791ae`.
A new source build can change wheel bytes and requires renewed validation; it is not
the hash-locked installation tested here. Until the bundle or newly validated builds
are available to a reproducer, exact Python-environment recreation is incomplete.
This limitation must remain visible in the final reproducibility review.

The scripts `export_runtime.py`, `verify_runtime_rebuild.py`,
`prepare_runtime_wheels.py`, `record_runtime_wheel_sources.py` and
`verify_python_rebuild.py` capture and test the installation. The completed
[Conda](results/runtime-conda-rebuild-20260919.json) and
[Python](results/runtime-python-rebuild-20260919.json) checks reconstructed seven
environments and passed 22 selected existing controls. Exomiser/LIRICAL dependencies
were inventoried, not freshly installed in that exercise.

## 2. Reproduce the public small-variant comparison

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

Obtain the challenge data through its access process and retain it outside the Git
checkout. The internal analysis used all eight source FASTQs from four lanes, with
publisher-checksum verification. Do not create an inferred clinical dataset or treat
the public HG002 reads as patient data.

- `01_align.sh` performs reference-aware BWA alignment and duplicate marking,
  preserves separate lane artifacts and validates the completed CRAM.
- `call_genome.py` runs native diploid GATK and bcftools discovery by whole-reference
  contig shards. `resume_reads.py` and `resume_calling.py` are dependency-aware
  controllers for the original setup; inspect their configuration before reuse.
- `prioritize.py` and the Exomiser configurations perform the phenotype-assisted
  analyses. `run_lirical.py` and related summary/resource helpers provide the separate
  disease comparison. `rare_biallelic_baseline.py` and transcript helpers implement
  the phenotype-independent review. Exact resource versions and paths must be supplied.
- Read, mate and WhatsHap helpers retain source genotype/phase context. Full original
  allele lists, GT, AD and PL are needed alongside split allele rows; a matching
  component must not hide an incompatible second allele or imply parental phase.
- The CNVpytor, structural and ROH workflows require final alignment/QC, appropriate
  masks, ploidy/callability review and their own interpretation. An ordinary diploid
  callset does not establish mosaic, mitochondrial, sex-chromosome or structural coverage.

At this source snapshot the four-lane patient alignment was still completing and
full-genome callsets were not available. The patient-derived ranking/read summaries,
draft submissions, final reports and video are not included here. They remain separate
work products awaiting their stated review; this source repository does not complete
those submission requirements.

## 4. Check components and preserve evidence

The `workflow/tests` directory contains unit and native integration checks. Some
require external runtimes or previously completed public-control artifacts. A skipped
test is not a passing native test. Run the tests relevant to the installed component,
inspect skips and native outputs, and keep the executed source/configuration hashes.

The source snapshot copy is byte-checked. Its only alignment-script change is a
documented comment replacement removing an unnecessary acquisition-specific header
and barcode example; executable lines are identical. Python compilation and shell
syntax checks validate the copy's syntax, not end-to-end scientific correctness.
Controlled sequencing data, indexes and genotype-scale intermediates are not inputs
to this Git repository and must not be committed to make tests convenient.
