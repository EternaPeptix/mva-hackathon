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
