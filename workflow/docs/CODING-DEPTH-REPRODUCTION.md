# Reproducing fixed-gene coding-depth measurements

This source method measures aligned-read depth across the union of annotated CDS
intervals for a declared gene list. It does not call variants, infer parental phase,
assign a locus of origin or establish clinical callability or a negative screen.
The exact study configuration is retained under the existing private source-export
authorisation; alignment data and case results are not included.

## Source and resources

- [Target and denominator construction](../scripts/coding_depth_targets.py)
- [Native depth workflow](../scripts/measure_coding_depth.py)
- [Exact configuration](../configs/coding-depth-v2.json)
- [Fictional controls](../tests/test_coding_depth.py)
- [Native qualification record](../source-snapshots/coding-depth-controls-20260920/controls.json),
  [original receipt](../source-snapshots/coding-depth-controls-20260920/native-controls-original.json)
  and [complete test output](../source-snapshots/coding-depth-controls-20260920/native-test-output.txt)

The workflow reuses the unchanged depth collector and alignment/header checks
already included in the repository. Keep all four scripts named in the
configuration's `source_sha256` map together in `workflow/scripts`. The configuration
pins the alignment/index, reference/index, GFF3 annotation and samtools by SHA256.
Access to the controlled alignment must be established separately; the source
repository does not grant or supply it.

Absolute resource paths identify the original owned-storage layout. For another
machine, create a new configuration with reviewed local paths and the same verified
input bytes. Preserve the original configuration. Recompute the configuration's
SHA256 and pass it explicitly; do not weaken input or source checks to relocate
files. `gene_list_source` is retained provenance for the prior review list, not an
input used to infer the list from case depth.

## Selection and quality policy

Gene names must resolve uniquely to protein-coding genes on primary nuclear
chromosomes, with valid transcript/CDS parentage, strand and bounds. CDS intervals
from all annotated transcripts are unioned per gene. Disjoint segments preserve
every overlapping gene membership, count each panel base once and retain reference
ACGT versus non-ACGT denominators. No promoter, UTR, splice-flank, repeat or
whole-genome completeness is implied.

Indexed selection preserves the original alignment records and read groups. The
three declared views are MQ0/BQ1, MQ20/BQ20 and MQ60/BQ20. All use explicit zeros,
raw qualities without BAQ, samtools overlapping-mate suppression and flag exclusion
`0xF04`. Deletions and reference skips do not contribute base depth. This policy
does not independently establish unique molecules or genomic origin. The first
versus second view changes both quality axes; the second versus third holds BQ
fixed. Views share observations and must not be summed as independent evidence.

Samtools can count unavailable base qualities at a positive BQ threshold. The
workflow therefore refuses flag-eligible records with MQ255, absent qualities or
any quality-255 sentinel before producing depth reports. Explicit quality zero is
valid input data and is excluded by the declared positive BQ gates. Every expected
position must appear exactly once in order. Means with an empty ACGT denominator
remain null; bins at 1, 10, 20 and 30 are descriptive, not validated clinical gates.

## Reproduce qualification and run

The retained native qualification ran all ten tests on Mac Studio 512S2, using
the project's pinned Python/pysam and samtools environment. It includes literal
expected BAM/CRAM counts, mate overlap, quality and flag exclusions, CIGAR deletion
and skip behavior, uncovered targets, ambiguous reference bases, overlapping genes,
duplicate isoforms, input/hash refusals and preservation of existing outputs.
The receipt in this source update is a labelled derivation of the original receipt
and test log; creating it did not rerun the tests or establish new performance.

The native test fixture intentionally inventories every Python file in its sibling
scripts directory. Reproduce the qualified five-file layout in a new temporary
directory, rather than adding unrelated repository scripts to that inventory.
From the repository root, with the reviewed native dependencies installed:

```sh
MVA_DEPTH_CONTROLS=$(mktemp -d)
mkdir -p "$MVA_DEPTH_CONTROLS/workflow/scripts" "$MVA_DEPTH_CONTROLS/workflow/tests"
for MVA_DEPTH_SCRIPT in coding_depth_targets.py measure_coding_depth.py \
  survey_chromosome_depth.py prepare_gatk_compatible_cram.py; do
  cp "workflow/scripts/$MVA_DEPTH_SCRIPT" "$MVA_DEPTH_CONTROLS/workflow/scripts/"
done
cp workflow/tests/test_coding_depth.py "$MVA_DEPTH_CONTROLS/workflow/tests/"
MVA_SAMTOOLS=/absolute/path/to/samtools python -B -m unittest discover \
  -s "$MVA_DEPTH_CONTROLS/workflow/tests" -p test_coding_depth.py -v
```

Require all ten tests to run without skips, errors or failures. The earlier
qualification before the unavailable-quality guard is not the accepted method.
Run the analysis only with the reviewed configuration digest and a new output
directory:

```sh
python -B workflow/scripts/measure_coding_depth.py \
  --config workflow/configs/coding-depth-v2.json \
  --expected-config-sha256 REVIEWED_CONFIGURATION_SHA256 \
  --output NEW_OUTPUT_DIRECTORY
```

The workflow preserves source/input fingerprints, selection and denominators,
original selected records, commands, every reported position, window/gene/panel
histograms and terminal status. Check completion, source identity, full position
accounting and independent reconstruction before interpreting a result. A repeat
using the same samtools engine is technical reproduction, not a separate biological
assay. No actual coverage counts, reads or individual clinical conclusions are
included in this source update.
