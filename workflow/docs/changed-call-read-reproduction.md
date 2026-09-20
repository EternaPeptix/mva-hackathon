# Review read evidence for changed original calls

`review_changed_call_reads.py` compares selected original caller hypotheses against
one existing alignment. It extracts unfiltered indexed regions, preserves caller
GT/ALT/AD/PL context, and reuses the qualified SNP-fragment and exact-window helpers.
It does not run a genotyper, coerce ploidy, infer reference from absence, or assign
causal probability. Case selections, raw reads and outputs are supplied separately
from the private source-only snapshot.

The executed runtime was the existing Apple Silicon Mac Studio environment:
Python 3.12.14, pysam 0.24.1, samtools 1.24 and the pinned GRCh38 reference. The
fictional [native-pysam controls](../source-snapshots/changed-call-read-controls-20260920/controls.json)
exercise competing alternates at one position, original haploid multi-allele
context, overlapping mates, read-group separation, filters and explicit rejection
of unknown mapping quality. They are method controls, not biological validation.

Run from an owned project root with software already installed. Supply UTF-8 JSON
with `sample`, `expected_read_groups`, `inputs`, `helper_sha256`, `caller_inputs`,
`targets`, `slice_padding_bp`, `caller_context_padding_bp`,
`mapping_quality_thresholds` and `window_flanks_bp`.

- Each `inputs` entry has a `path` and expected `sha256`. Required input keys are
  `alignment`, `reference`, `samtools` and each name listed in `caller_inputs`.
  Include alignment/reference/caller indexes and relevant source receipts as
  additional pinned inputs. Relative paths resolve under `--root`.
- `helper_sha256` maps the qualified imported helper filenames to their hashes.
  Use the source snapshot identities for the helpers imported by the executed
  version. The completed output records both helper and main-script identities.
- Each target has `chrom`, one-based `pos`, concrete distinct `ref`/`alt`, and
  `key` exactly `chrom:pos:ref:alt`. Multiple alternate hypotheses at one coordinate
  are permitted; duplicate keys and inconsistent reference representations are not.
- The executed sensitivity settings were MQ 0, 20, 40 and 60, with BQ20; indel
  flanks were 1/2, 10/10 and 25/25 bp. These overlapping views must not be summed.
  The case used 1,000 bp extraction padding and 50 bp caller-context padding.

```sh
MVA_ROOT=/absolute/path/to/project
MVA_REPO="$MVA_ROOT/repo"
MVA_PY="$MVA_ROOT/tools/bio/bin/python"
MVA_SELECTION=/absolute/path/to/authorized-selection.json
MVA_OUT="$MVA_ROOT/results/new-changed-call-review"

"$MVA_PY" "$MVA_REPO/workflow/tests/test_changed_call_reads.py" -v
"$MVA_PY" "$MVA_REPO/workflow/scripts/review_changed_call_reads.py" \
  --root "$MVA_ROOT" --config "$MVA_SELECTION" --output "$MVA_OUT"
```

Use a previously absent output directory outside the source tree, and stop after
any failed command. The runner creates `status.json` before verifying inputs, so
check its live PID/process when an observation is interrupted. Do not relaunch
solely because a connection failed. Completion requires `status: complete`, intact
input identities and verified `output_sha256` values. Retain failed stages intact.

The SNP path recomputes BAQ, excludes orphan pairs and collapses eligible reads by
read group/query name; discordant mates remain conflicts. Exact-window counts use
raw BQ across aligned window boundaries and retain otherwise eligible improper
pairs. Consequently their cohorts differ. Unknown MQ255 is rejected before using
the SNP helper. Raw profiles are read counts, not independent fragment counts or
caller depth. Nearby changes, clipping, incomplete windows and repeat-equivalent
alignment can prevent a literal single-edit match without proving a reference
allele. High reported mapping quality also does not independently prove locus
assignment in homologous regions.

The retained annotation-overlap script in
[the source snapshot](../source-snapshots/changed-call-read-controls-20260920/annotate-local-context.py)
documents the original root and Ensembl release-116 GFF3 pin. It must be adapted
with recorded path changes for another dataset/root. Its one-based interval
overlaps and ±100 kb gene-span context do not predict regulatory targets or
pathogenicity. The fixed X/Y comparison that motivated the review is retained in
[its separate source snapshot](../source-snapshots/fixed-ploidy-review-20260920/review-fixed-case.py).
