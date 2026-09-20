# Completed-shard handoff audit

`audit_discovery_handoff.py` checks a snapshot of completed autosomal shards from
`call_genome.py` before a downstream screen consumes them. It reads existing files
and writes a new audit directory. It does not restart a caller, rewrite a VCF, merge
callers, select candidates or estimate a causal probability.

The current calling workflow retains native multiallelic GATK VCF and bcftools BCF
files, but its published shard outputs split multiallelic records. Downstream work
must preserve the original full ALT list, genotype, allele depths and likelihoods.
The original source paths and newly computed identities are included in this audit.

## Preconditions and reproduction

Use the pinned native environment from `workflow/environments/20260919/` on a Mac
Studio. Supply the exact reference, its existing index, a known calling-definition
hash, and the existing quality-screen configuration. The implementation currently
requires the `call_genome.py` normalization command contract and 22 explicitly named
autosomes. It refuses a shard mixing autosomes and other contigs.

```sh
MVA_ROOT=/path/to/project
"$MVA_ROOT/tools/bio/bin/python" workflow/scripts/audit_discovery_handoff.py \
  --run "$MVA_ROOT/calls/a-recorded-discovery-run" \
  --reference "$MVA_ROOT/ref/reference.fa" \
  --bcftools "$MVA_ROOT/tools/bio/bin/bcftools" \
  --quality-config workflow/configs/rare-biallelic-baseline.json \
  --expected-definition-sha256 RECORDED_SHA256 \
  --output /path/to/a-new-audit-directory
```

The default names are `1` through `22`. For a compatible prefixed reference, pass
`--autosomes chr1,chr2,...,chr22` with all 22 names explicitly expanded. Add
`--require-complete` when a consumer requires every autosome; without it, pending
autosomes are recorded and the result is explicitly a partial snapshot. A terminal
`status: complete` means this audit finished over its stated snapshot, not that
whole-reference discovery or interpretation finished.

The quality configuration must name the same sample as the calling definition.
Its existing rules are counted on all original raw records **before annotation**:
explicit PASS, fully called diploid GT, at least one called ALT, and its recorded
minimum DP/GQ. Missing fields remain missing. FILTER `.` does not become PASS;
missing GQ is not inferred from PL. This diagnoses input compatibility and does not
approve those rules for new callers or count coding/causal candidates.

## What is checked

1. Check the provided definition hash, reference hash/dictionary, completed-shard
   hashes, command-history identities and normalized indexes. Traverse each declared
   whole contig with an explicit indexed query; do not use no-region `fetch()` as
   an index test. Snapshot receipt
   membership before processing; later completions belong to a new audit.
2. Fingerprint each retained native source. The older completed-shard receipt did
   not hash those native files, so their newly recorded identities are distinguished
   from the previously pinned normalized outputs.
3. Count native GT/ploidy, multiallelic and retained-but-uncalled ALT fields, quality
   missingness and the existing screen's nonexclusive rejection reasons. This
   retains nonconcrete ALT accounting without treating those symbols as sequence.
4. Replay the original reference-checked split normalization to a stream, using the
   installed bcftools. Require the count and SHA-256 of **every ordered rendered
   record** to match the frozen shard, including genotype/sample fields. Header
   command lines and compressed bytes are not expected to match a streaming replay.
5. Check that recorded input sizes and modification times remain unchanged at the
   end. Save complete/failure state explicitly; never overwrite an existing output.

This consistency check uses the same installed normalizer and decoder. It does not
independently validate every allele/haplotype, establish biological accuracy, or
resolve mapping, callability, phase, mosaicism or the differential. In particular,
the separate X/Y and mitochondrial analyses are not replaced by the all-contig
diploid discovery outputs. Absence of a record or unfinished chromosome is never
reference evidence.

## Native fictional controls

```sh
MVA_PROJECT_ROOT=/path/to/project \
  /path/to/project/tools/bio/bin/python -m unittest discover \
  -s workflow/tests -p test_discovery_handoff.py -v
```

The controls build a tiny reference and frozen calling receipts, preserve a raw
two-alternate genotype, check missing quality versus explicit PASS, refuse altered
raw genotypes, a readable wrong index and invalid/incomplete scope, and cover an empty valid shard. They
qualify the handoff logic; they do not measure variant-calling sensitivity.
