# Native coding review: explicit inputs and preserved genotypes

The [preparation workflow](../scripts/prepare_native_coding_review.py) turns a
verified snapshot of completed autosomal calling shards into a caller-specific
protein-altering allele inventory and descriptive recessive models. It retains
unfiltered records, missing quality fields and alleles absent from the native
genotype. It does not estimate EPCR, rank diagnoses or establish pathogenicity.

This document describes the method qualified on **20 September 2026**. It contains
no case results. A successful run means that the stated snapshot was processed;
it does not mean that genome-wide calling or interpretation is complete.

## Inputs and command

Use the native Mac Studio environment and a completed
[discovery handoff audit](discovery-handoff-reproduction.md). The handoff supplies
the exact native sources, their hashes, record counts, rendered-record payload
hashes, sample and selected/pending chromosomes. The preparation workflow accepts
autosomes named `1` through `22`; it requires their selected and pending lists to
be disjoint and exhaustive. It never discovers additional completed shards while
running. Produce a new handoff and new output directory to extend the snapshot.

The reference and existing index must match the handoff. Supply the exact GFF3
annotation hash and a reviewed configuration whose sample matches the handoff.
The configuration specifies consequence groups, frequency views and descriptive
quality handling. The tested consequence source is the Ensembl release 116 GRCh38
primary-chromosome GFF3. The output directory must not already exist.

From the repository root, adapt these paths and recorded hashes:

```sh
MVA_ROOT=/path/to/project
"$MVA_ROOT/tools/bio/bin/python" workflow/scripts/prepare_native_coding_review.py \
  --root "$MVA_ROOT" \
  --handoff /path/to/completed-handoff/summary.json \
  --expected-handoff-sha256 RECORDED_HANDOFF_SHA256 \
  --caller gatk \
  --reference /path/to/reference.fa \
  --gff /path/to/pinned-annotation.gff3.gz \
  --expected-gff-sha256 RECORDED_GFF_SHA256 \
  --config /path/to/reviewed-native-coding-config.json \
  --output /path/to/new-gatk-coding-review
```

Use `--caller bcftools` with a different new output directory for the second
caller. Add `--allow-partial-snapshot` explicitly when autosomes remain pending.
Without that flag, incomplete autosomal scope is refused. This flag acknowledges
scope; it does not establish callability on the included chromosomes.

The default frequency query uses these retained project-relative dependencies:

- `tools/bio/lib/jvm/bin/java`
- `lirical/lirical-cli-2.4.1/lirical-cli-2.4.1.jar`
- `exomiser/2512_hg38/2512_hg38_variants.mv.db`
- [AlleleFrequencyProbe.java](../scripts/AlleleFrequencyProbe.java)

For exact precomputed responses, supply both `--frequency-responses /path/to/responses.jsonl`
and `--expected-frequency-sha256 RECORDED_RESPONSE_SHA256`. There must be exactly
one response for each queried normalized allele key, with no extra keys. A file
hash establishes identity, not the provenance or correctness of supplied counts.

## Native provenance and annotation

Only a derived copy receives the one-based `NC_SOURCE_ORD` record ordinal. Before
annotation, its complete ordered rendered-record payload, after removing that
added field, must match the original source and the completed handoff census.
Original headers are saved separately. Native input files are not indexed,
rewritten or assigned new quality values by this workflow.

Reference-checked `bcftools norm -m -any --old-rec-tag NC_ORIG` creates normalized
allele rows. `bcftools csq --local-csq --samples -` annotates their local sequence
consequences. It does not infer compound coding haplotypes. A before/after payload
comparison requires every nonconsequence rendered record field to remain the
same across annotation; only BCSQ annotations may differ.

The workflow then reads the original native source sequentially again and recovers
each selected coding row by source ordinal. Its stable native record identity is
caller, source SHA-256 and ordinal. This distinguishes otherwise identical duplicate
records. The [record helper](../scripts/native_coding_records.py) preserves the full
original ALT list, all sample FORMAT values including GT/AD/PL, phase flag,
PGT/PID/PS when present, FILTER, QUAL and the original rendered record. Split GT,
AD or PL fields never replace that native evidence.

The original-record tag must match the source chromosome, position, REF, full ALT
list and original ALT index. A comma-containing String tag may arrive from pysam
as a tuple; its string elements are rejoined before strict validation. An absent
tag is accepted only for an unchanged biallelic source identity. For every selected
concrete coding allele, applying its original and normalized sequence edits to
the reference must produce the same sequence across their combined span. This
checks the selected allele independently of its origin tag, including left shifts.

Nonconcrete rows, including `*`, remain in a representation-review output with
their source ordinal. A spanning-deletion `*` is not an independently interpretable
local sequence deletion or reference allele. Concrete rows lacking BCSQ or the
configured protein-altering effects are accounted for separately; this is not an
inventory of every functional or disease-relevant variant class.

## Quality, genotype and frequency distinctions

Every selected concrete protein-altering row remains in `coding-alleles.jsonl.gz`,
including uncalled ALT, incomplete or non-diploid GT, named-filter records and
missing or low quality values. The native call has independent review axes:

| Field or property | Retained distinctions |
|---|---|
| FILTER | explicit PASS; unfiltered `.`; original named filters |
| DP and GQ | absent FORMAT field; present but missing; explicit zero; below or at/above the historical diagnostic threshold |
| Genotype | fully called diploid; partial/missing; haploid; other ploidy |
| Reviewed ALT | observed in GT; absent from a fully called GT; not observed in an incomplete/missing GT |
| Representation | nonconcrete alleles in GT; colliding normalized keys; original record identities |

DP 10 and GQ 20 are **historical diagnostic bins only**. No PASS label or GQ is
imputed, and PL is not converted to substitute GQ. FILTER, DP, GQ and QUAL do not
determine eligibility for the descriptive models. These choices do not validate a
new caller-specific quality threshold or establish precision, recall or sensitivity.

Public frequency is the maximum valid AC/AN fraction over available retained
gnomAD exome/genome population counts. Source counts and denominators remain
attached. An absent database allele, absent usable counts or AN=0 remains
frequency-unresolved with a null maximum. AC=0 with positive AN is an observed
zero within that dataset, not proof that the allele is absent from every population.
Invalid counts or missing query responses are refused rather than converted to zero.

The configured AF views are 0.001, 0.0001 and 0.01. Observed alleles above a view's
threshold are omitted from that model view but retained in the allele inventory.
Unresolved alleles remain in a separate frequency-unresolved model group. These
views do not establish population callability, pathogenicity or validated clinical
sensitivity. Querying an identical normalized key once does not merge its native
record identities.

## Descriptive model contract

Ordinary models require a called ALT, a fully called diploid native GT, concrete
sequence alleles for the genotype and no normalized-key collision. Thus `./1`,
haploid `1` and a concrete ALT paired with `*` remain review evidence rather than
ordinary heterozygous models. All colliding source records remain in the inventory;
the workflow does not choose one by quality or silently discard their evidence.

Models stay within one caller and chromosome and require a shared affected
transcript in the same gene. A homozygous ALT yields one model whose consequence
is counted twice. Distinct heterozygous alleles require their own source evidence.
Two concrete ALTs from the same native diploid `1/2` record are labelled
`opposite_alleles_in_source_GT`. This is an assertion of that recorded genotype,
not parental phase or independent validation of the call.

Source phase blocks are namespaced by caller, sample and source-file hash.
Identical PS or PID values in different callers or source streams cannot establish
phase. PGT must agree with the GT allele multiset. Comparable PS and PID signals
that conflict leave phase unknown with the conflict recorded. Compatible explicit
source-block signals may describe cis or opposite alleles; they do not establish
parental origin. Unphased calls stay unknown.

Explicit cis pairs and overlapping representations are retained in separate pair
review files and omitted from ordinary models. The overlap exception is the two
concrete ALT alleles of the same fully called diploid native record. Consequence
tiers are descriptive: two predicted-disruptive alleles, one predicted-disruptive
plus one other protein-altering allele, or two other protein-altering alleles.
Each model retains its supporting transcripts, full native calls, frequency state
and quality review axes. Display ordering does not create a clinical rank.

## Outputs and interpretation limits

Each new run writes status, summary, pinned input identities, commands and output
hashes. Per-shard files retain original headers, derived/normalized/annotated BCFs,
annotation logs and nonconcrete representation review. Top-level outputs include
the full coding allele inventory, normalized-key collisions, frequency queries and
responses, descriptive models for each AF view, and retained cis/overlap pair
reviews. Input sizes and modification times are checked again before completion.

Selected and pending autosomes remain explicit in status and summary. A record
absent from a variant-only file is not callable-reference evidence; an unfinished
chromosome is not a negative result. Separate callers must be reconciled before
claiming agreement, discordance or a unique biological event count. Partial model
lists must be extended and reviewed after the remaining calling scope finishes.

This method does not validate loss of function, transcript relevance, disease-gene
relationships, parental phase, mosaicism, mapping origin or a diagnosis. It does
not replace phenotype review, callability, X/Y, mitochondrial, structural,
copy-number, repeat or regulatory analyses. Neither the inventory nor its tiers
provide EPCR values or clinical validation.

## Native fictional qualification

The immutable `qualification-v1` control run completed on **Mac Studio 512S2 at
2026-09-20 13:50:48 UTC**: **18 tests passed, none skipped**, exit code 0. The
recorded runtime was Python 3.12.14, pysam 0.24.1 and macOS 26.4.1 arm64. The
receipt SHA-256 is
`c97db2bcf6ad0810a0a1717f77170a569724ee10e15bcb08d9666d302c07bcdd`.

The suite comprises 13 [record/model controls](../tests/test_native_coding_records.py)
and five [native preparation controls](../tests/test_prepare_native_coding_review.py).
It exercises native VCF/BCF decoding, comma-containing origin tags, actual
normalization/local consequence annotation, both ALTs and full PL in `1/2`,
duplicate identities, missing/zero quality, partial/haploid/star handling,
caller/source phase separation, conflicting phase evidence, frequency unknown
versus observed zero, cis/overlap review, explicit partial scope, an empty valid
snapshot, input hash refusal and repeat-shifted allele sequence equivalence.

Qualified source identities:

| File | SHA-256 |
|---|---|
| `prepare_native_coding_review.py` | `1fbf7335578b865df013f6299ebaa422d90387db99cd8e265c29bd787ccb6da4` |
| `native_coding_records.py` | `c395dd130a2577e7b25ba0ac62d600d7cb64c9240f3eb218536eff22c5c44cfe` |
| `test_prepare_native_coding_review.py` | `b1dfc5c5352781574b0039b6d3796200fcfba95966773826e3ee45ada6a959b7` |
| `test_native_coding_records.py` | `a96bff7c20f31fa4a120e32f6a542f7eb86874646e6c3c065a0d6cb97854ff48` |

Reproduce the combined controls from the repository root on the Studio runtime:

```sh
MVA_PROJECT_ROOT=/path/to/project \
  /path/to/project/tools/bio/bin/python -m unittest discover \
  -s workflow/tests -p 'test_*native_coding*.py' -v
```

Unindexed temporary native BCFs can emit an index-lookup diagnostic during these
sequential-read controls. Qualification completed with that diagnostic retained.
The fictional controls qualify representation and workflow behavior; they do not
measure biological sensitivity, frequency-database correctness, clinical accuracy
or treatment effects.
