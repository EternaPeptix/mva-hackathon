# Full-query clipped-read review

This workflow compares complete available seed-read sequences with bounded local
sequence hypotheses, then separately remaps canonical seed and mate ends against
the whole reference. It retains clipping, qualities, ambiguous fits and reported
alternative placements. It does not call a variant, infer a diploid haplotype
pair, establish genomic origin or parental phase, or estimate clinical benefit.

The implementation is [review_clipped_reads.py](../scripts/review_clipped_reads.py),
with exact sequence fitting in [clipped_read_helpers.py](../scripts/clipped_read_helpers.py).
It reuses raw-source collection from
[review_joint_haplotypes.py](../scripts/review_joint_haplotypes.py), per-caller
sequence enumeration from [joint_haplotype_helpers.py](../scripts/joint_haplotype_helpers.py),
and original caller fields from [native_coding_records.py](../scripts/native_coding_records.py).
The retained private [configuration](../configs/nabp2-clipped-read-review.json)
provides input identities and the intended selection. This document contains
method and fictional-control information only, with no case read identities,
read sequences or results.

## Input contract and output preservation

The runner accepts `--root`, `--config` and a new `--output` directory outside
the source tree. Relative input paths resolve against the project root. The
configuration pins:

- The original alignment and explicit alignment index, an existing regional
  seed BAM and its explicit index, and the reference and associated FAI.
- Both retained raw caller sources, complete source-record counts and ordered
  rendered-record SHA-256 values, the sample and expected read-group count.
- The BWA binary, its index receipt and every whole-reference index file.
- The runner and all four imported workflow helpers listed above.
- The seed interval, target cores, flanks, required source records and work caps.

Input hashes are read in 8 MiB blocks. Configured sizes are checked, and file
size/modification time must remain unchanged while hashing and through the run.
The reference FAI path must correspond to the pinned reference, and BWA index
identities must bind to that reference. No raw call, original alignment, index
or earlier output is rewritten. An existing output directory is refused.

The new output retains original ends and mate retrieval outcomes, the complete
regional native context and source census, sequence panels and cached fits,
FASTQ exports, BWA SAM records and logs, commands, status and summary. The
summary fingerprints outputs; completed status links to its summary hash.
Exceptions leave failed status and retained partial evidence. A corrected run
uses a new directory and its own source/configuration identities.

## Seed records and canonical ends

The exact configured seed interval is reread from both the original alignment
and pinned seed BAM. Their ordered decoded SAM records and SQ/RG headers must
match. Raw seed records remain preserved, including secondary, supplementary
and excluded-flag context. This proves extraction identity for the selection;
it does not validate the original mapping or enumerate every potentially
relevant read elsewhere in the dataset.

Primary ends are identified by `(RG,QNAME,end_number)`. Missing identities,
paired records without an end number, both end-number bits, or repeated primary
identities are refused. Mapped mate retrieval checks a unique primary end and
reciprocal coordinates. A mate that is already a seed must have the same full
SAM payload and is retained once with both roles. Unpaired/declared-unmapped,
failed/ambiguous lookup and nonreciprocal mates remain explicit; unsuccessful
coordinate lookup does not prove that a mate is absent from the original data.

The original alignment-oriented sequence and raw qualities, forward-oriented
sequence and corresponding forward qualities, CIGAR, terminal soft-clip slices,
flags, MQ, strand, read group and mate geometry are retained. Reverse-strand
forward conversion reverses the qualities as well as reverse-complementing the
sequence. No clipping, consensus correction, quality replacement, BAQ or
recalibration is applied to these retained inputs.

“Full query” means the complete **available** query. Soft-clipped bases are part
of that query; hard-clipped bases are absent. The CIGAR and available length
must be inspected before calling a record a complete sequenced read. The runner
does not recover missing hard-clipped sequence.

Only primary **seed** ends receive local fits. They exclude flags `0xF04`, MQ255
and missing/nonconcrete ACGT sequence. Low numeric MQ and absent/low BQ are kept
as recorded observations; qualities do not enter the distance score. Recovered
mate-only ends receive original-context review and whole-reference remapping,
not local panel scores.

## Local panels and fixed bounds

For target core `[a,b)` and flank `f`, the reference panel is `[a-f,b+f)`.
Coordinates are zero-based and half-open. An out-of-reference panel fails;
there is no silent coordinate clipping. The current configuration sets:

| Bound | Value |
| --- | ---: |
| Reference flanks per panel | 250 and 500 bases |
| Raw seed records | 512 |
| Canonical primary seed/mate ends | 1,024 |
| Attempted Cartesian hypothesis states across callers, per panel | 256 |
| Bases per materialized hypothesis | 100,000 |
| Nominal cells per unique full-query fit | 1,000,000 |
| Optimal CIGAR witnesses per fit | 64 |
| Nominal cells summed over uncached fits | 500,000,000 |
| BWA threads | 2 |

Both complete raw caller files are scanned sequentially. The expected record
counts and ordered rendered-record hashes must match. Maximum-window collection
uses reported-end/literal-REF overlap, including upstream long records that a
POS-only regional query could miss. The original caller, file identity, ordinal,
all ALTs, GT/phase and complete FORMAT, quality/filter fields and rendered record
remain available.

For each panel, complete concrete literal-REF-contained ALTs enter enumeration,
including uncalled ALTs. Crossing or nonconcrete alleles are retained as omitted
when selected by the literal-span rule; they are never clipped into a variant.
The wider context can contain an upstream symbolic END-spanning record absent
from an individual panel's literal-span omission list. That list alone is not a
structural-variant census.

Enumeration is separate within each caller, selecting no ALT or one ALT per
source record. It does not mix callers or take two ALTs from one multiallelic
record. Intersecting literal REF spans are rejected and retained as rejected
combinations; this conservative representation rule is not proof of biological
incompatibility. Rejected combinations still count against the state bound.
Overflows fail rather than drop an inconvenient subset.

Identical constructed sequences share one hypothesis with all source lineages
attached. A caller with no records contributes a reference sequence hypothesis,
not a reference genotype. Native zygosity and phase are retained as source
context; enumeration is not a probabilistic inference of the sample's two
haplotypes. The panel does not include every structural, assembly or genomic
origin hypothesis that might explain a clip.

## Exact full-query fitting

`fit_full_query(query,target,max_cells,max_paths)` implements query-global,
target-local unit edit distance: match costs zero, mismatch/insertion/deletion
cost one, and unused target flanks are free. Every query base is consumed.
Terminal query insertions have their ordinary cost; no query clipping is free.
An empty target substring is a legal costly insertion-only placement, but
empty input strings and non-uppercase-ACGT inputs are refused.

The runner fits original SAM alignment-oriented queries against the panel's
forward sequence. This restriction is explicit. Genome remapping separately
provides reported orientation alternatives; the local fit is not a genome-wide
or exhaustive two-orientation search.

Before allocating the DP grid, the helper searches for every overlapping exact
substring occurrence. A zero-cost alignment can contain only matches, so these
occurrences exhaust all optimal placements. This fast path changes neither the
score nor the ambiguity. Cell/path limits are enforced before it runs. For a
nonexact fit, the helper uses bounded DP and iterative traceback with saturated
path counting, avoiding recursive or unbounded exponential traceback.

The cell budget is `(query_length+1)*(target_length+1)`, including the boundary
row and column. Each distinct query/target string pair is charged once, even if
shared by several physical ends or panels. The runner reports nominal
`budgeted_fit_cells` separately from `evaluated_DP_cells`: exact-substring fits
charge their nominal grid but evaluate no DP grid. These counts do not by
themselves measure an end-to-end speed gain.

Each fit preserves the exact minimum cost and **every optimal target endpoint**.
Each returned CIGAR witness has a unique `(start0,end0,cigar)` identity, consumes
the complete query and can be replayed against its target segment using
`=`, `X`, `I` and `D`. The helper's general limits default to 1,000,000 cells and
256 witnesses, with hard ceilings of 10,000,000 cells and 4,096 witnesses; the
runner uses the stricter configured values above.

When more than `max_paths` alignments exist, the witness list is capped and
`optimal_path_count_capped` is `max_paths+1`, a lower bound.
`path_count_is_exact` and `traceback_complete` are false, `saturated` is true and
status is `capped`. Otherwise that count is exact and traceback is complete.
`unique_optimal_alignment` is true only when exactly one optimum exists.
The exact-substring path additionally supplies `optimal_path_count_exact` even
when its witnesses are capped; this is separate from the capped-count field.
For a saturated nonexact DP fit, that additional exact count is unavailable.

The full query-by-hypothesis cost matrix, reference cost and every tied best
hypothesis ID remain in the output. There is no automatic allele-support label,
acceptance threshold or likelihood. A small improvement against a proposed
insertion may still leave a large unexplained suffix. A long clip is not
discarded to make a small-variant hypothesis fit.

## Ties, coordinates and quality interpretation

Distinct hypotheses can tie because the query never reaches the bases at which
they differ. A unique best fit can also occur elsewhere within the local panel,
away from the candidate. Multiple callers can describe the same constructed
sequence. None of these cases establishes an observed multi-variant haplotype.

Fit `start0/end0` values are offsets in the **constructed target sequence**.
For an indel-containing hypothesis, adding the reference panel start is not a
valid general conversion to genomic coordinates. Reference coordinates and
allele-discriminating coverage require the hypothesis edit-coordinate mapping.
A capped witness set cannot establish the start/CIGAR/allele coverage of every
optimal path. One favorable traceback must not be selected as definitive.

Unit edit costs ignore base qualities. The saved raw-quality vectors and
summaries must be consulted before interpreting a fit driven by a low-quality
base. Missing qualities are not high qualities, and MQ255 is unavailable rather
than a high mapping score. The scores are not quality-weighted genotype scores,
probabilities, EPCR, or calibrated clinical evidence.

## Whole-reference remapping and fragment denominators

Each available canonical end with valid forward sequence and FASTQ-representable
qualities is exported once per applicable arm, with a reversible synthetic
query ID. Missing/unrepresentable qualities remain an explicit remap-ineligible
reason; no qualities are fabricated. Original flagged contextual ends can enter
the single-end arm, while the paired arm requires flag-eligible original pairs,
both original paired bits, neither mate declared unmapped and reciprocal
coordinates. Selection scopes, original-flag counts and exported flagged-end
counts are retained separately.

Both arms use the pinned whole-reference BWA `mem -a -Y` with their own exports,
commands, logs and SAM outputs. All reported primary, secondary, supplementary
and unmapped records remain available. Exactly one primary output is required
per exported end. Whenever a remapped record contains sequence/qualities, its
forward payload must match the original export exactly; omitted secondary
payloads remain explicitly linked to the canonical end instead of receiving an
invented payload-verification claim.

BWA's reported alternatives are nonexhaustive. The paired arm re-estimates
insert-size behavior from a selected small subset and does not reproduce the
original whole-library alignment. The single and paired arms reuse the same
ends and reference and are not independent validation. Whole-reference remapping
also does not insert every alternate hypothesis into the reference.

Compare full-query coverage, residual clips, original/remapped MQ, score and
mate geometry alongside alternative and split records. An unreported competitor
does not give an infinite score margin. A short supplementary segment is not
interchangeable with a full-query alignment score, and overlapping split query
segments must not be added into independent support.

The runner records unique ends and `(RG,QNAME)` fragment identities; it does not
assign allele-support counts to fragments. A later recount must count one
fragment across agreeing mates, panels, hypothesis lineages and remap arms.
Different full strings from two mates are normal and are not themselves a
conflict. Mate-only context is not a variant-spanning allele vote. Raw read
counts, unique ends, identifier-group fragments and independent biological
molecules are different denominators.

## Reproducibility and qualification status

The [pure helper controls](../tests/test_clipped_read_helpers.py) cover exhaustive
short-string substring/edit-script oracles, replayed CIGAR costs, exact matches,
query-end penalties, a fictional full-length two-base insertion, an unexplained
long suffix, repeated placements, caps and invalid inputs. The
[native integration controls](../tests/test_clipped_read_workflow.py) exercise
actual synthetic BAM/CRAM, VCF/BCF, BWA indexes and remapping, mate recovery,
query/quality preservation, alternate loci, ties, missingness and refusal paths.
They do not run a genotyper or read case inputs.

The retained [source-only native qualification receipt](../source-snapshots/clipped-read-controls-20260920/controls.json)
records **12 tests passed, none skipped**, on Mac Studio 512S2 at
`2026-09-20T15:21:44.485903+00:00`: seven pure-helper controls and five native
integration controls. It records Python 3.12.14, pysam 0.24.1 and macOS 26.4.1
on arm64, together with exact hashes of the five source files, private
configuration and two test files. The receipt's SHA-256 is
`8f0a7684cf32e657cedd85e4a67b3f479ec20bb054bd82197f309b700cdfa4d4`.
This qualifies that frozen method on fictional inputs; it is not a case-result,
biological sensitivity or clinical validation claim.

From the source repository, use the installed Mac Studio environment. Set the
root to the private runtime that supplies Python/pysam, samtools and BWA. This
strict control invocation rejects skipped or empty suites as well as failures:

```sh
export MVA_PROJECT_ROOT=/path/to/private/project
"$MVA_PROJECT_ROOT/tools/bio/bin/python" - <<'PY'
import unittest
suite = unittest.defaultTestLoader.discover(
    'workflow/tests', pattern='test_clipped_read*.py')
result = unittest.TextTestRunner(verbosity=2).run(suite)
if not result.wasSuccessful() or result.skipped or result.testsRun == 0:
    raise SystemExit(1)
PY
```

For a configured review after qualification, use a new output directory:

```sh
"$MVA_PROJECT_ROOT/tools/bio/bin/python" workflow/scripts/review_clipped_reads.py \
  --root "$MVA_PROJECT_ROOT" \
  --config workflow/configs/nabp2-clipped-read-review.json \
  --output "$MVA_PROJECT_ROOT/results/NEW-CLIPPED-READ-REVIEW"
```

Changing data or parameters requires a separately reviewed configuration and
matching source/input hashes. Record the exact source/config/test/runtime
identities and require the entire native suite to pass with no skips before
qualifying a revision. Preserve failed qualification directories.

A successful fictional qualification is not a case-result check or biological
validation. Each actual run requires its own completed status/summary linkage,
input/output hash verification, original-record/end/export census, cost and
CIGAR reconciliation, cap interpretation and fragment recount. No independent
specimen, parental phase, functional validation, causal ranking or annotation-
engine validation is provided by this workflow.
