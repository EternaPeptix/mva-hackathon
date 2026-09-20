# Local joint-sequence hypothesis review

This workflow compares constructed local sequences with original mapped read
fragments. It can reveal that different caller representations describe the same
sequence, or that nearby changes need to be interpreted together. It does not run
a genotyper, infer a diploid haplotype pair, assign parental phase, or estimate a
causal probability.

The implementation is [review_joint_haplotypes.py](../scripts/review_joint_haplotypes.py),
with sequence and fragment rules in
[joint_haplotype_helpers.py](../scripts/joint_haplotype_helpers.py) and original
source fields retained through
[native_coding_records.py](../scripts/native_coding_records.py).
This document describes the method and fictional qualification controls only;
case selection, reads, sequences and results require their own retained run
configuration and evidence directory.

## Inputs and identity checks

The command accepts `--root`, `--config` and a new `--output` directory outside
the source tree. The configuration supplies:

- Expected SHA-256 identities for the original alignment, reference, annotation,
  native caller files, tools and accompanying indexes. Configured sizes are also
  checked. Paths resolve against the project root.
- SHA-256 identities for the runner and its two source helpers.
- The exact sample, expected read-group count, caller-source labels, and each
  whole raw source's expected record count and ordered rendered-record SHA-256.
- Explicit chromosome and zero-based, half-open core coordinates for each
  target. Optional required source-record identities check membership in the
  collected regional context; they do not define the entire hypothesis list.
- Flank lengths, mapping/base-quality views and enumeration limits.

Files are read-only inputs. SHA-256 reads stream in 8 MiB chunks; size and
modification time must remain stable during hashing and until review ends. No
native call file is indexed, normalized in place, filtered or rewritten. The
configuration itself is fingerprinted in the output. A caller label is
provenance, not independent biological replication.

## Window and hypothesis construction

For core `[a,b)` and flank `f`, the comparison window is `[a-f,b+f)`. Both callers
use that same reference-coordinate window. Each requested flank remains a
separate view. The qualified configuration uses flanks 3, 10 and 25 bases, MQ
thresholds 0, 20 and 60, and raw BQ thresholds 20 and 30. These are sensitivity
views, not clinically calibrated acceptance thresholds.

The runner scans each entire retained raw VCF/BCF sequentially and verifies its
complete ordered record payload. Regional collection uses the maximum requested
flank and an overlap interval ending at the larger of reported record end and
literal REF end. This collects long upstream records that a POS-only query could
miss. Original caller, file hash, record ordinal, all ALTs, complete FORMAT,
GT/phase fields, quality/filter fields and the rendered source record are retained.

For each individual window, all overlapping concrete native ALTs whose complete
literal REF fits in the window enter enumeration. They enter even when absent
from the source GT or when quality information is missing. The `called` flag
records ALT participation in native GT; it does not make a constructed sequence
an observed genotype. Literal window-crossing or nonconcrete alleles are not
clipped or converted to reference.

Enumeration occurs separately within each caller. For every native record the
choices are no ALT or one of its ALTs. Combinations never mix caller labels and
never select two ALTs from the same source record. Intersecting literal REF spans
are conservatively rejected, with the rejected edit combinations retained.
This rule can reject representationally compatible changes sharing an anchor;
it is not proof that the underlying biological changes cannot coexist.

The qualified limits are 256 attempted Cartesian states across callers per
window and 100,000 bases per materialized sequence. Rejected overlap states count
toward the state cap. Exceeding a cap fails the run rather than silently keeping
only a convenient subset.

Exact equal sequences collapse to one sequence class with all source-combination
lineages attached. Thus an insertion/deletion pair and a substitution pattern
can describe one sequence and receive one fragment count. An empty caller also
contributes a reference-sequence hypothesis; it does not contribute a reference
genotype call. Enumeration does not enforce native zygosity or phase as a
probabilistic model, and does not infer which two classes belong to the sample.

## Extraction and fragment accounting

The largest windows are merged for extraction. `samtools view -M -L` creates a
regional BAM without realignment, duplicate remarking or quality recalibration.
The original and extracted reference/read-group headers must match. Ordered
decoded records from sequentially consumed original regional iterators must
match the extracted BAM's count and rendered-record SHA-256. This verifies the
regional extraction, not the original alignment's biological correctness or
genome-wide mapping uniqueness.

Each read must align a query base to both reference boundary positions. The
helper reconstructs the exact interior query sequence using its CIGAR, including
internal insertions and deletions. It rejects reference skips in the window,
edge insertions, unaligned boundaries, missing qualities, malformed CIGAR,
internal clips and nonconcrete extracted sequence. Terminal clips can remain
when both window boundaries are genuinely spanned.

Evidence uses raw base qualities with no BAQ. It excludes flags `0xF04` and
unknown MQ255, and applies the requested MQ and minimum BQ across the entire
extracted sequence. Eligible improper pairs are retained. No mate lookup or
joining of two individually nonspanning reads manufactures a joint span.

All observations remain in the evidence. Eligible observations are grouped by
original `(RG, QNAME)`:

- Equal exact sequences from agreeing mates count once for the fragment.
- Different eligible sequences produce a conflict, even when both are unmatched
  by every hypothesis. They are not collapsed into generic agreement on “other.”
- A fragment with no eligible spanning observation remains explicit.
- A fragment with one eligible observation can contribute even when another
  observation was rejected; the rejected observation and reason remain visible.

Read-group, strand, CIGAR, MQ, raw-BQ and observation/fragment denominators are
retained. Counts are alignment-identifier groups, not cell counts, validated
independent DNA-molecule counts, or allele fractions corrected for sampling and
mapping bias. Flank and quality views reuse reads and must not be added. Increasing
a threshold can remove one disagreeing mate, so a particular matched-class
fragment count need not decrease monotonically even though observation
eligibility becomes stricter.

## Constructed sequence annotation

Each non-reference sequence class becomes a parsimonious replacement in a
separate, explicitly labelled **sites-only** VCF. It has no sample columns and no
GT. Reference classes have no variant record. The workflow performs reference-
checked normalization, sorts the normalized records, then runs `bcftools csq
--local-csq --samples -` against the pinned GFF.

Every expected constructed record ID must survive exactly once. Before/after
annotation payloads must match after removing only INFO/BCSQ. Each normalized
constructed allele is also applied back to the reference over the union interval
and must reproduce its complete intended sequence. These checks protect sequence
identity and annotation bookkeeping; they do not validate the predicted effect.
No independent annotation engine is used to validate the bcftools consequences.

A consequence describes one constructed local replacement conditional on that
sequence. It is not a newly called variant, a proven joint cellular haplotype,
or demonstrated loss of function. Local frame restoration does not imply normal
protein function. Missing BCSQ is not a benign classification. Consequences still
depend on transcript relevance, annotation completeness, splice context, and
changes outside the comparison window. A replacement spanning unchanged sequence
or exon boundaries needs explicit transcript-level review before being treated
as equivalent to separate biological mutation events.

## Boundary and interpretation limits

The hypothesis universe is deliberately bounded. Wider windows can introduce
background alleles absent from a narrower view; they also require more of a read
to span with sufficient quality. Changed counts across windows cannot be assigned
to biological phase without inspecting both effects.

There are two distinct boundary ledgers. Full regional native context preserves
reported-end overlaps, including nonconcrete records. The per-window edit list
uses literal REF-span overlap. A symbolic END-spanning record anchored outside a
window can therefore be present in regional context without appearing in the
window's omitted-allele list. That list alone is not a structural-variant census.

An insertion anchored at the final reference base can fit the literal REF rule
and enter the hypothesis family, while a read carrying that insertion exactly
at the window edge is rejected by the read-boundary rule. Such zero support is
geometry-censored. Use an explicitly wider reviewed window to examine it; do not
reinterpret absence as reference evidence. The same caution applies to crossing
deletions and insufficient spans.

For reuse with distinct windows on the same chromosome, an alignment record that
spans multiple unmerged extraction intervals can appear twice in independent
original-region fetches but once in the `-M` extraction. The exact payload check
then fails conservatively. It must not be bypassed or described as completed
verification. The present selected targets use separate chromosomes.

Read support remains conditional on the existing alignment. Nearby paralogs,
repeats, alignment artifacts, original library duplication and locus assignment
are not resolved by local sequence matching. Shared counts or equal caller
representations provide no parental origin, phase between distant sites, clinical
diagnosis, EPCR or treatment evidence. The review is post hoc and is not a blind
discovery experiment or genome-wide exclusion analysis.

## Qualification and reproducibility

The retained `qualification-v1/controls.json` receipt, copied into the
[source-only qualification receipt](../source-snapshots/joint-haplotype-controls-20260920/controls.json), records **20 tests passed,
none skipped**, on Mac Studio 512S2 at `2026-09-20T14:33:05.462578+00:00`.
It records Python 3.12.14, pysam 0.24.1 and macOS 26.4.1 on arm64. The 17
[helper controls](../tests/test_joint_haplotype_helpers.py) and three
[native integration controls](../tests/test_joint_haplotype_workflow.py) cover:

- Multi-ALT exclusivity, uncalled lineage, caller separation, equivalent sequence
  collapse, overlap rejection, hypothesis limits and replacement anchoring.
- Native CIGAR/quality handling, clips and boundaries, MQ255/flags, missing
  identities/qualities, lane namespaces, strands, unmatched reads and conflicts.
- A fictional insertion/deletion versus four-SNP equivalence with three read
  observations contributing two fragments, including an agreeing mate pair.
- Native raw VCF/BCF payload preservation, CRAM-to-BAM decoded extraction
  comparison, sites-only normalization/annotation and original-input immutability.
- Long upstream REF overlap, crossing-window context, and raw/read hash or
  complete-source payload/count refusals.

Qualified source identities are:

| File | SHA-256 |
| --- | --- |
| `review_joint_haplotypes.py` | `39e2e566bcc1dc0938b1b4318e8a3aeb045c51fe175245d9dfaa110206732416` |
| `joint_haplotype_helpers.py` | `0edf1f484bbe2b5009b31f5f0205a96e617b11d1d0677c6e3bb0b59eed06363e` |
| `native_coding_records.py` | `c395dd130a2577e7b25ba0ac62d600d7cb64c9240f3eb218536eff22c5c44cfe` |
| `test_joint_haplotype_helpers.py` | `c6bea90b77a6ae498c4123bea1f06c94ec78572061f155ee01b894403ae14276` |
| `test_joint_haplotype_workflow.py` | `bfb6d50274db3f7d0df72330507aa5204869c69472049e27b842326c417e3f2f` |
| `joint-haplotype-review.json` configuration | `ff131f40eaf74aaf9aadb7c11a35df78c11a2a9c9c21c7964e3a6cd9d093cfb9` |

These are fictional method controls, not a biological sensitivity, precision,
phase-accuracy or clinical validation study. A successful run still requires
checking its completed status/summary linkage, output hashes, native context,
omissions, sequence lineages, fragment denominators and constructed annotation.

To repeat the controls in the installed native Studio environment, set
`MVA_PROJECT_ROOT` to its project root and run the two test files with that
environment's Python:

```sh
python -m unittest discover -s workflow/tests -p 'test_joint_haplotype*.py' -v
```

The private source configuration pins the intended production inputs and target
selection. Production case results require their own retained output directory,
provenance and review, independently of the fixture-qualification receipt.
Reusing this method on another dataset requires its own explicit configuration
and provenance; the fictional controls are not evidence about that dataset.
