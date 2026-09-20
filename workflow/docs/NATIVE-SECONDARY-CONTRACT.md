# Native secondary-finding adapter: proposed contract

Preparation only, 20 September 2026. No case screen is authorised by this artifact.
The implementation and fictional controls live separately from all completed
results. Runtime is Mac Studio 512S2; no caller, aligner or internet retrieval.

## Reuse and source identity

Pin exact unchanged copies of `screen_secondary_findings.py` (SHA256
`34b7a4bebfb48d1c460fd579d862f18ff802f9ed56c8a39a5475882bc3e1802d`)
and `native_coding_records.py` (SHA256
`c395dd130a2577e7b25ba0ac62d600d7cb64c9240f3eb218536eff22c5c44cfe`).
Reuse the former's germline ClinVar strata, concrete-allele test, exact keys,
panel-count check and reference-checked normalization stream. Reuse the latter's
source context, original ALT recovery and native allele call. Do not invoke either
production `run`, the coding consequence/model selection, or the historical quality
eligibility gate. Verify these module bytes before importing them.

## Frozen input contract

A separately SHA-pinned JSON configuration names the exact input paths, hashes and
sizes. All JSON is strict (no nonfinite constants or duplicate keys); required
booleans and integers are typed. The real-input draft is deliberately not ready.
Inputs include reference/index, bcftools, dated ClinVar VCF/index, 84-row ClinGen
panel, full-autosome handoff summary and native-verification receipt, fixed-X/Y
terminal status and independent review/verification receipts. Hashes are checked
before interpretation; wrong resolved paths cannot be replaced by same-byte files.
All input metadata is checked again at completion. No old input is indexed or edited.

Each caller is separate. The autosomal native-verification `files[]` must bind the
actual complete handoff summary and each original source, with exact paths/hashes/
sizes. The handoff must prove all 22 autosomes once, no pending autosomes, complete
split replay, and original record counts/payload fingerprints. Each of the 22
original sources is replayed in full against that census, including zero-record
sources. The X/Y source paths derive exclusively from the fixed run's directory and
its named `gatk.original.vcf.gz`/`bcftools.original.vcf.gz` outputs and pinned indexes.
The fixed independent review must bind the status hash, original record counts and
output hash map; its verification receipt must bind the actual review summary path.
The legacy verification schema has no status field. Require aware terminal times,
matching executed-source digests, all fixed outputs counted and mapped, and explicit
original-preservation/bcftools-reuse evidence from the linked completed review.
A header or absence of variant rows is never accepted as completed source coverage.
Repeated source paths or duplicate caller/scope membership fail closed.

The fixed X/Y configuration must cover X/Y without interval gaps/overlaps; X PAR is
diploid, remaining X/Y haploid according to the pinned model. Required Y PAR masks
must be N in the same reference. Native record REF spans must remain within one
ploidy interval. Informative or partially called GT must have the declared shape. Completely
unavailable GT, including absent FORMAT GT, literal `.`, and `./.`, remains unknown
without asserting haploidy or rejecting its missing encoding. Encoded slot count is
retained separately from observed ploidy. Missing and partial states stay explicit. Autosomes use diploid shape checks separately. The adapter never applies
the autosomal handoff's scope assumptions to X/Y or substitutes default whole-run
X/Y calls. This computational model does not establish demographics or karyotype.

## Normalization and full native provenance

Stream every original record into a new ordinal-tagged copy, preserving and
replaying all original rendered records after removing only the added ordinal.
Normalize that copy with pinned bcftools `norm --check-ref e -m -any`, reference and
old-record tag. No atomization, quality exclusion, AF threshold, coding restriction
or consequence inference is applied. Reserve tags and refuse preexisting collisions.
For every normalized exact panel match, recover the original by source ordinal and
ALT index; verify literal reference/edit equivalence independently of the old tag.

Every occurrence retains caller, sample, source SHA, ordinal, source namespace,
full original rendered record/FORMAT, original ALT list/index, GT including ploidy/
phase flag, AD/PL/DP/GQ and physical-phase fields. Identical duplicate originals
remain distinct occurrences. Repeated normalized identities are grouped for unique
allele accounting, with all source lineages retained and duplicate/conflicting
native evidence explicitly visible. No duplicate row is silently dropped, assigned
consensus GT or counted as another unique allele. Multiple source ALT lineages that
converge on one key remain separately identifiable. Per-caller groups do not imply
cross-caller genotype agreement.

## Quality and exact-match reporting

Retain all exact concrete matches, including ALT unobserved in GT and partial or
absent GT. Report called-ALT-present, fully-called-ALT-present and uncalled/partial
categories separately. Retain explicit PASS, unfiltered and named filters; DP/GQ
FORMAT absence, present-missing, explicit zero and descriptive historical bins
remain distinct. Do not synthesize PASS, impute GQ from PL, or apply a combined
quality-pass eligibility gate. Haploid calls remain eligible for exact-allele
inventory; no diploid recessive model or parental phase is inferred.

Public selection retains the frozen 84 genes and every concrete exact ClinVar
allele with the existing germline P/LP signal logic, regardless of coding annotation.
Clean aggregate P/LP, conflicts/mixed assertions, zero-star and haplotype/genotype-only
records remain separate. A multiallelic public original enters a distinct
`multiallelic_assertion_mapping_review` stratum even when its record-wide classification
is P/LP. Original and normalized CLNALLE, ALT lineage and original record stratum are
retained; the adapter does not claim proven per-ALT assertion mapping. Original
record strata come from the original INFO, while normalized-record strata are
separate and may be null. For example, an original Number=A P/Benign assertion
keeps both split ALT lineages in the ambiguous review stratum, including the
normalized benign ALT without a P/LP signal. Biallelic source and normalized
assertions must agree. Selected
biallelic/multiallelic source counts make this scope visible. Somatic assertions do
not become germline classifications;
unknown review labels fail for review. Preserve every full original ClinVar rendered record/INFO context, original ordinal
and ALT index, normalized record and normalization lineage multiplicity, alongside
the exact dated resources. Reportability remains pending condition/variant/
inheritance review. No new clinical classification is performed.

Output is a new fail-closed directory: source scope/censuses, all exact occurrence
rows, unique identity groups with their lineages, caller-separated counts/quality
axes/ClinVar strata, full input/helper/command/output pins and terminal status.
No exact matching occurrence and no exact called-ALT identity are separate results;
uncalled/reference-only matching records remain visible even when called-ALT count
is zero. Zero called-ALT identities means only none within the explicit sources and
dated resources. It does not establish a clinical negative, absent-record 0/0,
callability, novel/SV/repeat/mosaic absence, or coverage of the 84-gene panel. The
17-gene CDS-depth run remains a separate scope.

## Fictional qualification before any case proposal

Real bcftools normalization controls must cover: raw FILTER '.', missing GQ and
complete PL/custom FORMAT preservation; multi-ALT 1/2, haploid ALT and partial/missing
GT; duplicate originals and convergent normalized keys without multiplicity loss;
noncoding exact matches; conflict/haplotype/somatic distinctions; absent source
chromosome versus complete zero-record chromosome; stale hashes and wrong-path
same-byte proofs; X/Y interval/GT/boundary/mask failures; unchanged inputs and refusal
of existing output directories. Root reviews this contract before qualification.
Any actual case execution requires a later root review of exact qualified source,
complete real-input pins and output scope.
