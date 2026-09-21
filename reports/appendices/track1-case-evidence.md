# Track 1: case evidence and reconstruction map

**21 September 2026.** This appendix connects the submitted BUB1B hypothesis to
existing results. The [selected numerical evidence](../../results/track1-selected-case-evidence.json)
contains candidate-level and gene-level aggregates recovered from completed
analyses, with original summary names. It is a curated projection, not a new run
or an export of restricted case callsets. The [interpretation](causal-interpretation.md)
explains the classification and subjective EPCR; the numbers below do not supply
a calibrated causal probability.

## Ranking and predictor sensitivity

The 19 September genotype-assisted rankings used the **organizer-supplied VCF**,
not our later four-lane GATK/bcftools callsets. The primary phenotype input contains
the seven supplied subject features in the [LIRICAL configuration](../../workflow/configs/lirical-subject.json).
No unrecorded MVA-typical feature or negative finding was added. The previously
nominated candidate was known. No completed ranking rerun on the final native
callsets is claimed.

| Completed comparison | Retained result | What it shows |
|---|---|---|
| Exomiser, seven subject features and primary predictor set | BUB1B AD rank 1, score 0.5464; **AR rank 2, score 0.4877**, both nominated alleles contribute | Rank is gene × inheritance model; the AD entry does not establish dominant MVA. |
| Exomiser, REVEL only, `keepNonPathogenic: false` | AD rank 1, score 0.5464; **AR entry absent** | The missense is removed. Gene prominence is more stable than support for the exact pair. |
| Exomiser, REVEL only, `keepNonPathogenic: true` | AD rank 1, score 0.5464; **AR rank 3, score 0.1806**, both alleles retained | Pair retention depends on the declared predictor/filter policy. |
| Exomiser, added parental-history term | AD rank 1, score 0.5860; AR rank 2, score 0.5565 | Exploratory HP:0200067 addition; this is not a feature of the child's primary input. |
| LIRICAL, seven features, phenotype only | MVA1 **rank 4/8,525**, model probability 0.0178585 | Bartter 2, Bartter 1 and Costello rank above MVA1. |
| LIRICAL, phenotype only, omit rhabdomyosarcoma | MVA1 **rank 38/8,525**, model probability 0.0000981050 | The cancer phenotype materially affects the disease ranking. |
| LIRICAL, rhabdomyosarcoma alone | MVA1 **rank 10/8,525**, model probability 0.0212790 | One feature is not sufficient for the primary phenotype ranking. |
| LIRICAL, full global JSON, phenotype plus supplied VCF | MVA1 **rank 1/8,525**, model probability 0.9575393 | Genotype integration strongly promotes MVA1 within this model; this is not EPCR or observed trans phase. |

Exomiser is **15.1.0**, with **2512 hg38 and phenotype bundles**. The
[primary configuration](../../workflow/configs/proband01-exomiser.yml) uses REVEL,
MVP, AlphaMissense and SpliceAI; the [strict](../../workflow/configs/revel-only-strict.yml)
and [permissive](../../workflow/configs/revel-only-permissive.yml) comparisons use
REVEL only. LIRICAL is **2.4.1**, HPO **2026-09-01**, with the Ensembl transcript
database and shared Exomiser allele databases. The historical report dates
phenotype annotations to 2026-09-02, but that raw date header was not independently
recovered. The numerical evidence file retains the actual public-resource hashes.

LIRICAL ranks above come from its complete JSON, not the filtered native TSV
(only 2,168 diseases in the genotype-assisted run). Its parental-history run
failed strict validation and supplies no successful result. Model probabilities
are not a calibrated, mutually exclusive distribution of case diagnoses.

For **N1002K**, the retained database values are **REVEL 0.472** and
**AlphaMissense 0.9229**. Historical provenance identifies REVEL v1.3 and the primary
AlphaMissense hg38 release; an exact AlphaMissense release date is not retained.
The older ranking transcript is ENST00000287598.10, while current sequence
reconstruction uses .11. The adopted classification policy treats REVEL 0.472 as
indeterminate: **neither PP3 nor BP4**. AlphaMissense is a disclosed alternative
predictor, not a second independent observation or a reason to choose the more
favorable classification. LIRICAL's maximum-predictor rule can display
pathogenicity 1.0 from another predictor; that is not measured variant function.

## Final alignment and native-call evidence

Here, **native** means each caller's original records, including its actual
genotypes, filters and missing fields before split-allele recoding. Agreement
between callers or rankers reuses the same biological sample.

| Question | Selected completed result | Scope |
|---|---|---|
| Are both nominated alleles supported? | Leu737Ter: **21 REF / 26 ALT fragments**; N1002K: **15 / 12**, both with ALT on both strands | Final four-lane alignment; MQ≥20, BAQ-adjusted BQ≥20, flags 0xF04 excluded, orphan pairs excluded; one observation per read-group/query-name pair. These are not caller AD values. |
| Is there another emitted coding/splice candidate in BUB1B? | Whole-chr15 GATK/bcftools outputs have **17/16 records** in the gene-plus-5 kb window and **14/13** in the gene; each has only the **two nominated exon/CDS overlaps** and **zero additional screened splice overlaps** | GRCh38 15:40155984–40226137, inclusive; Ensembl 116, MANE ENST00000287598.11. Screen includes canonical two intronic bases and a broader three-exonic/eight-intronic-base edge; not every splice mechanism. Native FILTER is unfiltered `.`. |
| Is a coding alternative present in the ten comparison genes? | **Zero AF <1% protein-altering/canonical-splice alleles**, with caller agreement; three common alleles across CENATAC, SLF2 and SMC5 | Completed chr1–22 inventories, no predictor or FILTER/DP/GQ admission threshold. Exact gene counts and retained AFs are in the numerical evidence and [interpretation table](causal-interpretation.md#what-the-completed-analysis-supports). |
| Could coding depth explain the negative inventory? | BUB1B **3,549 CDS-union bases, minimum 24×** at both MQ20/BQ20 and MQ60/BQ20; all ten comparison genes have **zero zero-depth coding bases** | SMC5 still has 15/17 bases below 20× and BUB3 4/5, at MQ20/MQ60 respectively. Overlapping mates suppressed, flags 0xF04 excluded, no BAQ. Depth is not clinical callability. |
| Was trans phase observed? | **No.** Vendor maximum block 288 bp; earlier two-lane 100 kb review maximum 2,327 bp; final four-lane bounded review leaves the pair unknown in every run | Candidate distance 10,911 bp. Final review is 15:40204701–40225612; maximum local block 77 bp in the all-variant bcftools run, none in the other three runs; intervening marker gap 6,769 bp. WhatsHap 2.8 uses MQ≥20 and its quality model, without a hard BQ20 filter. |

The inventory AF is the maximum valid positive-denominator AC/AN among retained
gnomAD exome/genome sources in Exomiser 2512. Missing frequency is unresolved, not
zero; none of the five inventory alleles across BUB1B plus these ten genes has
missing frequency. The table uses **AF <0.01**; the native model sensitivity uses
**AF ≤0.01**. No relevant allele equals that boundary. This discovery threshold
is distinct from the candidate-specific PM2 judgment.

## Authorized reconstruction map

Obtain the controlled FASTQs, supplied VCF and clinical inputs through the
organizers; keep them and generated case outputs outside Git. Provision the
recorded reference and tool versions from the [environment/reproduction guide](../../REPRODUCIBILITY.md).
The public HG002 replay validates a separate component and does not reproduce the
case. The following maps the actual evidence stages to retained source:

| Stage and required inputs | Retained workflow/configuration | Output to inspect |
|---|---|---|
| Supplied VCF + seven HPO features + Exomiser 2512 resources | [prioritize.py](../../workflow/scripts/prioritize.py), primary/strict/permissive configurations linked above | Gene × inheritance scores, contributing alleles and sensitivity summaries; use the supplied VCF to reconstruct the reported ranks. |
| Same phenotype + LIRICAL resources; supplied VCF and Exomiser databases for genotype mode | [run_lirical.py](../../workflow/scripts/run_lirical.py), [prepare_lirical_data.py](../../workflow/scripts/prepare_lirical_data.py), subject configuration | Complete JSON disease results; genotype comparison uses `--global-analysis`, without `--show-all-diseases` or `--strict-genotypes`. |
| Eight checksum-verified FASTQs + matching GRCh38 reference | [01_align.sh](../../workflow/scripts/01_align.sh), [call_genome.py](../../workflow/scripts/call_genome.py), [CRAM compatibility helper](../../workflow/scripts/prepare_gatk_compatible_cram.py) | Completed four-lane alignment/state, GATK and bcftools contig outputs; retain original native genotypes. |
| Completed discovery outputs + Ensembl 116 + retained frequency database | [handoff audit](../../workflow/scripts/audit_discovery_handoff.py), [prepare_native_coding_review.py](../../workflow/scripts/prepare_native_coding_review.py), [native policy](../../workflow/configs/native-coding-review.json) | Each caller's `summary.json` and `coding-alleles.jsonl.gz`; group by the 11 stated genes and apply AF <0.01 without a new quality or predictor cutoff. |
| Final alignment/state + selected targets | [audit_final_alignment_targets.py](../../workflow/scripts/audit_final_alignment_targets.py), [target configuration](../../workflow/configs/bub1b-noncoding-review.json) | Candidate fragment summaries and regional callsets; source configuration's historical lane wording is superseded by the runner's required completed four-lane state. |
| Final regional callsets/alignment + exact reference | [run_read_phasing.py](../../workflow/scripts/run_read_phasing.py), [four-lane configuration](../../workflow/configs/bub1b-four-lane-phasing.json) | Four phase summaries and pair relation; do not substitute the earlier two-lane configuration. |
| Final alignment + Ensembl CDS targets | [measure_coding_depth.py](../../workflow/scripts/measure_coding_depth.py), [coding-depth configuration](../../workflow/configs/coding-depth-v2.json) | Per-gene depth summaries in the declared MQ/BQ modes. |
| Alternative sequence or copy-number hypotheses | [joint haplotype review](../../workflow/scripts/review_joint_haplotypes.py), [joint configuration](../../workflow/configs/joint-haplotype-review.json), [copy-number review](../../workflow/scripts/review_four_lane_copy_number.py) | The curated interpretation records ADAMTS1's net in-frame change, FAM136A's combined missense and retraction of the apparent chromosome-wide mosaic signal. Full case result ledgers remain restricted. |

These entry points require local paths and retained prerequisites, not just a
fresh clone. In particular, the full original BUB1B window census/geometry ledger,
historical case execution state, per-read evidence and genome-scale outputs are
not exported. The selected aggregate makes the claim inspectable; it does not
make those omitted results independently replayable here. Reconstruct the window
comparison from native chr15 records, the stated coordinates and transcript/edge
rules, rather than treating the coding inventory as an identical window census.

Example for a provisioned authorized workspace (paths are placeholders):

```sh
python3 workflow/scripts/prioritize.py \
  --vcf /authorized/organizer-supplied.vcf.gz \
  --exomiser /tools/exomiser-cli-15.1.0 --data /resources/exomiser \
  --analysis workflow/configs/proband01-exomiser.yml \
  --data-version 2512 --output /authorized/new-exomiser-primary

python3 workflow/scripts/run_lirical.py \
  --root /authorized/workspace --data /resources/lirical-compatible \
  --config workflow/configs/lirical-subject.json \
  --vcf /authorized/organizer-supplied.vcf.gz \
  --exomiser-data /resources/exomiser/2512_hg38 --global-analysis \
  --output /authorized/new-lirical-global
```

The LIRICAL root must contain the Java and 2.4.1 JAR paths required by the runner;
its prepared compatible resource directory must match the retained identities.
These are reconstruction examples, not newly executed runs. No fresh end-to-end
case runtime, Linux case replay, segregation result or functional assay is claimed.
