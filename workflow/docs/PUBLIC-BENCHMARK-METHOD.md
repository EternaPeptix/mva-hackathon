# Public small-variant benchmark: what the numbers measure

The existing comparison already used **RTG Tools 3.13 `vcfeval`**, not an ad hoc exact-position comparator. It evaluated public HG002 NovaSeq PCR-free **nominal 35×**, prealigned/deduplicated chromosome 20 reads against GIAB **v4.2.1** GRCh38 small-variant truth and its official confident-region BED. Nominal 35× is the public input's label; a measured mean-depth or depth-distribution analysis is not supplied by these results. The completed pilot records 19,430,746 mapped chr20 alignment records, which is not itself a coverage estimate.

The BED contains 10,192 intervals spanning 56,916,239 of 64,444,167 chr20 bases (88.3187%). Full chr20 VCF context is retained and the BED is passed through `--evaluation-regions`, allowing boundary-crossing matching. Of 85,951 prepared truth records, 82,589 are scored: 71,333 SNV, 11,227 indel and 29 mixed records. The prior scope audit found all 3,362 unscored truth records outside the BED by REF-span overlap; one scored boundary deletion has its anchor outside the BED. This is not a whole-genome result or a claim of accuracy in excluded regions.

The frozen command uses `--sample HG002,HG002 --sample-ploidy 2 --vcf-score-field QUAL --at-score 0 --output-mode split --threads 2`; it does not request ploidy squashing, decomposition or reference-overlap mode. The fixed-score counts agree with the native unthresholded summary. The callers used raw discovery without BQSR or an added post-call quality filter. [Recorded commands and outputs](../../results/germline-public-comparison-20260919.json) preserve the exact options.

RTG compares genotype-derived haplotype sequences and can match equivalent representations; both alleles must agree under this diploid comparison. Consequently, matched baseline and call record counts can differ. Recall uses baseline TP/(baseline TP+FN); precision uses call TP/(call TP+FP). These are not one-to-one error counts. See the [RTG command reference](https://realtimegenomics.github.io/rtg-tools/rtg_command_reference.html).

| Arm | Caller | Record class | Baseline TP / truth | Recall | Call TP / scored calls | Precision |
|---|---|---|---:|---:|---:|---:|
| Original split-normalized | gatk | SNV | 70,887 / 71,333 | 99.3748% | 70,915 / 71,664 | 98.9548% |
| Original split-normalized | gatk | INDEL | 10,934 / 11,227 | 97.3902% | 11,426 / 11,900 | 96.0168% |
| Original split-normalized | bcftools | SNV | 70,876 / 71,333 | 99.3593% | 70,951 / 71,703 | 98.9512% |
| Original split-normalized | bcftools | INDEL | 10,274 / 11,227 | 91.5115% | 10,693 / 11,727 | 91.1827% |
| Retained-genotype sensitivity | gatk | SNV | 70,909 / 71,333 | 99.4056% | 70,909 / 71,614 | 99.0156% |
| Retained-genotype sensitivity | gatk | INDEL | 11,130 / 11,227 | 99.1360% | 11,129 / 11,208 | 99.2951% |
| Retained-genotype sensitivity | bcftools | SNV | 70,898 / 71,333 | 99.3902% | 70,973 / 71,681 | 99.0123% |
| Retained-genotype sensitivity | bcftools | INDEL | 10,452 / 11,227 | 93.0970% | 10,442 / 11,107 | 94.0128% |

These classes describe each full REF/all-ALT VCF record, not an atomic SNP/INDEL decomposition. The 29 mixed truth records remain separate: original and retained arms match 28 for GATK and 24 for bcftools. The complete class table also retains two GATK call-side non-sequence records matched in the retained-genotype arm; these are not silently reassigned to the SNV or indel table. Class-specific baseline and call totals therefore need not be equal.

The retained-genotype arm is explicitly **post-observation sensitivity analysis** of the same discoveries, not another sample, an independently tuned validation set, or replacement of the original result. Relative to the original, GATK gains 219 truth matches and loses one; bcftools gains 204 and loses four. Its improved indel result shows that representation and genotype retention materially affect the original aggregate. [Reconciliation](../../results/germline-representation-reconciliation-20260919.json) and [all four-arm counts](../../results/public-small-variant-review-20260921.json) preserve the distinction.

A new read-only check on the Mac Studio independently parsed the 16 saved RTG TP/FN/FP VCFs, reproduced every record-class count and ratio, verified the BED geometry, and freshly matched 30 file hashes. It did not rerun calling or matching. A second `vcfeval` run would be a reproducibility test, not an orthogonal comparison engine; no new hap.py/xcmp validation is claimed. The public reads' upstream alignment/duplicate marking differs from the private workflow, so this benchmark does not validate raw-read alignment, patient callability, clinical sensitivity or causal probability. [Review receipt](../../results/public-small-variant-review-20260921.json).
