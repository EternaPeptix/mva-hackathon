# Evidence and component index

The tables distinguish public human benchmark data from synthetic controls and
runtime metadata. Counts belong to their individual receipts and are **not summed**
into a validation score. Historical tests qualify the recorded source/input bytes;
this documentation revision does not claim that every test was rerun.

A receipt records commands, versions, hashes and outcomes. “Native” means the actual
installed tool ran; “fictional” in historical files means synthetic test data.

## Track 1: small variants

| Component | Data type | Control count / scope | Receipt or guide |
|---|---|---|---|
| Portable public replay | Public HG002 / synthetic comparator inputs | Linux x86-64 emulated: 12 controls and 100 kb calling/comparison passed; full Mac replay pending | [Current status](public-benchmark-replay-20260921.json); [Command](../workflow/public-benchmark/README.md) |
| Small-variant primary comparison | HG002 chr20 | 82,589 scored truth records | [Record](../results/germline-public-comparison-20260919.json) |
| Truth regions and input identities | Public HG002/GIAB | Not a control suite | [Downloads](../results/germline-public-inputs-20260919.json); [Reference/BED](../results/germline-public-reference-20260919-v2.json) |
| Native GATK/bcftools execution | HG002 chr20 | Both original callers | [Record](../results/germline-public-pilot-complete-20260919.json) |
| Post-observation multiallelic sensitivity | HG002 chr20 | Same scored truth universe | [Record](../results/germline-retained-multiallelic-comparison-20260919.json); [Reconciliation](../results/germline-representation-reconciliation-20260919.json); [All losses](../results/germline-representation-loss-review-20260919.json) |
| Comparator behavior | Synthetic | 12 cases | [Record](../results/germline-public-comparator-controls-20260919.json) |
| Independent comparison audit | Synthetic / HG002 | 6 synthetic checks; actual data reviewed separately | [Controls](../results/germline-public-audit-controls-20260919.json); [HG002 audit](../results/germline-public-comparison-audit-20260919.json) |
| Multiallelic representation mechanisms | Synthetic | 6 cases | [Record](../results/germline-multiallelic-controls-20260919.json) |
| Metadata and interval checks | Synthetic / saved public metadata | 8 separate checks | [Record](../results/germline-representation-validation-20260919.json) |

## Structural variants

The public structural comparison uses the same 802 HG002 chr20 truth alleles in
both caller arms. Shared misses include repeat-flagged deletions and larger
insertions. Annotation stratification describes the gap but cannot attribute every
miss to short-read technology, truth composition or caller behavior.

| Component | Data type | Control count / scope | Receipt or guide |
|---|---|---|---|
| DELLY and SvABA comparison replay | HG002 chr20 | 802 truth alleles; saved calls reused | [Record](../results/runtime-published-source-public-sv-20260920.json) |
| Structural comparator | Synthetic | 6 cases | [Record](../results/runtime-structural-comparator-20260919.json) |
| Truth-annotation residual analyzer | Synthetic; intended for saved HG002 decisions | 29 producer controls | [Record](../workflow/source-snapshots/public-sv-residual-review-20260920/producer-preparation/qualification-v1/controls.json) |
| Independent residual reader/helpers | Synthetic | 9 helper checks; not full-checker qualification | [Record](../workflow/source-snapshots/public-sv-residual-review-20260920/independent-execution-review/qualification-v2/controls.json) |
| Structural residual method and observed summary | HG002 chr20 | No new caller | [Reproduction](../workflow/source-snapshots/public-sv-residual-review-20260920/README.md); [Results](../reports/appendices/public-sv-residuals.md) |

## Track 1: software controls

These checks exercise representation, boundaries, preservation and refusal of
invalid inputs. Passing them is evidence about the tested software behavior,
not an estimate of clinical sensitivity.

| Component | Data type | Control count / scope | Receipt or guide |
|---|---|---|---|
| CRAM compatibility and record preservation | Synthetic | 195 planted records; mutation challenges, not a numbered unit suite | [Record](../results/gatk-cram-compatibility-native-control-20260920.json) |
| Structural CRAM 3.0 compatibility | Synthetic | 6 compared native calls; full record equality | [Record](../results/structural-cram30-control-20260920.json) |
| Copy-number representation | Synthetic | 9 | [Record](../results/cnv-context-method-controls-20260920.json) |
| Structural/CNV selection | Synthetic | 2 | [Record](../results/structural-cnv-context-controls-20260920.json) |
| Tandem junction/transcript geometry | Synthetic | 6 | [Record](../results/tandem-junction-method-controls-20260920.json) |
| Nearby-pair phase/mapping/mates | Synthetic | 6 | [Record](../results/paralog-pair-method-controls-20260920.json) |
| SNV clusters/fragments | Synthetic | 5 | [Record](../results/snv-cluster-method-controls-20260920.json) |
| Local sequence comparison | Synthetic | 4 | [Record](../results/public-clone-method-controls-20260920.json) |
| Target-base orientation | Synthetic | 2 | [Record](../results/cluster-target-projection-controls-20260920.json) |
| Four sequencing-lane indel review | Synthetic | 4 | [Record](../results/four-lane-indel-controls-20260920.json) |
| Indel window/orientation | Synthetic | 2 | [Record](../results/four-lane-indel-context-controls-20260920.json) |
| Fixed-region depth/missingness | Synthetic | 4 | [Record](../results/chromosome-depth-method-controls-20260920.json) |
| Regional-ploidy interval handling | Synthetic | 3 | [Record](../results/ploidy-region-method-controls-20260920-v2.json) |
| Diploid/haploid native calls | Synthetic | 3 planted alleles × 2 callers | [Record](../results/ploidy-calling-native-controls-20260920.json) |
| Ploidy-boundary genotype preservation | Synthetic | 4 | [Record](../results/ploidy-boundary-method-controls-20260920.json) |
| Ploidy-boundary native recovery | Synthetic | 12 planted SNVs; 24 fixed-ploidy caller/allele checks | [Record](../results/ploidy-boundary-native-controls-20260920.json) |
| Mitochondrial geometry/read preservation | Synthetic | 7 | [Record](../results/mitochondrial-method-controls-20260920.json) |
| Mitochondrial original-allele review | Synthetic | 2 | [Record](../results/mitochondrial-review-method-controls-20260920.json) |
| Circular-read calling | Synthetic | 4 planted alleles | [Record](../results/mitochondrial-native-control-20260920.json) |
| Mitochondrial read-context evidence | Synthetic | 7 | [Record](../results/mitochondrial-context-method-controls-20260920.json) |
| Mixture caller/importer contract | Synthetic | Separate scenario families; see record | [Record](../results/mutserve-contract-controls-20260920.json) |
| Mixture selection/allele rank | Synthetic | Separate scenario families; see record | [Record](../results/mitochondrial-mixture-method-controls-20260920.json) |
| Changed-call read review | Synthetic | 3 | [Record](../workflow/source-snapshots/changed-call-read-controls-20260920/controls.json) |
| Completed-autosome handoff | Synthetic | 5 | [Record](../workflow/source-snapshots/discovery-handoff-controls-20260920/controls.json) |
| Coding inventory and descriptive models | Synthetic | 18 | [Record](../workflow/source-snapshots/native-coding-controls-20260920/controls.json) |
| Joint sequence and fragment evidence | Synthetic | 20 | [Record](../workflow/source-snapshots/joint-haplotype-controls-20260920/controls.json) |
| Full-query clipping/remapping | Synthetic | 12 | [Record](../workflow/source-snapshots/clipped-read-controls-20260920/controls.json) |
| Coding-region depth | Synthetic | 10 | [Record](../workflow/source-snapshots/coding-depth-controls-20260920/controls.json) |
| Normalized-representation overlay | Synthetic | 21 | [Record](../workflow/source-snapshots/adamts1-normalized-overlay-20260920/controls.json) |
| Independent overlay checker | Synthetic | 23 | [Record](../workflow/source-snapshots/adamts1-normalized-overlay-20260920/independent-verification/controls.json) |
| Secondary inventory producer | Synthetic | 35 | [Record](../workflow/source-snapshots/native-secondary-controls-20260920/controls.json) |
| Independent secondary checker | Synthetic | 36 | [Record](../workflow/source-snapshots/native-secondary-independent-20260920/controls.json) |
| Whole-reference record preservation | Synthetic | 6 | [Record](../workflow/source-snapshots/whole-reference-completion-20260920/qualification-native.json) |
| Public secondary-panel target construction | Synthetic | 16 | [Record](../workflow/source-snapshots/secondary-panel-depth-20260920/public-target-qualification/controls.json) |
| Depth producer input guards | Synthetic | 21 | [Record](../workflow/source-snapshots/secondary-panel-depth-20260920/producer-guard-qualification/controls.json) |
| Independent depth checker | Synthetic | 40 | [Record](../workflow/source-snapshots/secondary-panel-depth-20260920/independent-qualification/controls.json) |
| Independent depth input guards | Synthetic | 8 | [Record](../workflow/source-snapshots/secondary-panel-depth-20260920/independent-guard-qualification/controls.json) |
| Public nuclear-panel transform | Synthetic | 25 | [Record](../workflow/source-snapshots/nuclear-mixture-20260920/public-transform-qualification/controls.json) |
| Nuclear-mixture native pilot | Synthetic | 25 | [Record](../workflow/source-snapshots/nuclear-mixture-20260920/fictional-pilot-qualification/controls.json) |
| Nuclear assay/saved-output wrapper | Synthetic | 34 | [Record](../workflow/source-snapshots/nuclear-mixture-20260920/assay-qualification/controls.json) |
| Saved-depth interval/transcript map | Synthetic | 39 | [Record](../workflow/source-snapshots/secondary-depth-gap-map-20260920/qualification-v1/controls-run/controls.json) |
| Depth appendix adaptation | Synthetic | 6 | [Record](../workflow/source-snapshots/secondary-depth-gap-map-20260920/appendix-adaptation/controls.json) |
| Independent interval-union helper | Synthetic | 7 examples; not full-checker qualification | [Record](../workflow/source-snapshots/secondary-depth-gap-map-20260920/independent-audit-adaptation/union-helper-example-controls.json) |

## Track 2 and pitch

| Component | Data type | Control count / scope | Receipt or guide |
|---|---|---|---|
| Prospective founder-cell ledger | Synthetic | 21 software controls; no experimental cells | [Record](../workflow/source-snapshots/founder-ledger-controls-20260920/controls.json); [Method](../workflow/docs/FOUNDER-LEDGER.md) |
| Assay capacity | Planning assumptions | Calculation; not an experiment or power estimate | [Record](../results/assay-capacity-planning-20260920.json) |
| Sonic narration assembly | Synthetic media/control inputs | 14 software tests | [Record](../workflow/source-snapshots/sonic36-20260920/assembly-controls.json) |
| Pitch reproduction | Source and production recipe | Not a biological validation | [Earlier local recipe](../workflow/docs/pitch-v2-reproduction.md); [Sonic source mapping](../workflow/docs/V27-SOURCE-REPRODUCTION.md) |

## Environments and public resources

| Component | Data type | Control count / scope | Receipt or guide |
|---|---|---|---|
| Runtime inventory | Software metadata | Original analysis environment | [Record](../workflow/environments/20260919/runtime.json) |
| Python packages built from published sources | Public software sources | 40 publisher wheels; 5 source builds | [Record](../results/runtime-published-source-build-20260920.json); [Source identities](../results/runtime-wheel-sources-20260919.json) |
| Fresh reconstructed environments | Synthetic/runtime | 22 reused checks across 7 environments; 6 extension cases | [Record](../results/runtime-published-source-validation-20260920.json) |
| Public CNV context resources | Public gnomAD/ClinGen | Acquisition identities; not case matches | [Record](../results/cnv-context-public-sources-20260920.json) |
| Mitochondrial references/indexes | Public references | Acquisition/sequence verification | [Record](../results/mitochondrial-public-resource-verification-20260920.json) |
| Mixture-screen runtime sources | Public software sources | Acquisition/version review | [Record](../results/mutserve-runtime-source-review-20260920.json) |
| GATK regional-ploidy behavior | Upstream source review | Pinned version semantics | [Record](../results/ploidy-boundary-public-source-review-20260920.json) |

## Reproduction guides and historical records

The root [reproduction guide](../REPRODUCIBILITY.md) separates the portable public
benchmark from the controlled-access historical case workflow. Historical path and
source maps remain linked below for auditors; their old host/version labels are not
prerequisites for understanding the project.

| Component | Data type | Control count / scope | Receipt or guide |
|---|---|---|---|
| Frozen small-variant comparison | Documentation / recorded evidence | See linked scope | [public-small-variant-comparison.json](../workflow/configs/public-small-variant-comparison.json) |
| changed-call read-review recipe | Documentation / recorded evidence | See linked scope | [changed-call-read-reproduction.md](../workflow/docs/changed-call-read-reproduction.md) |
| completed-shard handoff audit | Documentation / recorded evidence | See linked scope | [discovery-handoff-reproduction.md](../workflow/docs/discovery-handoff-reproduction.md) |
| native coding review | Documentation / recorded evidence | See linked scope | [NATIVE-CODING-REVIEW.md](../workflow/docs/NATIVE-CODING-REVIEW.md) |
| joint-haplotype review | Documentation / recorded evidence | See linked scope | [JOINT-HAPLOTYPE-REVIEW.md](../workflow/docs/JOINT-HAPLOTYPE-REVIEW.md) |
| clipped-read workflow | Documentation / recorded evidence | See linked scope | [CLIPPED-READ-REVIEW.md](../workflow/docs/CLIPPED-READ-REVIEW.md) |
| reproduction guide | Documentation / recorded evidence | See linked scope | [CODING-DEPTH-REPRODUCTION.md](../workflow/docs/CODING-DEPTH-REPRODUCTION.md) |
| overlay method | Documentation / recorded evidence | See linked scope | [ADAMTS1-NORMALIZED-OVERLAY.md](../workflow/docs/ADAMTS1-NORMALIZED-OVERLAY.md) |
| secondary contract | Documentation / recorded evidence | See linked scope | [NATIVE-SECONDARY-CONTRACT.md](../workflow/docs/NATIVE-SECONDARY-CONTRACT.md) |
| secondary reproduction guide | Documentation / recorded evidence | See linked scope | [NATIVE-SECONDARY-REPRODUCTION.md](../workflow/docs/NATIVE-SECONDARY-REPRODUCTION.md) |
| v26 reproduction guide | Documentation / recorded evidence | See linked scope | [V26-SOURCE-REPRODUCTION.md](../workflow/docs/V26-SOURCE-REPRODUCTION.md) |
| native layout map | Documentation / recorded evidence | See linked scope | [V26-NATIVE-RESTORATION.json](../workflow/docs/V26-NATIVE-RESTORATION.json) |
| nine-suite qualification map | Documentation / recorded evidence | See linked scope | [V26-QUALIFICATION-MAP.json](../workflow/docs/V26-QUALIFICATION-MAP.json) |
| source and qualification map | Documentation / recorded evidence | See linked scope | [V27-SOURCE-MAP.json](../workflow/docs/V27-SOURCE-MAP.json) |

Prepared case findings are in [reports](../reports/README.md); they are separate from
these public and synthetic validation records. [Change history](../CHANGELOG.md).
