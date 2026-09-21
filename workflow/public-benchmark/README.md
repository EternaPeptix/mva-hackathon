# Public HG002 chromosome 20 replay

From a fresh clone, with Bash, the standard `ps` utility (Linux `procps`) and Python 3.9 or newer already installed:

```bash
bash workflow/public-benchmark/run.sh ../hg002-public-replay
```

Choose a **new output directory outside the clone**, under an existing parent directory. The command downloads public data and isolated dependencies, prepares the reference, runs both callers over the complete public chromosome, compares both callsets with GIAB truth, and writes a self-contained `summary.json` plus logs and exact input/runtime/source receipts. It never reads the private case. Existing output directories are refused rather than overwritten. A failed run retains its evidence; there is no automatic retry.

The public input is the **prealigned**, nominal-35× HG002 NovaSeq PCR-free chr20 BAM. This replays a calling/comparison component; it does not reproduce FASTQ alignment, whole-genome accuracy, measured 35× coverage, clinical sensitivity, parental phase or causal ranking. The full 195-sequence reference dictionary and unchanged contig sharding are retained; mapped reads occur only on chr20. All official chr20 GIAB v4.2.1 confident regions are evaluated. No tiny demonstration substitutes for the full default run.

## What runs

1. Install the checksum-pinned Micromamba bootstrap into the output tree, then install the platform-specific **explicit package lock** without changing the user's shell, default environment or package configuration.
2. Fetch RTG Tools 3.13; verify its exact wrapper and jar. Relocate only its Java path in a new `rtg.cfg`, preserving memory/thread/telemetry settings. Run the existing twelve comparator semantic controls. A failing control stops the replay.
3. Download and verify the six public BAM/index/reference/truth/index/BED files, plus the GATK 4.6.1.0 distribution. The exact extracted GATK jar must match before it executes. The GATK ZIP has no retained publisher digest; the executable jar has the historical local SHA-256 pin.
4. Decompress the unchanged complete reference, verify its SHA-256, sample/read-group/dictionary and mapped-record scope, and run the existing reference preparation including the chr20 reference round trip.
5. Run the unchanged `call_genome.py` with GATK 4.6.1.0 and bcftools 1.24. The default is two workers; `--workers 1` through `--workers 4` are supported. Each JVM has a 6-GiB heap and two advertised CPUs. Historical calling settings, ploidy, whole-contig geometry and split normalization are unchanged.
6. Run the unchanged genotype-aware RTG comparison and record-accounting audit. Preserve full chromosome context; the BED is supplied as evaluation regions. No threshold is tuned to maximize F-measure.

On Linux, the launcher sets `JDK_JAVA_OPTIONS=-Dsun.zip.disableMemoryMapping=true` to avoid the observed JDK class-loading failure on an emulated/virtiofs mount. This changes jar I/O rather than caller/comparator algorithms. Inherited `JAVA_TOOL_OPTIONS`, `_JAVA_OPTIONS`, `JDK_JAVA_OPTIONS`, `JAVA_OPTS`, `JAVA_OPTIONS` and `JVM_OPTS` are cleared on every platform before the disclosed Linux option is set. The generated Java path is shell-quoted, including paths containing shell metacharacters.

The new adapter changes downloads, runtime location and platform verification. It does not edit the historical scripts or their receipts. It replaces the historical Mac-binary-only `runtime()` locator in memory; all caller arguments, comparison arguments, allele checks, reference checks and control expectations remain unchanged. `profile.json` pins the reused sources. New runtime/source receipts distinguish this replay from the historical machine-bound execution.

## Outputs and interpretation

- `status.json`: guard status, platform, process-group identity, elapsed time and sampled limits.
- `summary.json`: final raw TP-baseline/TP-call/FP/FN counts, baseline recall and call precision for each caller, plus exact provenance links. It exists with a complete benchmark result only after all stages pass.
- `comparison/`: complete native RTG outputs and per-class summaries; the same scored truth universe is required for both callers.
- `audit/`: original-record accounting and bounded examples, not a new matching algorithm.
- `pilot/calls/`: raw and normalized caller outputs and original contig/shard completion receipts.
- `controls/`, `logs/`, `input-receipt.json`, `portable-runtime.json`, `source-provenance.json`: retained control, execution and identity evidence.

The frozen primary representation uses `bcftools norm -m -any`. The separately reported post-observation unsplit-multiallelic sensitivity is **not** silently substituted. Record-level recall and precision can have different truth/call representation denominators. The historical run scored 82,589 of 85,951 prepared truth records; the new run reports its own observed counts instead of copying those numbers.

## Platforms and resources

As of 21 September 2026, the new full macOS arm64 replay is running after verified public downloads, reference preparation and twelve passing RTG controls. Linux x86-64 under QEMU completed the pinned environment, all twelve controls, and actual GATK/bcftools calling plus both RTG comparisons on a predefined public 100-kb interval. That is an operational smoke, not a full Linux chromosome benchmark. Linux aarch64 has a resolved lock only. The final launcher adds shell-safe Java-path quoting and inherited-JVM-option sanitation after the Mac run started; its launched source remains recorded separately. The [current platform record](../../results/public-benchmark-replay-20260921.json) distinguishes these stages and retained compatibility failures.

Explicit SHA-256 package locks are supplied for `linux-64` (x86-64, glibc 2.28 or newer), `linux-aarch64` and `osx-arm64`. They pin Python 3.12.14, pysam 0.24.1, samtools/bcftools 1.24, OpenJDK 21.0.10 and all resolved dependencies. They were solved separately; a solved environment is not evidence that the workflow ran on that platform. Current execution status is recorded separately in the review results; do not infer Linux verification from a Mac run.

Public data transfers total **2,149,884,390 bytes**. GATK/RTG, bootstrap and Conda packages add about another 1 GB; compressed archives and the installed environment are retained. Prefer at least **24 GiB RAM** for the default two-worker profile; `--workers 1` reduces concurrent JVM heaps to 6 GiB. The small Linux smoke does not establish that the full replay fits its 6-GiB container. Plan roughly **15–25 GiB** of working storage; this is an estimate, not a measured maximum. At least **40 GiB free** is required. The default guard stops only its owned process group after 6 hours, over 32 GiB sampled aggregate RSS, over 30 GiB output, or below 40 GiB free. Checks occur approximately every five seconds and are not instantaneous kernel caps. The job runs at reduced scheduling priority. It does not start Docker, VMs, cloud jobs or other services.

The historical full caller run took **6,969.6 seconds (about 116 minutes)** on its Mac Studio; reference preparation, downloads and comparison are additional. This is not a runtime guarantee for another computer or architecture. Dependencies and public resources must remain available over HTTPS. Existing primary GIAB truth hashes are locally recorded integrity pins, not publisher-published checksums; the profile states that distinction.

For a deliberately limited installation/comparator check:

```bash
bash workflow/public-benchmark/run.sh ../hg002-controls-only --controls-only
```

This installs the pinned runtime and runs twelve fictional RTG semantic controls. It does not download the BAM or truth, call public variants, or establish benchmark performance. It must not be reported as full replay verification.

Dependency sources: [Micromamba manual installation](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html), [GATK 4.6.1.0](https://github.com/broadinstitute/gatk/releases/tag/4.6.1.0), [RTG Tools 3.13](https://github.com/RealTimeGenomics/rtg-tools/releases/tag/3.13). Exact public data URLs and hashes are in `profile.json`.
