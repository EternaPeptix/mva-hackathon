# Saved public SV residual sources

This snapshot preserves the exact code, approved input configuration, historical launch plans and fictional-control receipts used for the completed saved-VCF public HG002 chr20 residual review. It adds no caller, rematching, read analysis or annotation policy. Actual record ledgers, summaries, audit results, launch logs, raw inputs, media and patient data are omitted.

[SOURCE-MAP.json](SOURCE-MAP.json) records all 12 historical files, their native restoration paths, six qualification source bindings, two exact dependencies already in v27 and 17 external prerequisite pins referenced by the plans. The base is private source commit `400f7991116f087c4c64b4a9e6ab74467b877eb0`. The map documents source completeness and dependencies; it does not make a fresh clone able to reproduce the benchmark result.

## Files and distinct control scopes

| Component | Preserved files | Evidence scope |
| --- | --- | --- |
| Producer | [Analyzer](producer-preparation/qualification-v1/source/analyze_sv_residuals.py), [tests](producer-preparation/qualification-v1/source/test_analyze_sv_residuals.py), [controls](producer-preparation/qualification-v1/controls.json) | 29 native fictional parsing, identity, annotation, accounting and producer-integration controls; two exact source pins. |
| Execution | [Approved configuration](root-execution-v1/files/residual-approved-v1.json), [producer plan](root-execution-v1/files/execution-plan.json), [parent guard](root-execution-v1/files/run_saved_review_guard.py) | Exact historical selectors and executed guard; these three files are outside both fictional suites' source membership. |
| Independent readback | [Checker](independent-execution-review/qualification-v2/source/verify_saved_residuals.py), [separate parser](independent-execution-review/qualification-v2/source/independent_vcf_reader_v2.py), [six parser tests](independent-execution-review/qualification-v2/source/test_reader_v2.py), [three helper tests](independent-execution-review/qualification-v2/source/test_independent_readback.py), [controls](independent-execution-review/qualification-v2/controls.json) | Nine native fictional parser/identity/missingness/aggregate-helper checks, with four exact source pins. The full checker was imported, not exercised end-to-end on a fictional fixture. This is not a nine-test guard suite. |
| Independent execution | [Independent plan](independent-execution-review/execution-v1/plan.json) | Exact historical selectors, including expected producer-result paths and hashes; no actual result payloads. |

Both receipts record Python 3.12.14 on Mac Studio 512S2. The producer represents no skips as `[]`; the independent receipt represents no skips as integer `0`. Preserve those schemas. These are separate 29- and nine-control scopes, not a combined accuracy or biological-validation score. No controls were rerun to prepare this export.

Keep each source directory intact: the producer test loads its adjacent analyzer; the checker and its two tests import the adjacent `independent_vcf_reader_v2.py`. To rerun only the fictional suites using a deliberately selected Python runtime:

```text
python3 -B -m unittest discover -s workflow/source-snapshots/public-sv-residual-review-20260920/producer-preparation/qualification-v1/source -p test_analyze_sv_residuals.py -v
python3 -B -m unittest discover -s workflow/source-snapshots/public-sv-residual-review-20260920/independent-execution-review/qualification-v2/source -p 'test_*.py' -v
```

The scripts use Python's standard library. The historical native runtime executable is separately pinned in both receipts and plans; it is not bundled here. A new environment's test pass would be new qualification evidence, not an alteration of the original receipts.

## Reused dependencies

The guard imports the existing [guard helper](../secondary-depth-gap-map-20260920/qualification-v1/workflow/scripts/vendor/guard_helpers.py), SHA256 `e30e5ad4819d488e2509feffff1027bab6102e473bb44ed33bcf4dbfc734ffb6`. Its historical native path is `/Users/jeweled/Documents/mva-hackathon/methods/secondary-depth-gap-map-20260920/qualification-v1/workflow/scripts/vendor/guard_helpers.py`. It is reused by hash and restored path, not duplicated in this snapshot. Run the new parent guard, not the helper's unrelated retained command-line entry point.

The original [public comparison configuration](../../configs/public-sv-pilot-comparison.json), SHA256 `86b2e380e5dd390b05ee7deb3393171976ebe84f828bb556e3e4f5ef8f667e80`, is also already present in v27. Its historical path is under the native `repo/workflow/configs` directory. All other plan-pinned source/configuration/receipt files in this snapshot map to their exact native `methods/public-sv-residual-review-20260920` locations in SOURCE-MAP.

## Historical plans and deliberate reproduction

The two plans are **nonportable historical selectors and hash metadata**, preserved byte-for-byte. Their `ready:true` values record the original approvals; they do not authorize another run. The producer config retains documentary public benchmark-count gates and input pins. Its inclusion supplies neither the twelve saved inputs nor an actual ledger/result.

The producer needs the pinned refined baseline VCFs for the two callers, the original no-star chr20 full-context comparison VCF, the official benchmark BED, prior configuration/audit/count/terminal receipts and the exact runtime. The two refined baseline VCFs and original chr20 context are saved derivatives of the earlier public experiment, not resources created by these source files. The independent checker additionally needs the completed producer outputs and terminal launch record. Raw VCF/BED payloads, actual producer outputs and process-snapshot-bearing launch status are intentionally absent. The plans preserve only their path/byte/hash selectors.

For an authorized reproduction, restore the qualified source layout and exact available inputs, preserve all existing outputs, and create **new separately reviewed plans**. Choose new result and launch directories, update command arguments and the matching plan fields together, and freshly pin every changed plan, path, source, runtime and control receipt. If input locations change, create a new pinned configuration and update all dependent references rather than editing this archival snapshot. Only after a new producer run completes can a new independent plan bind its fresh output/status hashes and corresponding `--run` and `--launch` paths.

The historical command arrays in the plans describe the exact producer and checker CLIs. The intended parent invocation is:

```text
<qualified-python> <guard-path>/run_saved_review_guard.py --plan <new-reviewed-plan.json> --expected-plan-sha256 <new-plan-sha256> --expected-guard-sha256 b8b582639e425599c588d563fc2dc7d87909657660e69fcb4882ccbf2c758e49
```

The guard requires 512S2, native arm64 and Mac15,14; a new exclusive output/launch pair; the pinned helper and command inputs; and fixed sampled limits of 300 seconds, 2 GiB process-group RSS, 64 MiB output and at least 40 GiB free. It applies nice+10 and starts a separate process group. Limits are sampled, not instantaneous kernel-enforced caps. Existing output directories are refused. There is no automatic retry, and unrelated work must remain untouched.

The analyzer's output contract is `ledger.jsonl`, `summary.json`, `headers.json` and terminal `status.json`. The checker reconstructs the ledger and tables through its separate parser and emits verification/summary/status files. None of those actual payloads is part of this source-only snapshot. Literal annotation absence, genotype order/phase encoding and descriptive reference-span containment retain the source's stated limitations; neither source availability nor passing fictional controls establishes clinical sensitivity or a negative patient finding.
