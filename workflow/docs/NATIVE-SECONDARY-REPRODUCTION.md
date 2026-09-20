# Reproducing the native secondary inventory controls

The final producer is the exact qualified version 5, SHA256
`5e190155a01a9dc6144c123f79e6630a4110e8d062862c5b247115189d7e3786`.
Its 35 fictional controls passed on 512S2. The independent checker is SHA256
`ff438cd0c5aa14b066b83645ea09b56bcff2dc198dd8198191f3c48556413185`;
its 36 fictional controls passed separately. These are software and provenance
checks, not biological validation or a clinical reporting protocol.

The qualified `NATIVE-SECONDARY-CONTRACT.md` is preserved byte for byte as
historical method documentation, including its preparation wording. This guide
describes the final layout and reproduction prerequisites without rewriting that
qualified document. The frozen ClinVar filename, snapshot label and literal
header date are distinct fields, and an undeclared optional CLNALLE remains null.

## Source layout

- `workflow/scripts/review_native_secondary.py` is the producer.
- Its two immutable dependencies are adjacent under `workflow/scripts/vendor/`:
  `screen_secondary_findings.py` and `native_coding_records.py`.
- `workflow/tests/test_native_secondary.py` generates fictional inputs only.
- `workflow/scripts/qualify_native_secondary.py` snapshots the entire workflow
  into a new version directory and executes only that producer test module.
- `workflow/source-snapshots/native-secondary-independent-20260920/source/`
  keeps `verify-secondary.py`, `test_verify_secondary.py` and
  `qualify-verifier.py` together. The checker imports neither producer nor vendor
  code. Its test module imports a qualified producer fixture solely to generate
  fictional integration inputs and intentional corruptions.
- Historical producer and independent `controls.json` receipts remain unchanged.
  They are evidence of their original runs, not relocation-aware qualification
  receipts for a new checkout.

Python must have pysam and the pinned bioenvironment must provide bcftools. No
automatic package or database installation is part of these commands. Native
controls require zero skips; a standard-library-only local run explicitly skips
the integration controls and does not substitute for native qualification.

## Exact native control commands

From the integrated repository root on the authorized Mac Studio, with a new
version name that does not already exist:

```sh
env PATH="/Users/jeweled/Documents/mva-hackathon/tools/bio/bin:$PATH" PYTHONDONTWRITEBYTECODE=1 \
  /Users/jeweled/Documents/mva-hackathon/tools/bio/bin/python -B \
  workflow/scripts/qualify_native_secondary.py \
  --version secondary-reproduction-producer-v1 --require-no-skips
```

That command creates `secondary-reproduction-producer-v1/controls.json` beside
the repository's workflow directory and a matching adjacent workflow snapshot.
It does not load an approved case configuration or open original patient sources.
Preserve the printed receipt hash and never replace an existing version directory.

The following exact independent-control command uses the original private,
immutable producer qualification, which must still exist with its entire adjacent
`workflow/` snapshot. It generates only fictional inputs and audit outputs in a separate new archive,
so it does not mutate the sealed completed method archive:

```sh
env PATH="/Users/jeweled/Documents/mva-hackathon/tools/bio/bin:$PATH" PYTHONDONTWRITEBYTECODE=1 \
  /Users/jeweled/Documents/mva-hackathon/tools/bio/bin/python -B \
  workflow/source-snapshots/native-secondary-independent-20260920/source/qualify-verifier.py \
  --output /Users/jeweled/Documents/mva-hackathon/methods/native-secondary-control-reproduction-20260920-v1 \
  --producer-qualification /Users/jeweled/Documents/mva-hackathon/methods/native-secondary-review-20260920/qualification-v5-512s2/controls.json \
  --expected-producer-qualification-sha256 9c88971c0f73a52287687e067e07dcc08d798f34c3c94cd983208c3562ed805a
```

To qualify a newly relocated producer snapshot instead, replace the last two
arguments with the fresh producer receipt path and its exact SHA256 from the first
command. The entire adjacent `workflow/` tree is required: the independent
fixture and checker deliberately verify that source inventory. A copied
historical receipt alone cannot satisfy that requirement.

## Private case prerequisites

The approved version 2 configuration contains the sample identifier, absolute
input paths, source sizes/SHA256 values and explicit resource metadata. It contains
no actual matching-allele evidence. Preserve it as an execution record; its
`ready: true` is not blanket authorization for a new run. A new case run needs
separate authorization, a fresh output directory, fresh source/receipt checks and
an appropriately qualified guarded launcher.

Executing the actual inventory additionally requires the original private
46-source callset, complete 22-autosome handoff/verification receipts, separately
completed fixed-X/Y status/review/verification and exact reference/tool/public
files. The source-only repository and historical control receipts do not contain
those inputs and cannot reconstruct them. No variant caller, alignment extraction
or newer database should be substituted to make missing prerequisites pass.

The producer CLI takes `--config`, `--expected-config-sha256` and `--output`.
The independent checker takes the exact config, completed producer status and
producer qualification hashes, the producer directory and a separate new audit
directory. The private completed artifacts document the specific guarded
invocations. Their resource limits were 90 minutes, 16 GiB group RSS, 20 GiB output
and at least 40 GiB free disk, sampled every two seconds, with nice +10 and signals
restricted to each newly created process group.

Preserve missing versus zero FORMAT values, every original/normalized lineage,
caller separation and the stated variant-only scope. Successful reproduction does
not establish complete panel callability, reference genotypes at absent records,
parental phase, reportability or a clinical negative.
