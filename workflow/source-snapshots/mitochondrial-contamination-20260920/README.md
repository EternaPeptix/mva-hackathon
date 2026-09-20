# Exact mitochondrial caller/importer control sources

These are the executed fictional-control and native Java probe sources from
Mac Studio 512S2. They contain no patient reads or result tables. The two Python
control scripts record the study root `/Users/jeweled/Documents/mva-hackathon`;
review and change that root deliberately when reproducing elsewhere. Preserve
existing output directories rather than rerunning into them.

1. Acquire the official archive recorded in
   `../../configs/mitochondrial-contamination-runtime.json` and check its SHA-256.
   Extract `mutserve.jar` and `rCRS.fasta` into
   `ROOT/tools/mutserve-2.0.0-rc12`, and `haplocheck.jar` into
   `ROOT/tools/haplocheck-1.3.3`. Check the separate file identities. Native Java
   21 and Python/pysam use the recorded `ROOT/tools/bio` environment.
2. Place the reviewed repository under `ROOT/repo`, and copy these five source
   files into `ROOT/methods/mutserve-qualification-20260920`.
3. Run `run-contract-controls.py` with the pinned Python. It creates a new
   `ROOT/results/mutserve-contract-controls-20260920` directory, synthesizes known
   pileups, executes the distributed caller and probes both distributed JARs.
4. Run `verify-contract-controls.py`. It checks the observed failure modes and
   verifies the selected caller records through the actual haplocheck importer.
   Its qualification receipt is a required input to the case wrapper.
5. `MutserveAlleleRankProbe.java` is an additional fictional rank-switch control.
   Compile against the pinned mutserve JAR, then run its class with a new JSON
   output path. A 2% G allele on both strands loses eligibility when an additional
   3% forward-only A allele changes the minor-allele ranking. This is a software
   mechanism test, not an empirical error model.

The case wrapper is `../../scripts/run_mitochondrial_contamination.py`; its
`--help` lists required paths. It requires input/reference/sample identities from
the reviewed configuration, a verified native-control receipt and
`HaplocheckParserProbe.java`. Its new output directory and all parent directories
must have no periods because the pinned caller derives sidecar prefixes at the
first period. It preserves native outputs and validates selected SNP semantics
before generating hotspot-inclusive and hotspot-excluded reports.

The selection helper is `../../scripts/prepare_haplocheck_input.py`. No GT/AF
recoding is performed. An empty qualified input is reported as uninformative.
NO/ND must not become zero contamination, a nuclear contamination estimate or a
clinical negative. Error-free synthetic reads do not validate a 1% detection limit.
