# Native joint-sequence review with a separate normalized overlay

This method preserves full native caller records and evaluates explicit constructed sequences in bounded reference windows. It is a technical reconciliation of alternative representations, not a genotype caller, parental phasing procedure, clinical validation or measure of causal probability. It builds on the unchanged qualified `review_joint_haplotypes.py`, `joint_haplotype_helpers.py` and `native_coding_records.py` workflow.

The private configuration contains exact target/sample/input-hash selection. Original alignments, native variant files, source manifests and completed case evidence are separate private prerequisites. Source publication does not include individual reads, read names, case result payloads, clinical reports or media. A reproduction must obtain authorized inputs that match the configuration; a missing input is not a reference or negative finding.

## Preparation and qualification

Place the exact source snapshot files in `ROOT/methods/adamts1-joint-review-20260920` and the unchanged qualified native sources in their recorded qualification directory. Preserve the original configuration, overlay specification and source hashes. The `reference-preflight.py` command opens variant/reference inputs only. It replays original caller streams, verifies exact full source counts and rendered payload hashes, preserves native ordinals/full ALT/GT/FORMAT/FILTER/QUAL, and validates the declared normalized allele against its original context. It explicitly does not hash or open the alignment or its index.

Run `test_normalized_overlay.py --output NEW_CONTROLS_JSON` with `MVA_QUALIFIED_SCRIPTS` pointing to the unchanged qualified native script directory. This fictional suite covers literal equality of long raw REF and canonical representations, disjoint versus overlapping substitutions, preserved `GT1/2` spanning-star ALT2 and six PL values, equivalent class collapse, repeated observations versus fragments, mate conflicts, missingness, exact source and output refusal, and a completed-input fixture with the real native evidence schema. Terminal completion belongs to native summary/status files; evidence is independently pinned and does not require an invented status field.

Run the reference-only preparation with `--root ROOT --config CONFIG --spec SPEC --output NEW_SELECTION_JSON`. Keep every prior preparation or failed run and use new output paths. Its source/input pins must match the subsequently qualified overlay. If source bytes change, retain the old receipt and create a newly bound preparation rather than rewriting it.

## Native observation review

The native workflow uses full declared cores with flanks 3/10/25, mapping-quality sensitivities 0/20/60 and raw base-quality sensitivities 20/30. It excludes flags 0xF04 and unknown mapping quality 255, applies no BAQ, and requires each accepted original observation to span the complete window with valid CIGAR boundaries. Partial mates are not stitched. It verifies original decoded extraction, retains unmatched sequences and mate conflicts, and aggregates fragments by original read group and query name. Caller hypotheses remain separate during edit enumeration; overlapping raw REF edits are explicitly rejected. Every concrete ALT remains a labelled hypothesis regardless of native genotype participation. Nonconcrete stars retain their full native context but do not become invented concrete sequence edits.

The configuration's state and sequence-size guards are engineering limits. They do not establish biological sensitivity. Use an isolated low-priority process with finite resource limits and a new result directory. Preserve unrelated workloads and original caller state. Do not automatically rerun a failed case.

Before accepting an overlay input, independently verify original source ordinals/payloads, reference windows, extraction, sequence classes/lineages, overlap rejections, annotation bindings and all observation/fragment quality arms. Retain the complete verifier receipt and fresh hashes for native summary, status, evidence and regional context.

## Separate normalized representation overlay

`normalized-overlay.py` consumes the completed independently verified native case, source-bound fictional controls and reference-only selection. Its command requires `--root`, `--native-run`, `--output`, `--selection`, `--expected-selection-sha256`, `--controls`, `--expected-controls-sha256`, `--native-verification` and `--expected-native-verification-sha256`. All outputs must be new and separate from the native case. The overlay performs no read extraction, alignment, genotype calling or annotation.

The overlay proves full-window raw/canonical deletion equality before composing the explicitly declared disjoint substitution. It preserves rejected overlapping raw and canonical combinations. Derived lineages must match already enumerated native sequence classes; otherwise this scoped method refuses rather than silently extending the panel. Original native records, full star/ALT/GT/PL context and observer evidence remain unchanged. A sequence class receives one observation/fragment count regardless of how many native or derived lineages describe it.

An independent overlay audit must reconstruct the sequence algebra and regroup the pinned original observations without importing the overlay implementation. Match the entire output class/lineage and quality-arm schema, all source/receipt bindings and zero-support classes. Keep source/control qualification distinct from actual-data reconciliation.

## Interpretation and preservation

Annotations describe the constructed combined local sequence. Net-length arithmetic, an in-frame annotation or local read support does not prove protein function, rescue, locus origin, diploid genotype, parental phase, clinical reportability or causality. A zero accepted count is bounded by the exact observation gates and hypothesis scope. Nested windows and quality sensitivities reuse fragments; do not sum them or infer relative allele fractions.

Preserve original native and derived representations together. Keep prior source/preparation versions, completed and failed statuses, command logs and independent receipts. Final protection should use a canonical unique list of absolute paths with bytes, modification time and SHA256, with a separate terminal status linked to the seal hash. Scientific evidence remains separate from source-only publication.
