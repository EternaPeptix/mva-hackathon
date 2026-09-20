# Saved-depth independent checker source derivatives

These are separate source derivatives, not the original executed sources. [Source provenance](source-provenance.json) records original and derivative hashes without exporting observed case outputs.

[audit_saved_gaps.py](audit_saved_gaps.py) independently reconstructs the declared saved-depth tables from the original GFF and pinned saved outputs. Its only change derives a receipt label from the computed interval list. The changed full checker has not been executed or qualified. The original mapper's 39-control receipt is an input gate; it is not qualification of this checker.

[check_report_tables.py](check_report_tables.py) reconciles compatible saved report tables and command metadata. Numeric comparisons are derived from the separately completed independent verification, not embedded observed counts. It accepts `--base`, `--original-depth`, and required `--output`; output creation is exclusive. Local comparison on the same saved inputs reproduced the original checker's numeric fields, input pins, and command metadata; the original checker was replayed with only its output destination redirected in memory. This is saved-output parity, not fictional qualification or new read-level validation. It does not programmatically assess scientific prose or clinical implications.

[replay_union_helper_examples.py](replay_union_helper_examples.py) reconstructs seven retained examples and executes only isolated `need` and `merged` definitions. [The receipt](union-helper-example-controls.json) explicitly preserves list/tuple type labels. The original inline runner was not retained, so this is not a byte-exact historical replay. No full-checker fixture suite is claimed.

Use standard Python 3 without `-O` or `PYTHONOPTIMIZE`. All three sources use only the standard library. The helper replay can run after relocation:

```text
python3 replay_union_helper_examples.py --source audit_saved_gaps.py --output new-helper-receipt.json
```

The main checker intentionally retains its Mac15,14/Darwin/arm64 host gate, historical root paths, fixed source/configuration/output identity hashes and public target geometry (84 genes, 2,578 transcripts, 2,112 windows, 393,880 reference bases). These are nonportable scope constraints, not measured outcome assertions. Restoring these source files alone does not restore the separately retained original GFF, target files, runtime, saved depth arrays, terminal receipts or reports. The report checker likewise needs its supplied saved-data directories and linked files. Do not overwrite any original source or receipt when restoring a derivative.

No alignment, read sequence, genotype, case-result receipt, observed depth table, private report, execution log or media is included here. Identity hashes only identify those separately retained artifacts; they do not reproduce their contents.
