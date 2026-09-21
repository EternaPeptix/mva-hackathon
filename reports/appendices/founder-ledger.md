> Package copy: local links/navigation are adapted for this finite review package. Labels identify underlying files that are not included. Scientific text is retained; approved status/source metadata changes are recorded separately; original and packaged hashes are recorded in PROVENANCE.

> **Retained dated method detail.** This appendix preserves the original analysis or proposal at its stated date. The current [Track 1 report](../jeweled_track1_report_submission.md), [Track 2 report](../jeweled_track2_report_submission.md) and [methods form](../jeweled_methods_form_answers.md) govern interpretation; see the [revision and evidence scope](../PROVENANCE.md#evidence-scope). Historical delivery or ranking statements are not current status.

# Descriptive founder-ledger analysis

`workflow/scripts/analyze_founder_ledger.py` checks an adjudicated first-mitosis
ledger against a separately frozen pretreatment roster. It reports the linked
readouts in the prospective pilot design. It does not track images, adjudicate
biology, verify true pretreatment freezing, prove the roster is complete, or
determine whether a drug should advance. No laboratory data are included in the
fictional controls.

## Four pinned inputs

All inputs are UTF-8 JSON. Each file requires an externally supplied lowercase
SHA256 digest. The design also embeds the founder-roster digest. Verify and freeze
the roster independently before treatment; an internally consistent file created
after outcomes were known would still pass mechanical validation. All input
digests are checked before analysis. Duplicate JSON keys, NaN/Infinity/overflow,
unknown fields and implicit conversions from numbers to booleans are rejected.

The three tables have exactly `schema_version: 1` and one array: `wells`,
`founders`, or `outcomes`. Every listed row field is required; use explicit
`null` where its value is unavailable or inapplicable. Empty strings are invalid
identifiers or reasons. The design is a single object with the fields below.

| Design field | Meaning |
|---|---|
| `schema_version` | Integer 1 |
| `data_kind` | `synthetic_control` or `observed`; fictional controls stay labeled |
| `acute_end_hours`, `secondary_end_hours` | Finite positive acute duration and secondary boundary at least as late |
| `time_unit` | Exactly `hours_after_time_zero` |
| `roster_frozen_at_utc`, `scoring_frozen_at_utc` | Explicit UTC timestamps no later than any founder's time zero |
| `roster_sha256` | Exact digest of the separate founders input |
| `adjudication_frozen` | Explicit `true`; identifies the supplied records as the frozen primary adjudication |
| `scoring_rule_id` | Frozen scoring-rule identity |
| `error_criteria` | Nonempty object mapping each prespecified error label to its definition; this content is pinned with the design |
| `source_acquisition_id` | Identity of the acquisition or explicitly fictional source |

## Planned wells and actual hierarchy

Each well has exactly: `experiment_id`, `well_id`, `block_id`,
`masked_condition_id`, `assay_id`, `assay_scope`, `model_id`, `material_origin`,
`source_background_id`, `editing_family_id`, `clone_id`, `preparation_id`,
`experimental_day`, `expected_founders`, `acquisition_status`, `failure_reason`.

The well key is `(experiment_id, well_id)`. `assay_scope` is `spontaneous`,
`challenge`, or `assay_control`. `material_origin` is `engineered`,
`patient_derived`, or `other`. Engineered material needs its actual editing
family; other origins must use a null family. A clone ID maps consistently to
origin, source background, editing family and model throughout these inputs.
Preparations may cross editing families; no false nesting is imposed. IDs must
refer to real design units, not invented replicates. Keep the unmasked assignment
key outside these scoring records.

`expected_founders` is the independently declared starting count: a nonnegative
integer when known, otherwise null. The roster must match it exactly. There are
three acquisition statuses:

- `acquired`: known count and null failure reason. A genuinely observed empty
  well has count zero and empty roster/outcome sets; its ratios are undefined.
- `failed`: a nonempty failure reason. If a frozen starting census exists, retain
  its known count, full roster and every outcome. An interrupted acquisition may
  retain earlier A/U/S and classify remaining founders as L with time/reason.
  If no starting census is available, use a null count and no primary rows; that
  planned well is explicitly unanalysed.
- `unrun`: null starting count, nonempty reason and no primary rows. Planned
  cells are not an observed starting census.

Missing acquisitions and unknown starting counts are never analyzed as zero
founders. No failed well with a supplied roster is silently discarded.

## Pretreatment founders

Each roster row has exactly: `experiment_id`, `well_id`, `founder_id`, `field_id`,
`source_cell_id`, `raw_image_id`, `parent_id`, `time_zero_utc`,
`baseline_mitotic_state`.

The primary key is `(experiment_id, well_id, founder_id)`. `parent_id` must be
null: descendants are not primary founders. `raw_image_id` identifies a specific
source image/acquisition object across the submitted dataset; `field_id` is a
field within that object and `source_cell_id` is a cell/track within that field.
The tuple `(raw_image_id, field_id, source_cell_id)` must be unique even across
well or experiment labels. Reusing a field-local cell number in a different field
is allowed. Relabeling one source cell under another well is rejected. These
opaque identifiers preserve declared provenance, not proof of image identity.

`time_zero_utc` must have an explicit zero UTC offset. Baseline state is
`interphase`, `mitotic`, or `unknown`. An already-mitotic founder has a
left-censored entry time; do not invent an entry at zero.

## One primary adjudicated outcome per founder

Each outcome has exactly: `experiment_id`, `well_id`, `founder_id`, `parent_id`,
`status`, `classification_hours`, `observation_end_hours`,
`observed_through_window`, `entry_hours`, `entry_unknown_reason`,
`anaphase_hours`, `exit_without_anaphase_hours`, `death_hours`, `loss_hours`,
`death_stage`, `reason`, `visible_error`, `error_types`, `secondary`.

`parent_id` is null. The outcome identity set must equal the frozen roster set;
missing, extra, duplicate or descendant primary rows fail. `classification_hours`
is the biological first-event or administrative endpoint time, **not** the time
when a scorer later reviewed the image. All primary times are finite numbers in
the acute window. Entry precedes classification; observation continues through
classification. `observation_end_hours` describes the retained primary follow-up;
separately recorded later fate belongs in the secondary ledger.

`observed_through_window` is explicit true/false. True requires primary follow-up
through the acute endpoint. B/O require true, classification at the endpoint and
no previous first-mitosis outcome. B has no entry and requires an adjudicated
interphase baseline. An unknown baseline cannot establish B; review the source
and adjudicate it rather than automatically relabeling it L, which specifically
means tracking loss. O
asserts that the first mitosis remains active at that endpoint; it is not a
permanent-arrest claim.

- A/U require their first anaphase at classification, and no competing primary
  slippage/death/loss field. A needs a strict boolean `visible_error` and unique
  known error labels: true iff the label list is nonempty. U needs null error
  status, empty definite-error labels, and an unscorable reason.
- S requires exit without anaphase at classification and an adjudication reason.
- D requires death at classification, primary follow-up ending there and an
  explicit stage (`interphase`, `mitotic`, or `unknown`). Interphase death cannot
  contradict entry or a mitotic baseline. Do not fabricate a stage.
- L requires loss at classification, primary follow-up ending there and a reason.
- Only A can have scored primary error labels; all other error indicators are
  null and their definite-error lists empty. U, S and L require reasons.

For an interphase baseline, A/U/S/O and mitotic D require an observed entry time.
For mitotic or unknown baseline, a missing entry needs an explicit reason.
Postbaseline entry is rejected for an already-mitotic baseline. D/L cannot be
combined with supplied earlier anaphase/slippage evidence. Missing evidence that
was never entered cannot be detected by this validator. Ambiguous event ordering
must be adjudicated under the frozen rules before this input is produced; the
tool does not invent a tie-break from coarse image times.

## Optional secondary observations

`secondary` is null when unavailable, or a list (possibly empty) of separately
supplied records. Each has `event_id`, `subject_id`, `parent_id`, `kind`,
`time_hours`, `observed`, `note`. Kinds are `death`, `tracking_loss`, `division`,
`cytokinesis_failure`, `other`. `observed` is a strict boolean; missing is not
false. The time must follow primary classification and not exceed the secondary
boundary. These are supplied event/absence observations, not automatic tracking
or a complete daughter census. A founder subject has null parent; a descendant
has a nonempty parent ID. The tool does not reconstruct a full lineage graph.

Later death/loss never overwrites fixed A/U/S. Secondary follow-up cannot silently
resume a primary D/L record. A daughter division never increases N0, A or E.
Secondary summaries count supplied records only, retain unavailable founder
ledgers, and use null for event kinds with no supplied observations. They are not
incidence estimates or evidence of normal daughter function.

## Descriptive outputs and invocation

The output directory must be new and its parent must already exist. All inputs
are validated before the directory is created. `summary.json` retains every
well, roster/outcome row, frozen design, input hashes and analyzer source hash.
`status.json` pins the completed summary. No existing output is overwritten.

```sh
python3 workflow/scripts/analyze_founder_ledger.py \
  --design design.json --design-sha256 DESIGN_SHA256 \
  --wells wells.json --wells-sha256 WELLS_SHA256 \
  --founders founders.json --founders-sha256 FROZEN_ROSTER_SHA256 \
  --outcomes outcomes.json --outcomes-sha256 OUTCOMES_SHA256 \
  --output-dir new-ledger-result
```

Per-well counts report all A/U/B/S/D/O/L and require their sum to equal roster N0.
E counts each A founder once even with several error types. E/A is null at A=0.
The observed yield `(A−E)/N0`, A/N0 and all status proportions are null at N0=0.
Yield without a scored visible error is not a normal-karyotype result. The U
bounds are E/(A+U) and (E+U)/(A+U), null at A+U=0; L is excluded and remains
separate. Primary classification completeness is (N0−L)/N0, not complete later
lineage surveillance. Individual primary timing records remain available.

Block-condition descriptive sums retain **every** group field: experiment,
block, masked condition, assay ID/scope, model, origin, source background, family,
clone, preparation and day. Thus only technical wells with the same declared
stratum are combined. Their denominators are summed, not averaged cell rates.
Planned/analysed/acquired/failed/unrun well counts remain visible. There are no
cross-origin or cross-assay pooled estimates, confidence intervals, p-values,
power calculations, automatic advancement or efficacy/clinical conclusions.

Run fictional controls with `python3 -B workflow/tests/test_founder_ledger.py`.
The tests include the protocol's 80/16 versus 20/2 counterexample, unavailable and
empty censuses, interrupted acquisition, strict errors/times/hashes, retained
first outcomes, crossed preparations, unequal wells, and source-cell identity.
