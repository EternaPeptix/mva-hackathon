#!/usr/bin/env python3
"""Validate a frozen founder ledger and report descriptive, denominator-safe counts."""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import sys


STATUSES = ('A', 'U', 'B', 'S', 'D', 'O', 'L')
WELL_FIELDS = ('experiment_id', 'well_id', 'block_id', 'masked_condition_id', 'assay_id',
    'assay_scope', 'model_id', 'material_origin', 'source_background_id', 'editing_family_id',
    'clone_id', 'preparation_id', 'experimental_day', 'expected_founders',
    'acquisition_status', 'failure_reason')
GROUP_FIELDS = ('experiment_id', 'block_id', 'masked_condition_id', 'assay_id', 'assay_scope',
    'model_id', 'material_origin', 'source_background_id', 'editing_family_id', 'clone_id',
    'preparation_id', 'experimental_day')
FOUNDER_FIELDS = ('experiment_id', 'well_id', 'founder_id', 'field_id', 'source_cell_id',
    'raw_image_id', 'parent_id', 'time_zero_utc', 'baseline_mitotic_state')
OUTCOME_FIELDS = ('experiment_id', 'well_id', 'founder_id', 'parent_id', 'status',
    'classification_hours', 'observation_end_hours', 'observed_through_window', 'entry_hours', 'entry_unknown_reason',
    'anaphase_hours', 'exit_without_anaphase_hours', 'death_hours', 'loss_hours',
    'death_stage', 'reason', 'visible_error', 'error_types', 'secondary')
TIME_FIELDS = ('classification_hours', 'observation_end_hours', 'entry_hours',
    'anaphase_hours', 'exit_without_anaphase_hours', 'death_hours', 'loss_hours')
SECONDARY_KINDS = ('death', 'tracking_loss', 'division', 'cytokinesis_failure', 'other')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def fields(value, expected, context):
    require(type(value) is dict and set(value) == set(expected),
            context + ': missing or unknown fields')


def text(value, context):
    require(type(value) is str and bool(value.strip()), context + ': expected nonempty text')


def nullable_text(value, context):
    if value is not None:
        text(value, context)


def number(value, context, minimum=0):
    require(type(value) in (int, float) and math.isfinite(value) and value >= minimum,
            context + ': expected finite nonnegative number')


def instant(value, context):
    text(value, context)
    try:
        result = datetime.fromisoformat(value.replace('Z', '+00:00'))
    except ValueError as error:
        raise ValueError(context + ': invalid timestamp') from error
    require(result.tzinfo is not None and result.utcoffset().total_seconds() == 0,
            context + ': timestamp must explicitly use UTC')
    return result


def well_key(row):
    return (row['experiment_id'], row['well_id'])


def founder_key(row):
    return well_key(row) + (row['founder_id'],)


def analysable(well):
    return well['expected_founders'] is not None and well['acquisition_status'] in ('acquired', 'failed')


def strict_json(data):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, 'duplicate JSON object key: ' + key)
            result[key] = value
        return result

    def bad_constant(value):
        raise ValueError('nonfinite JSON constant: ' + value)

    value = json.loads(data, object_pairs_hook=pairs, parse_constant=bad_constant)

    def finite(item):
        if type(item) is float:
            require(math.isfinite(item), 'nonfinite JSON number')
        elif type(item) is dict:
            for child in item.values():
                finite(child)
        elif type(item) is list:
            for child in item:
                finite(child)
    finite(value)
    return value


def load_pinned(path, expected):
    require(type(expected) is str and re.fullmatch('[0-9a-f]{64}', expected), 'invalid SHA256 pin')
    path = Path(path)
    require(path.is_file() and not path.is_symlink(), 'input must be a regular nonsymlink file: ' + str(path))
    before = path.stat()
    data = path.read_bytes()
    after = path.stat()
    require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), 'input changed during read')
    digest = hashlib.sha256(data).hexdigest()
    require(digest == expected, 'input SHA256 mismatch: ' + str(path))
    return strict_json(data), dict(path=str(path.resolve()), bytes=len(data), sha256=digest)


def table(value, name):
    fields(value, ('schema_version', name), name)
    require(type(value['schema_version']) is int and value['schema_version'] == 1, 'unsupported schema version')
    require(type(value[name]) is list, name + ': expected list')
    return value[name]


def validate_design(design, roster_sha):
    fields(design, ('schema_version', 'data_kind', 'acute_end_hours', 'secondary_end_hours',
        'time_unit', 'roster_frozen_at_utc', 'roster_sha256', 'scoring_frozen_at_utc',
        'adjudication_frozen', 'scoring_rule_id', 'error_criteria', 'source_acquisition_id'), 'design')
    require(type(design['schema_version']) is int and design['schema_version'] == 1, 'unsupported design schema')
    require(design['data_kind'] in ('synthetic_control', 'observed'), 'unknown data kind')
    require(design['time_unit'] == 'hours_after_time_zero', 'unsupported time unit')
    number(design['acute_end_hours'], 'acute_end_hours')
    require(design['acute_end_hours'] > 0, 'acute window must be positive')
    number(design['secondary_end_hours'], 'secondary_end_hours')
    require(design['secondary_end_hours'] >= design['acute_end_hours'], 'secondary window precedes acute endpoint')
    require(design['adjudication_frozen'] is True, 'adjudication must be explicitly frozen')
    require(design['roster_sha256'] == roster_sha, 'design does not pin the separately supplied roster')
    for name in ('scoring_rule_id', 'source_acquisition_id'):
        text(design[name], name)
    require(type(design['error_criteria']) is dict and bool(design['error_criteria']), 'error criteria must not be empty')
    for label, description in design['error_criteria'].items():
        text(label, 'error label')
        text(description, 'error definition')
    return (instant(design['roster_frozen_at_utc'], 'roster freeze'),
            instant(design['scoring_frozen_at_utc'], 'scoring freeze'))


def validate_outcome(row, founder, design):
    fields(row, OUTCOME_FIELDS, 'outcome')
    require(row['parent_id'] is None, 'descendant cannot be a primary outcome')
    status = row['status']
    require(status in STATUSES, 'unknown primary status')
    end = design['acute_end_hours']
    for name in TIME_FIELDS:
        value = row[name]
        if value is not None:
            number(value, name)
            require(value <= end, name + ': primary time beyond acute endpoint')
    t, tracked = row['classification_hours'], row['observation_end_hours']
    require(t is not None and tracked is not None and tracked >= t, 'missing or inconsistent primary follow-up')
    require(type(row['observed_through_window']) is bool, 'endpoint observation must be explicit true/false')
    require(not row['observed_through_window'] or tracked == end, 'endpoint observation contradicts follow-up time')
    entry, anaphase, slip = row['entry_hours'], row['anaphase_hours'], row['exit_without_anaphase_hours']
    death, loss = row['death_hours'], row['loss_hours']
    baseline = founder['baseline_mitotic_state']
    nullable_text(row['reason'], 'reason')
    nullable_text(row['entry_unknown_reason'], 'entry_unknown_reason')
    require(row['death_stage'] in (None, 'interphase', 'mitotic', 'unknown'), 'unknown death stage label')
    require(type(row['error_types']) is list and all(type(x) is str for x in row['error_types']), 'error types must be a list of labels')
    require(len(set(row['error_types'])) == len(row['error_types']), 'duplicate error types')
    require(all(x in design['error_criteria'] for x in row['error_types']), 'unknown error type')
    if entry is not None:
        require(entry <= t and row['entry_unknown_reason'] is None, 'entry timing or unknown-entry reason inconsistent')
        require(baseline != 'mitotic', 'already-mitotic baseline cannot have an invented postbaseline entry time')
    entered = status in ('A', 'U', 'S', 'O') or status == 'D' and row['death_stage'] == 'mitotic'
    if entered and entry is None:
        text(row['entry_unknown_reason'], 'unknown entry time needs a reason')
        require(baseline in ('mitotic', 'unknown'), 'interphase baseline needs an observed entry for a mitotic outcome')
    if not entered and entry is None and baseline == 'mitotic' and status in ('D', 'L'):
        text(row['entry_unknown_reason'], 'baseline-mitotic founder has left-censored entry')
    elif not entered and entry is None:
        require(row['entry_unknown_reason'] is None, 'entry reason supplied without a classified mitotic outcome')
    for event in (anaphase, slip):
        if event is not None and entry is not None:
            require(event >= entry, 'outcome precedes mitotic entry')
    if status in ('A', 'U'):
        require(anaphase == t and anaphase is not None and slip is None and death is None and loss is None,
                'A/U needs first anaphase as the classification event; later fate belongs in secondary')
    elif status == 'S':
        require(slip == t and slip is not None and anaphase is None and death is None and loss is None,
                'S needs exit without anaphase as the classification event')
        text(row['reason'], 'S adjudication reason')
    elif status == 'D':
        require(death == t and death is not None and tracked == t and loss is None,
                'D requires death as terminal primary follow-up')
        require(anaphase is None and slip is None, 'death cannot overwrite an earlier A/U/S event')
        require(row['death_stage'] in ('interphase', 'mitotic', 'unknown'), 'D needs an explicit death stage, including unknown')
        if row['death_stage'] == 'interphase':
            require(entry is None and baseline != 'mitotic', 'interphase death contradicts mitotic entry/baseline')
    elif status == 'L':
        require(loss == t and loss is not None and tracked == t and death is None,
                'L requires loss time as terminal primary follow-up')
        require(anaphase is None and slip is None, 'loss cannot overwrite an earlier A/U/S event')
        text(row['reason'], 'L loss reason')
    else:
        require(t == end and tracked == end and row['observed_through_window'] is True and anaphase is None and slip is None and death is None and loss is None,
                'B/O requires actual endpoint follow-up without a prior classified outcome')
        if status == 'B':
            require(entry is None and baseline == 'interphase', 'B requires an adjudicated interphase baseline and no entry')
    if status != 'D':
        require(row['death_stage'] is None, 'death stage is only a primary D field')
    if status == 'A':
        require(type(row['visible_error']) is bool, 'A requires explicit true/false visible-error status')
        require(row['visible_error'] == bool(row['error_types']), 'A error boolean and labels disagree')
    else:
        require(row['visible_error'] is None and not row['error_types'], 'only A may have definite scored primary errors')
    if status == 'U':
        text(row['reason'], 'U unscorable reason')
    secondary = row['secondary']
    require(secondary is None or type(secondary) is list, 'secondary must be null or a list')
    event_ids = set()
    for event in secondary or []:
        fields(event, ('event_id', 'subject_id', 'parent_id', 'kind', 'time_hours', 'observed', 'note'), 'secondary event')
        for name in ('event_id', 'subject_id', 'note'):
            text(event[name], 'secondary ' + name)
        nullable_text(event['parent_id'], 'secondary parent_id')
        require(event['event_id'] not in event_ids, 'duplicate secondary event ID')
        event_ids.add(event['event_id'])
        require(event['kind'] in SECONDARY_KINDS and type(event['observed']) is bool, 'invalid secondary kind/boolean')
        number(event['time_hours'], 'secondary time')
        require(t < event['time_hours'] <= design['secondary_end_hours'], 'secondary observation must follow primary classification within its window')
        require(status not in ('D', 'L'), 'secondary follow-up cannot silently resume primary death or tracking loss')
        if event['subject_id'] == row['founder_id']:
            require(event['parent_id'] is None, 'founder secondary event cannot have a parent')
        else:
            text(event['parent_id'], 'descendant secondary event needs a parent ID')
    return row


def fraction(numerator, denominator):
    return None if denominator == 0 else numerator / denominator


def summarize(rows, n0):
    counts = {status: sum(row['status'] == status for row in rows) for status in STATUSES}
    require(sum(counts.values()) == n0, 'founder census does not reconcile')
    errors = sum(row['status'] == 'A' and row['visible_error'] is True for row in rows)
    a, u = counts['A'], counts['U']
    secondary_records = [event for row in rows for event in row['secondary'] or []]
    return dict(N0=n0, counts=counts, E=errors,
        error_fraction_among_scorable_anaphases=fraction(errors, a),
        observed_yield_without_scored_error=fraction(a-errors, n0),
        scorable_anaphase_fraction=fraction(a, n0),
        status_proportions={status: fraction(count, n0) for status, count in counts.items()},
        observed_anaphase_error_bounds=dict(lower=fraction(errors, a+u), upper=fraction(errors+u, a+u),
            denominator=a+u, excludes_tracking_loss=True),
        primary_classification_completeness=fraction(n0-counts['L'], n0),
        primary_times=[dict(founder_id=row['founder_id'], well_id=row['well_id'], status=row['status'],
            entry_hours=row['entry_hours'], classification_hours=row['classification_hours']) for row in rows],
        secondary=dict(founders_with_unavailable_secondary=sum(row['secondary'] is None for row in rows),
            founders_with_secondary_ledger=sum(row['secondary'] is not None for row in rows),
            observed_event_records={kind: (sum(e['kind'] == kind and e['observed'] is True for e in secondary_records)
                if any(e['kind'] == kind for e in secondary_records) else None)
                                    for kind in SECONDARY_KINDS},
            explicit_absence_records={kind: (sum(e['kind'] == kind and e['observed'] is False for e in secondary_records)
                if any(e['kind'] == kind for e in secondary_records) else None)
                                      for kind in SECONDARY_KINDS},
            scope='Counts of supplied secondary observation records, not incidence or complete ascertainment.'))


def analyze(design, wells_input, founders_input, outcomes_input, roster_sha):
    roster_freeze, scoring_freeze = validate_design(design, roster_sha)
    wells, founders, outcomes = (table(wells_input, 'wells'), table(founders_input, 'founders'), table(outcomes_input, 'outcomes'))
    registry, clone_map = {}, {}
    for row in wells:
        fields(row, WELL_FIELDS, 'well')
        for name in WELL_FIELDS:
            if name not in ('editing_family_id', 'expected_founders', 'failure_reason'):
                text(row[name], 'well ' + name)
        key = well_key(row)
        require(key not in registry, 'duplicate well key')
        registry[key] = row
        require(row['assay_scope'] in ('spontaneous', 'challenge', 'assay_control'), 'unknown assay scope')
        require(row['material_origin'] in ('engineered', 'patient_derived', 'other'), 'unknown material origin')
        nullable_text(row['editing_family_id'], 'editing family')
        if row['material_origin'] == 'engineered':
            text(row['editing_family_id'], 'engineered model needs actual editing family')
        else:
            require(row['editing_family_id'] is None, 'non-engineered material must not acquire an editing family')
        if row['expected_founders'] is not None:
            require(type(row['expected_founders']) is int and row['expected_founders'] >= 0, 'invalid expected founder count')
        require(row['acquisition_status'] in ('acquired', 'failed', 'unrun'), 'unknown acquisition status')
        if row['acquisition_status'] == 'acquired':
            require(row['expected_founders'] is not None and row['failure_reason'] is None,
                    'acquired well needs known starting count and no failure reason')
        else:
            text(row['failure_reason'], 'failed/unrun well needs a reason')
            if row['acquisition_status'] == 'unrun':
                require(row['expected_founders'] is None, 'unrun well cannot have an observed starting census')
        identity = tuple(row[name] for name in ('material_origin', 'source_background_id', 'editing_family_id', 'model_id'))
        require(row['clone_id'] not in clone_map or clone_map[row['clone_id']] == identity, 'inconsistent clone hierarchy')
        clone_map[row['clone_id']] = identity
    require(registry, 'well registry must retain at least one planned well')
    roster, source_cells, by_well = {}, set(), defaultdict(list)
    for row in founders:
        fields(row, FOUNDER_FIELDS, 'founder')
        for name in FOUNDER_FIELDS:
            if name not in ('parent_id',):
                text(row[name], 'founder ' + name)
        require(row['parent_id'] is None, 'descendant cannot be a pretreatment founder')
        key, wk = founder_key(row), well_key(row)
        require(wk in registry and analysable(registry[wk]), 'founder refers to unknown well or unavailable starting census')
        require(key not in roster, 'duplicate roster founder')
        source_key = (row['raw_image_id'], row['field_id'], row['source_cell_id'])
        require(source_key not in source_cells, 'duplicate source cell represented by multiple founder IDs')
        source_cells.add(source_key)
        require(row['baseline_mitotic_state'] in ('interphase', 'mitotic', 'unknown'), 'unknown baseline mitotic state')
        start = instant(row['time_zero_utc'], 'time zero')
        require(roster_freeze <= start and scoring_freeze <= start, 'roster/scoring freeze must precede time zero')
        roster[key] = row
        by_well[wk].append(row)
    for key, row in registry.items():
        if analysable(row):
            require(len(by_well[key]) == row['expected_founders'], 'roster and independent declared starting count differ')
    scored, outcome_by_well = {}, defaultdict(list)
    for row in outcomes:
        require(type(row) is dict and all(name in row for name in ('experiment_id', 'well_id', 'founder_id')), 'outcome identity missing')
        key = founder_key(row)
        require(key in roster, 'extra outcome or descendant absent from frozen roster')
        require(key not in scored, 'duplicate primary outcome')
        scored[key] = validate_outcome(row, roster[key], design)
        outcome_by_well[well_key(row)].append(row)
    require(set(scored) == set(roster), 'missing outcome for a frozen pretreatment founder')
    per_well, groups = [], defaultdict(list)
    for key, well in registry.items():
        summary = summarize(outcome_by_well[key], len(by_well[key])) if analysable(well) else None
        per_well.append(dict(well=well, analysis=summary))
        groups[tuple(well[name] for name in GROUP_FIELDS)].append(well)
    block_conditions = []
    for group, members in groups.items():
        analyzed = [well for well in members if analysable(well)]
        rows = [row for well in analyzed for row in outcome_by_well[well_key(well)]]
        n0 = sum(len(by_well[well_key(well)]) for well in analyzed)
        block_conditions.append(dict(stratum=dict(zip(GROUP_FIELDS, group)),
            planned_wells=len(members), analysed_wells=len(analyzed),
            acquired_wells=sum(w['acquisition_status'] == 'acquired' for w in members),
            failed_wells=sum(w['acquisition_status'] == 'failed' for w in members),
            unrun_wells=sum(w['acquisition_status'] == 'unrun' for w in members),
            well_ids=[well['well_id'] for well in members],
            analysis=summarize(rows, n0) if analyzed else None,
            scope='Descriptive sum of wells with known starting censuses in this exact retained hierarchy, including interrupted acquisitions; not independent cells or efficacy inference.'))
    return dict(schema_version=1, data_kind=design['data_kind'], design=design,
        roster_founders=len(roster), planned_wells=len(registry), per_well=per_well,
        block_condition_descriptive_sums=block_conditions,
        adjudicated_founder_records=[dict(roster=roster[key], outcome=scored[key]) for key in sorted(roster)],
        limitations=['Validates supplied records, not raw-image truth, actual pretreatment freezing, biological independence or roster completeness.',
            'Yield without a scored visible error is not a normal-karyotype or daughter-function result.',
            'U bounds exclude L and do not solve informative loss, survival or entry selection.',
            'No pooled-cell efficacy, confidence intervals, p-values, power, automatic advance or clinical inference.'])


def run_cli(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for label in ('design', 'wells', 'founders', 'outcomes'):
        parser.add_argument('--' + label, type=Path, required=True)
        parser.add_argument('--' + label + '-sha256', required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(argv)
    require(not args.output_dir.exists(), 'output directory already exists; refusing overwrite')
    data, provenance = {}, {}
    for label in ('design', 'wells', 'founders', 'outcomes'):
        data[label], provenance[label] = load_pinned(getattr(args, label), getattr(args, label + '_sha256'))
    result = analyze(data['design'], data['wells'], data['founders'], data['outcomes'], provenance['founders']['sha256'])
    source = Path(__file__).resolve()
    result['provenance'] = dict(created_at_utc=datetime.now(timezone.utc).isoformat(), inputs=provenance,
        source=dict(path=str(source), bytes=source.stat().st_size, sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
    encoded = (json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()
    args.output_dir.mkdir(parents=False, exist_ok=False)
    with (args.output_dir / 'summary.json').open('xb') as stream:
        stream.write(encoded)
    receipt = dict(status='complete', summary_sha256=hashlib.sha256(encoded).hexdigest(),
        data_kind=result['data_kind'], roster_founders=result['roster_founders'], planned_wells=result['planned_wells'])
    with (args.output_dir / 'status.json').open('x') as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    print(json.dumps(receipt, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    try:
        run_cli()
    except (ValueError, OSError, TypeError, KeyError) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        raise SystemExit(2)
