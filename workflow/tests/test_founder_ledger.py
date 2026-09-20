"""Fictional records only; no patient or laboratory data are used."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / 'scripts' / 'analyze_founder_ledger.py'
SPEC = importlib.util.spec_from_file_location('founder_ledger', SCRIPT)
ledger = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ledger)


def fictional_inputs(counts=None):
    counts = counts or {'A': 1}
    well = dict(experiment_id='fictional-experiment', well_id='W1', block_id='B1',
        masked_condition_id='masked-01', assay_id='spontaneous-imaging-v1', assay_scope='spontaneous',
        model_id='masked-model-01', material_origin='engineered', source_background_id='fictional-background',
        editing_family_id='F1', clone_id='C1', preparation_id='P1', experimental_day='fictional-day-1',
        expected_founders=sum(counts.values()), acquisition_status='acquired', failure_reason=None)
    founders, outcomes = [], []
    for status, count in counts.items():
        for _ in range(count):
            ident = 'founder-%03d' % (len(founders)+1)
            founders.append(dict(experiment_id=well['experiment_id'], well_id='W1', founder_id=ident,
                field_id='field-01', source_cell_id='cell-'+ident, raw_image_id='fictional-image-01', parent_id=None,
                time_zero_utc='2026-01-01T01:00:00Z', baseline_mitotic_state='interphase'))
            entered = status in ('A', 'U', 'S', 'O')
            outcomes.append(dict(experiment_id=well['experiment_id'], well_id='W1', founder_id=ident,
                parent_id=None, status=status, classification_hours=48 if status in ('B', 'O') else 12,
                observation_end_hours=48 if status in ('B', 'O') else 12,
                observed_through_window=status in ('B', 'O'), entry_hours=6 if entered else None,
                entry_unknown_reason=None, anaphase_hours=12 if status in ('A', 'U') else None,
                exit_without_anaphase_hours=12 if status == 'S' else None,
                death_hours=12 if status == 'D' else None, loss_hours=12 if status == 'L' else None,
                death_stage='interphase' if status == 'D' else None,
                reason='fictional adjudication reason' if status in ('U', 'S', 'L') else None,
                visible_error=False if status == 'A' else None, error_types=[], secondary=None))
    design = dict(schema_version=1, data_kind='synthetic_control', acute_end_hours=48,
        secondary_end_hours=120, time_unit='hours_after_time_zero',
        roster_frozen_at_utc='2026-01-01T00:00:00Z', roster_sha256='a'*64,
        scoring_frozen_at_utc='2026-01-01T00:00:00Z', adjudication_frozen=True,
        scoring_rule_id='fictional-criteria-v1',
        error_criteria={'lagging':'Fictional lagging-material rule', 'bridge':'Fictional bridge rule'},
        source_acquisition_id='fictional-acquisition-only')
    return dict(design=design, wells=dict(schema_version=1, wells=[well]),
        founders=dict(schema_version=1, founders=founders), outcomes=dict(schema_version=1, outcomes=outcomes))


def analyze(data):
    return ledger.analyze(data['design'], data['wells'], data['founders'], data['outcomes'], 'a'*64)


def summary(data):
    return analyze(data)['per_well'][0]['analysis']


class FounderLedgerTests(unittest.TestCase):
    def test_protocol_counterexample_does_not_hide_division_loss(self):
        values = []
        for a, e in ((80, 16), (20, 2)):
            data = fictional_inputs({'A':a, 'B':100-a})
            for row in data['outcomes']['outcomes'][:e]:
                row.update(visible_error=True, error_types=['lagging'])
            result = summary(data)
            self.assertEqual(result['N0'], 100)
            values.append((result['error_fraction_among_scorable_anaphases'], result['observed_yield_without_scored_error']))
        self.assertEqual(values, [(0.2, 0.64), (0.1, 0.18)])

    def test_all_u_all_l_and_empty_are_distinct(self):
        u, lost, empty = (summary(fictional_inputs(counts)) for counts in ({'U':3}, {'L':3}, {'A':0}))
        self.assertIsNone(u['error_fraction_among_scorable_anaphases'])
        self.assertEqual((u['observed_anaphase_error_bounds']['lower'], u['observed_anaphase_error_bounds']['upper']), (0, 1))
        self.assertIsNone(lost['observed_anaphase_error_bounds']['lower'])
        self.assertEqual(lost['primary_classification_completeness'], 0)
        self.assertEqual(empty['N0'], 0)
        self.assertIsNone(empty['observed_yield_without_scored_error'])
        self.assertTrue(all(value is None for value in empty['status_proportions'].values()))

    def test_all_statuses_and_multierror_count_once(self):
        data = fictional_inputs(dict.fromkeys(ledger.STATUSES, 1))
        data['outcomes']['outcomes'][0].update(visible_error=True, error_types=['lagging','bridge'])
        result = summary(data)
        self.assertEqual(result['counts'], dict.fromkeys(ledger.STATUSES, 1))
        self.assertEqual(result['E'], 1)
        self.assertEqual(result['observed_anaphase_error_bounds']['lower'], .5)
        self.assertEqual(result['observed_anaphase_error_bounds']['upper'], 1)

    def test_tracking_loss_does_not_enter_u_bounds(self):
        results = []
        for losses in (0, 10):
            data = fictional_inputs({'A':2,'U':2,'L':losses})
            data['outcomes']['outcomes'][0].update(visible_error=True,error_types=['bridge'])
            results.append(summary(data))
        self.assertEqual(results[0]['observed_anaphase_error_bounds'], results[1]['observed_anaphase_error_bounds'])
        self.assertGreater(results[0]['primary_classification_completeness'], results[1]['primary_classification_completeness'])

    def test_missing_extra_duplicate_and_descendant_rows_refused(self):
        for kind in ('missing','extra','duplicate','descendant','roster_duplicate','source_duplicate'):
            with self.subTest(kind=kind):
                data = fictional_inputs({'A':2})
                if kind == 'missing':data['outcomes']['outcomes'].pop()
                if kind == 'extra':
                    row = copy.deepcopy(data['outcomes']['outcomes'][0]); row['founder_id']='daughter-new'
                    data['outcomes']['outcomes'].append(row)
                if kind == 'duplicate':data['outcomes']['outcomes'].append(copy.deepcopy(data['outcomes']['outcomes'][0]))
                if kind == 'descendant':data['founders']['founders'][0]['parent_id']='some-parent'
                if kind == 'roster_duplicate':data['founders']['founders'].append(copy.deepcopy(data['founders']['founders'][0]))
                if kind == 'source_duplicate':data['founders']['founders'][1]['source_cell_id']=data['founders']['founders'][0]['source_cell_id']
                with self.assertRaises(ValueError):analyze(data)

    def test_declared_starting_count_is_independent_authority(self):
        data = fictional_inputs({'A':2})
        data['outcomes']['outcomes'].pop(); data['founders']['founders'].pop()
        with self.assertRaisesRegex(ValueError, 'starting count'):analyze(data)

    def test_endpoint_tracking_required_for_b_and_o(self):
        for status in ('B','O'):
            for mutation in ({'observation_end_hours':47},{'observed_through_window':False},{'classification_hours':47}):
                with self.subTest(status=status,mutation=mutation):
                    data=fictional_inputs({status:1}); data['outcomes']['outcomes'][0].update(mutation)
                    with self.assertRaises(ValueError):analyze(data)

    def test_baseline_mitotic_entry_is_left_censored(self):
        for status in ('A','U','S','O','D','L'):
            with self.subTest(status=status):
                data=fictional_inputs({status:1}); row=data['outcomes']['outcomes'][0]
                data['founders']['founders'][0]['baseline_mitotic_state']='mitotic'
                row.update(entry_hours=None, entry_unknown_reason='Already mitotic at time zero')
                if status=='D':row['death_stage']='mitotic'
                self.assertEqual(summary(data)['counts'][status],1)
                row['entry_hours']=0
                with self.assertRaises(ValueError):analyze(data)

    def test_unknown_baseline_cannot_establish_no_entry_b(self):
        data=fictional_inputs({'B':1})
        data['founders']['founders'][0]['baseline_mitotic_state']='unknown'
        with self.assertRaisesRegex(ValueError,'interphase baseline'):analyze(data)
        self.assertEqual(data['outcomes']['outcomes'][0]['status'],'B')
        for status in ('A','U','S','O'):
            data=fictional_inputs({status:1})
            data['founders']['founders'][0]['baseline_mitotic_state']='unknown'
            data['outcomes']['outcomes'][0].update(entry_hours=None,
                entry_unknown_reason='Fictional baseline uncertain; category directly adjudicated')
            self.assertEqual(summary(data)['counts'][status],1)

    def test_early_death_loss_and_unknown_death_stage(self):
        for status in ('D','L'):
            data=fictional_inputs({status:1}); row=data['outcomes']['outcomes'][0]
            row.update(classification_hours=1,observation_end_hours=1)
            row['death_hours' if status=='D' else 'loss_hours']=1
            self.assertEqual(summary(data)['counts'][status],1)
            row['entry_hours']=2
            with self.assertRaises(ValueError):analyze(data)
        data=fictional_inputs({'D':1});data['outcomes']['outcomes'][0]['death_stage']='unknown'
        self.assertEqual(summary(data)['counts']['D'],1)

    def test_fixed_outcome_cannot_be_overwritten_by_death_or_loss(self):
        for status in ('D','L'):
            for evidence in ('anaphase_hours','exit_without_anaphase_hours'):
                data=fictional_inputs({status:1});data['outcomes']['outcomes'][0][evidence]=6
                with self.assertRaises(ValueError):analyze(data)

    def test_secondary_death_retains_a_and_descendant_is_not_new_founder(self):
        data=fictional_inputs({'A':1});row=data['outcomes']['outcomes'][0]
        row['secondary']=[dict(event_id='s1',subject_id=row['founder_id'],parent_id=None,
            kind='death',time_hours=24,observed=True,note='Fictional later death'),
            dict(event_id='s2',subject_id='daughter-1',parent_id=row['founder_id'],kind='division',
            time_hours=36,observed=True,note='Fictional daughter division')]
        result=summary(data)
        self.assertEqual((result['N0'],result['counts']['A'],result['counts']['D']),(1,1,0))
        self.assertEqual(result['secondary']['observed_event_records']['death'],1)
        row['secondary'][0]['time_hours']=6
        with self.assertRaises(ValueError):analyze(data)

    def test_error_booleans_labels_and_reasons_are_strict(self):
        for status,mutation in [('A',{'visible_error':None}),('A',{'visible_error':0}),
            ('A',{'visible_error':True}),('A',{'error_types':['unknown']}),
            ('U',{'visible_error':False}),('U',{'reason':None}),('S',{'reason':None}),('L',{'reason':None})]:
            with self.subTest(status=status,mutation=mutation):
                data=fictional_inputs({status:1});data['outcomes']['outcomes'][0].update(mutation)
                with self.assertRaises(ValueError):analyze(data)

    def test_crossed_preparations_are_valid_but_clone_mapping_is_consistent(self):
        data=fictional_inputs({'A':0});first=data['wells']['wells'][0]
        for family,clone in [('F2','C2'),('F3','C3')]:
            second=copy.deepcopy(first);second.update(well_id='W'+family,editing_family_id=family,clone_id=clone)
            data['wells']['wells'].append(second)
        result=analyze(data)
        self.assertEqual(len(result['block_condition_descriptive_sums']),3)
        self.assertEqual({r['stratum']['preparation_id'] for r in result['block_condition_descriptive_sums']},{'P1'})
        data['wells']['wells'][1]['clone_id']='C1'
        with self.assertRaisesRegex(ValueError,'clone hierarchy'):analyze(data)

    def test_patient_origin_and_assay_scope_are_not_pooled(self):
        data=fictional_inputs({'A':0});first=data['wells']['wells'][0]
        second=copy.deepcopy(first);second.update(well_id='W2',material_origin='patient_derived',editing_family_id=None,clone_id='patient-clone')
        third=copy.deepcopy(first);third.update(well_id='W3',assay_scope='challenge',assay_id='challenge-01')
        data['wells']['wells'].extend([second,third])
        self.assertEqual(len(analyze(data)['block_condition_descriptive_sums']),3)
        second['editing_family_id']='invented-family'
        with self.assertRaises(ValueError):analyze(data)

    def test_unequal_technical_wells_keep_real_denominators(self):
        data=fictional_inputs({'A':2,'B':1});other=fictional_inputs({'A':1,'B':3})
        data['outcomes']['outcomes'][0].update(visible_error=True,error_types=['bridge'])
        other['wells']['wells'][0]['well_id']='W2'
        for row in other['founders']['founders']:
            row.update(well_id='W2',raw_image_id='fictional-image-02')
        for row in other['outcomes']['outcomes']:row['well_id']='W2'
        for name in ('wells','founders','outcomes'):data[name][name].extend(other[name][name])
        result=analyze(data)
        self.assertEqual([row['analysis']['N0'] for row in result['per_well']],[3,4])
        self.assertEqual(len(result['block_condition_descriptive_sums']),1)
        combined=result['block_condition_descriptive_sums'][0]['analysis']
        self.assertEqual((combined['N0'],combined['counts']['A'],combined['E']),(7,3,1))
        self.assertAlmostEqual(combined['error_fraction_among_scorable_anaphases'],1/3)
        self.assertAlmostEqual(combined['observed_yield_without_scored_error'],2/7)

    def test_source_identity_is_field_local_but_not_relabelled_by_well(self):
        data=fictional_inputs({'A':2})
        first,second=data['founders']['founders']
        second.update(source_cell_id=first['source_cell_id'],field_id='field-02')
        self.assertEqual(summary(data)['N0'],2)
        second.update(field_id=first['field_id'],well_id='W2')
        data['wells']['wells'][0]['expected_founders']=1
        well=copy.deepcopy(data['wells']['wells'][0]);well['well_id']='W2'
        data['wells']['wells'].append(well)
        data['outcomes']['outcomes'][1]['well_id']='W2'
        with self.assertRaisesRegex(ValueError,'source cell'):analyze(data)

    def test_failed_and_unrun_are_not_valid_empty_censuses(self):
        data=fictional_inputs({'A':0})
        for status in ('failed','unrun'):
            well=copy.deepcopy(data['wells']['wells'][0]);well.update(well_id=status,
                expected_founders=None,acquisition_status=status,failure_reason='Fictional unavailable acquisition')
            data['wells']['wells'].append(well)
        result=analyze(data)
        self.assertEqual(result['per_well'][0]['analysis']['N0'],0)
        self.assertIsNone(result['per_well'][1]['analysis'])
        self.assertIsNone(result['per_well'][2]['analysis'])
        data['wells']['wells'][2]['expected_founders']=0
        with self.assertRaises(ValueError):analyze(data)

    def test_failed_after_known_roster_keeps_a_and_l(self):
        data=fictional_inputs({'A':1,'L':2})
        data['wells']['wells'][0].update(acquisition_status='failed',failure_reason='Fictional acquisition interrupted after census')
        result=summary(data)
        self.assertEqual(result['N0'],3)
        self.assertEqual((result['counts']['A'],result['counts']['L']),(1,2))

    def test_nonfinite_boolean_and_disordered_times_refused(self):
        for value in (float('nan'),float('inf'),True,-1,49):
            data=fictional_inputs();data['outcomes']['outcomes'][0]['classification_hours']=value
            with self.assertRaises(ValueError):analyze(data)
        data=fictional_inputs();data['outcomes']['outcomes'][0]['entry_hours']=13
        with self.assertRaises(ValueError):analyze(data)
        data=fictional_inputs();data['design']['roster_frozen_at_utc']='2026-01-02T00:00:00Z'
        with self.assertRaises(ValueError):analyze(data)
        for raw in ('{"x":NaN}','{"x":Infinity}','{"x":1e9999}','{"x":1,"x":2}'):
            with self.assertRaises(ValueError):ledger.strict_json(raw)

    def test_cli_hash_provenance_and_new_output_fail_closed(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); data=fictional_inputs({'A':1,'U':1,'L':1});paths={};pins={}
            for name in ('founders','wells','outcomes'):
                paths[name]=root/(name+'.json');paths[name].write_text(json.dumps(data[name]))
                pins[name]=hashlib.sha256(paths[name].read_bytes()).hexdigest()
            data['design']['roster_sha256']=pins['founders']
            paths['design']=root/'design.json';paths['design'].write_text(json.dumps(data['design']))
            pins['design']=hashlib.sha256(paths['design'].read_bytes()).hexdigest()
            command=[sys.executable,str(SCRIPT)]
            for name in ('design','wells','founders','outcomes'):
                command.extend(['--'+name,str(paths[name]),'--'+name+'-sha256',pins[name]])
            output=root/'result';command.extend(['--output-dir',str(output)])
            run=subprocess.run(command,capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stderr)
            result=json.loads((output/'summary.json').read_text())
            self.assertEqual(result['provenance']['inputs']['founders']['sha256'],pins['founders'])
            self.assertEqual(result['roster_founders'],3)
            self.assertEqual(subprocess.run(command,capture_output=True).returncode,2)
            bad=list(command);bad[bad.index('--founders-sha256')+1]='0'*64;bad[-1]=str(root/'rejected')
            self.assertEqual(subprocess.run(bad,capture_output=True).returncode,2)
            self.assertFalse((root/'rejected').exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
