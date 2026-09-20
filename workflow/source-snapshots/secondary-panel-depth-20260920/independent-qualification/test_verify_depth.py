"""Fictional checks only; never reads project case data."""
from array import array
from collections import Counter
import copy
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

SOURCE = Path(__file__).with_name('verify-depth.py')
spec = importlib.util.spec_from_file_location('depth_audit', SOURCE)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)

GFF = '''##gff-version 3
1\tx\tgene\t2\t12\t.\t+\t.\tID=gene:g1;Name=ONE;biotype=protein_coding
1\tx\tmRNA\t2\t12\t.\t+\t.\tID=transcript:t1;Parent=gene:g1
1\tx\tCDS\t2\t5\t.\t+\t0\tParent=transcript:t1
1\tx\tCDS\t8\t10\t.\t+\t2\tParent=transcript:t1
1\tx\tmRNA\t2\t12\t.\t+\t.\tID=transcript:t2;Parent=gene:g1
1\tx\tCDS\t2\t5\t.\t+\t0\tParent=transcript:t2
1\tx\tgene\t4\t6\t.\t-\t.\tID=gene:g2;Name=TWO;biotype=protein_coding
1\tx\tmRNA\t4\t6\t.\t-\t.\tID=transcript:t3;Parent=gene:g2
1\tx\tCDS\t4\t6\t.\t-\t1\tParent=transcript:t3
'''


def write_gff(path):
    with gzip.open(path, 'wt') as stream:
        stream.write(GFF)


def fixture(directory):
    path = directory / 'annotation.gff3.gz'
    write_gff(path)
    design = audit.independent_targets(path, ['ONE', 'TWO'], {'1': 60})
    return design


class TargetsAndRows(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.design = fixture(Path(self.temp.name))

    def test_duplicate_isoform_CDS_and_overlap_count_panel_once(self):
        self.assertEqual(self.design['panel_union_bases'], 8)
        self.assertEqual(self.design['gene_union_bases_sum'], 10)
        self.assertEqual([g['CDS_union_bases'] for g in self.design['genes']], [7, 3])
        self.assertEqual(self.design['windows'], [
            dict(chrom='1', start=1, end=3, gene_ids=['gene:g1']),
            dict(chrom='1', start=3, end=5, gene_ids=['gene:g1', 'gene:g2']),
            dict(chrom='1', start=5, end=6, gene_ids=['gene:g2']),
            dict(chrom='1', start=7, end=10, gene_ids=['gene:g1'])])

    def test_annotation_missing_gene_refused(self):
        with self.assertRaisesRegex(ValueError, 'missing or ambiguous'):
            audit.independent_targets(Path(self.temp.name) / 'annotation.gff3.gz', ['ABSENT'], {'1': 60})

    def lines(self):
        return [f"1\t{p+1}\t{(p % 3)}\n" for w in self.design['windows'] for p in range(w['start'], w['end'])]

    def test_explicit_zero_and_non_ACGT_denominators(self):
        depths = audit.positions_from_lines(self.lines(), self.design['windows'])
        rows, groups = audit.aggregate(depths, self.design, ['AN', 'GT', 'N', 'AAA'])
        self.assertEqual(groups['panel_union']['span']['bases'], 8)
        self.assertEqual(groups['panel_union']['ACGT']['bases'], 6)
        self.assertEqual(groups['panel_union']['non_ACGT_reference_bases'], 2)
        self.assertEqual(groups['panel_union']['span']['zero_depth_bases'], 2)
        self.assertEqual(groups['gene:g1']['span']['bases'], 7)
        self.assertEqual(groups['gene:g2']['span']['bases'], 3)
        self.assertEqual(rows[0]['span_depth']['histogram'], {'1': 1, '2': 1})

    def test_missing_zero_row_refused(self):
        lines = self.lines()
        del lines[2]
        with self.assertRaises(ValueError):
            audit.positions_from_lines(lines, self.design['windows'])

    def test_extra_position_refused(self):
        with self.assertRaisesRegex(ValueError, 'extra'):
            audit.positions_from_lines(self.lines() + ['1\t40\t0\n'], self.design['windows'])

    def test_reordered_position_refused(self):
        lines = self.lines()
        lines[0], lines[1] = lines[1], lines[0]
        with self.assertRaisesRegex(ValueError, 'position/order'):
            audit.positions_from_lines(lines, self.design['windows'])

    def test_negative_position_depth_refused(self):
        lines = self.lines()
        lines[0] = '1\t2\t-1\n'
        with self.assertRaisesRegex(ValueError, 'negative'):
            audit.positions_from_lines(lines, self.design['windows'])

    def test_all_nonACGT_has_null_mean_not_zero_claim(self):
        d = audit.positions_from_lines(self.lines(), self.design['windows'])
        _, groups = audit.aggregate(d, self.design, ['NN', 'NN', 'N', 'NNN'])
        self.assertIsNone(groups['panel_union']['ACGT']['mean'])
        self.assertEqual(groups['panel_union']['ACGT']['bases'], 0)

    def test_cannot_hide_overlap_double_count(self):
        d = audit.positions_from_lines(self.lines(), self.design['windows'])
        _, groups = audit.aggregate(d, self.design, ['AA', 'AA', 'A', 'AAA'])
        bad = copy.deepcopy(groups)
        bad['panel_union']['span']['bases'] = 10
        with self.assertRaises(ValueError):
            audit.same(bad, groups, 'panel double counted')

    def test_first_overlap_assignment_preserves_gapped_window_reads(self):
        w = self.design['windows']
        self.assertEqual(audit.first_overlapping_window(0, 11, w), 0)
        self.assertEqual(audit.first_overlapping_window(5, 9, w), 2)
        self.assertIsNone(audit.first_overlapping_window(6, 7, w))


class PinsAndQuality(unittest.TestCase):
    def test_pin_content_mutation_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'input'
            p.write_bytes(b'one')
            pin = audit.Pins().add(p)
            p.write_bytes(b'two')
            with self.assertRaisesRegex(ValueError, 'sha256'):
                audit.Pins().add(p, pin)

    def test_pin_wrong_path_same_bytes_refused(self):
        with tempfile.TemporaryDirectory() as d:
            a, b = Path(d) / 'a', Path(d) / 'b'
            a.write_bytes(b'a'); b.write_bytes(b'a')
            with self.assertRaisesRegex(ValueError, 'path'):
                audit.Pins().add(b, audit.Pins().add(a))

    def test_input_modification_after_hash_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'input'; p.write_text('x')
            pins = audit.Pins(); pins.add(p); p.write_text('longer')
            with self.assertRaisesRegex(ValueError, 'changed during audit'):
                pins.stable()

    def test_nonfinite_duplicate_and_scalar_type_rejected(self):
        for text in ['{"x":NaN}', '{"x":1,"x":2}']:
            with self.assertRaises(ValueError): audit.parse(text)
        with self.assertRaises(ValueError): audit.same(True, 1, 'typed mismatch')

    def read(self, **kwargs):
        values = dict(flag=0, mapping_quality=60, query_qualities=[30, 30], query_length=2)
        values.update(kwargs)
        return types.SimpleNamespace(**values)

    def test_absent_quality_refused(self):
        with self.assertRaisesRegex(ValueError, 'unknown BQ'):
            audit.quality_available(self.read(query_qualities=None))

    def test_partial_unknown_quality_refused(self):
        with self.assertRaisesRegex(ValueError, 'unknown BQ'):
            audit.quality_available(self.read(query_qualities=[30, 255]))

    def test_unknown_mapping_quality_refused(self):
        with self.assertRaisesRegex(ValueError, 'unknown MQ'):
            audit.quality_available(self.read(mapping_quality=255))

    def test_known_zero_quality_and_excluded_unknown_allowed(self):
        audit.quality_available(self.read(query_qualities=[0, 0]))
        audit.quality_available(self.read(flag=1024, query_qualities=None, mapping_quality=255))


def save_fixture(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, indent=2)+'\n')
    return audit.Pins().add(path)


def public_fixture(root, config, design, sequences):
    directory = root / audit.PUBLIC_REL
    panel = {'total': len(config['genes']), 'rows': [{'symbol': s} for s in config['genes']]}
    panelpin = save_fixture(root/'panel.json', panel)
    inputs = {k:config['inputs'][k] for k in ('reference','reference_index','gff')}
    inputs['panel'] = panelpin
    public = copy.deepcopy(design)
    public.update(status='prepared',reads_inspected=False,inputs=inputs,panel_rows=panel['rows'])
    by_gene = {}
    for w,seq in zip(public['windows'],sequences):
        acgt = sum(x in 'ACGT' for x in seq)
        w.update(reference_ACGT_bases=acgt,reference_non_ACGT_bases=len(seq)-acgt,reference_sequence_sha256=hashlib.sha256(seq.encode()).hexdigest())
        for g in w['gene_ids']: by_gene[g]=by_gene.get(g,'')+seq
    for g in public['genes']:
        seq=by_gene[g['id']]; acgt=sum(x in 'ACGT' for x in seq)
        g.update(reference_ACGT_bases=acgt,reference_non_ACGT_bases=len(seq)-acgt,concatenated_CDS_union_sha256=hashlib.sha256(seq.encode()).hexdigest())
    summary=dict(status='complete',reads_inspected=False,genes=len(design['genes']),windows=len(design['windows']),panel_union_bases=design['panel_union_bases'],gene_union_bases_sum=design['gene_union_bases_sum'],panel_ACGT_bases=sum(w['reference_ACGT_bases'] for w in public['windows']),panel_non_ACGT_bases=sum(w['reference_non_ACGT_bases'] for w in public['windows']),gene_summary=[{k:g[k] for k in ('name','id','chrom','CDS_union_bases','reference_ACGT_bases','reference_non_ACGT_bases')} for g in public['genes']])
    outputs={n:save_fixture(directory/n,x) for n,x in [('reference-design.json',public),('target-summary.json',summary)]}
    (directory/'regions.bed').write_text(''.join(f"{w['chrom']}\t{w['start']}\t{w['end']}\n" for w in design['windows']))
    outputs['regions.bed']=audit.Pins().add(directory/'regions.bed')
    state=dict(status='complete',reads_inspected=False,inputs=inputs,outputs=outputs)
    pin=save_fixture(directory/'status.json',state)
    return public,summary,pin['sha256']


def launch_fixture(root, config, producer_state, producer_summary):
    run=root/audit.RUN_REL; launch=run.with_name(run.name+'-launch')
    inner=audit.Pins().add(root/audit.CONFIG_REL)
    outer=save_fixture(launch/'outer.json',dict(producer_configuration=config,producer_config_sha256=inner['sha256']))
    source_files={}
    for name in ['scripts/launch_secondary_depth_v2.py','scripts/qualify_secondary_depth_guard_v2.py','scripts/vendor/guard_helpers.py','tests/test_secondary_depth_guard_v2.py']:
        f=launch/'workflow'/name;f.parent.mkdir(parents=True,exist_ok=True);f.write_text('# fictional source only '+name+'\n');source_files[name]=audit.Pins().add(f)
    source=source_files['scripts/launch_secondary_depth_v2.py']
    qual=save_fixture(launch/'controls.json',dict(status='complete',tests_run=21,skipped=[],errors=0,failures=0,exit_code=0,source_files=source_files))
    results=dict(status=producer_state,summary=producer_summary)
    binding=dict(status='complete',before_case_access=True,producer_config=inner,inputs=dict(outer_config=outer,wrapper_source=source),producer_results=results)
    bound=save_fixture(launch/'worker-binding.json',binding)
    state=dict(status='complete',exit_code=0,worker_binding=bound,config=outer,source=source,qualification=qual,producer_status=producer_state)
    return save_fixture(launch/'status.json',state)['sha256']


class PublicAndLaunchBindings(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.design=fixture(self.root)
        self.sequences=['AN','GT','N','AAA']
        cfg=dict(genes=['ONE','TWO'],inputs={})
        for name in ('reference','reference_index','gff'):
            f=self.root/(name+'.txt');f.write_text(name);cfg['inputs'][name]=audit.Pins().add(f)
        self.cfg=cfg
        self.public,self.summary,self.sha=public_fixture(self.root,cfg,self.design,self.sequences)

    def test_public_reference_payload_reconciled(self):
        audit.reconcile_public_design(self.public,self.summary,self.design,self.sequences)
        with patch.object(audit,'PUBLIC_STATUS_SHA',self.sha),patch.object(audit,'EXPECTED_GENES',2):
            audit.public_bindings(self.root,audit.Pins(),self.cfg)

    def test_public_wrong_panel_scope_refused(self):
        with patch.object(audit,'PUBLIC_STATUS_SHA',self.sha):
            with self.assertRaisesRegex(ValueError,'panel count differs'):
                audit.public_bindings(self.root,audit.Pins(),self.cfg)

    def test_public_gene_subset_or_order_refused(self):
        for names in (['ONE'], ['TWO','ONE']):
            cfg=copy.deepcopy(self.cfg);cfg['genes']=names
            with patch.object(audit,'PUBLIC_STATUS_SHA',self.sha),patch.object(audit,'EXPECTED_GENES',2):
                with self.assertRaisesRegex(ValueError,'membership/order'):
                    audit.public_bindings(self.root,audit.Pins(),cfg)

    def test_public_status_wrong_bytes_refused(self):
        with patch.object(audit,'PUBLIC_STATUS_SHA','0'*64):
            with self.assertRaisesRegex(ValueError,'sha256'):
                audit.public_bindings(self.root,audit.Pins(),self.cfg)

    def test_public_sequence_changed_same_length_refused(self):
        seqs=self.sequences.copy();seqs[0]='AC'
        with self.assertRaises(ValueError):audit.reconcile_public_design(self.public,self.summary,self.design,seqs)

    def test_public_gene_identity_mutation_refused(self):
        bad=copy.deepcopy(self.public);bad['genes'][0]['id']='gene:OTHER'
        with self.assertRaisesRegex(ValueError,'gene/transcript'):
            audit.reconcile_public_design(bad,self.summary,self.design,self.sequences)

    def test_public_window_membership_mutation_refused(self):
        bad=copy.deepcopy(self.public);bad['windows'][0]['gene_ids']=['gene:g2']
        with self.assertRaisesRegex(ValueError,'window membership'):
            audit.reconcile_public_design(bad,self.summary,self.design,self.sequences)

    def test_public_summary_denominator_mutation_refused(self):
        bad=copy.deepcopy(self.summary);bad['panel_union_bases']+=1
        with self.assertRaisesRegex(ValueError,'denominator'):
            audit.reconcile_public_design(self.public,bad,self.design,self.sequences)

    def prepare_launch(self):
        inner=save_fixture(self.root/audit.CONFIG_REL,self.cfg)
        run=self.root/audit.RUN_REL
        state=save_fixture(run/'status.json',dict(status='complete'))
        summary=save_fixture(run/'summary.json',dict(status='complete'))
        return inner['sha256'],launch_fixture(self.root,self.cfg,state,summary)

    def test_launch_complete_bindings_accepted(self):
        sha,launch=self.prepare_launch()
        with patch.object(audit,'CONFIG_SHA',sha):audit.bind_launch(self.root,audit.Pins(),self.cfg,launch)

    def test_launch_wrong_explicit_hash_refused(self):
        sha,launch=self.prepare_launch()
        with patch.object(audit,'CONFIG_SHA',sha):
            with self.assertRaisesRegex(ValueError,'sha256'):
                audit.bind_launch(self.root,audit.Pins(),self.cfg,'0'*64)

    def test_launch_incomplete_refused(self):
        sha,_=self.prepare_launch();run=self.root/audit.RUN_REL
        statepath=run.with_name(run.name+'-launch')/'status.json';s=json.loads(statepath.read_text());s['status']='running'
        launch=save_fixture(statepath,s)['sha256']
        with patch.object(audit,'CONFIG_SHA',sha):
            with self.assertRaisesRegex(ValueError,'launch incomplete'):
                audit.bind_launch(self.root,audit.Pins(),self.cfg,launch)

    def test_launch_unbound_output_mutation_refused(self):
        sha,launch=self.prepare_launch();(self.root/audit.RUN_REL/'summary.json').write_text('{}')
        with patch.object(audit,'CONFIG_SHA',sha):
            with self.assertRaises(ValueError):audit.bind_launch(self.root,audit.Pins(),self.cfg,launch)

    def mutate_guard_controls(self, mutator):
        sha,_=self.prepare_launch();run=self.root/audit.RUN_REL;launch=run.with_name(run.name+'-launch')
        c=json.loads((launch/'controls.json').read_text());mutator(c)
        pin=save_fixture(launch/'controls.json',c)
        state=json.loads((launch/'status.json').read_text());state['qualification']=pin
        return sha,save_fixture(launch/'status.json',state)['sha256']

    def test_guard_failed_qualification_refused(self):
        sha,digest=self.mutate_guard_controls(lambda q:q.update(failures=1,status='failed'))
        with patch.object(audit,'CONFIG_SHA',sha):
            with self.assertRaisesRegex(ValueError,'guard qualification'):
                audit.bind_launch(self.root,audit.Pins(),self.cfg,digest)

    def test_guard_stale_source_qualification_refused(self):
        def mutate(q):q['source_files']['scripts/launch_secondary_depth_v2.py']['sha256']='0'*64
        sha,digest=self.mutate_guard_controls(mutate)
        with patch.object(audit,'CONFIG_SHA',sha):
            with self.assertRaisesRegex(ValueError,'sha256'):
                audit.bind_launch(self.root,audit.Pins(),self.cfg,digest)

    def test_launch_other_source_same_bytes_wrong_path_refused(self):
        sha,_=self.prepare_launch();run=self.root/audit.RUN_REL;launch=run.with_name(run.name+'-launch')
        other=launch/'other.py';shutil.copyfile(launch/'workflow/scripts/launch_secondary_depth_v2.py',other)
        st=json.loads((launch/'status.json').read_text());st['source']=audit.Pins().add(other)
        digest=save_fixture(launch/'status.json',st)['sha256']
        with patch.object(audit,'CONFIG_SHA',sha):
            with self.assertRaisesRegex(ValueError,'source path differs'):
                audit.bind_launch(self.root,audit.Pins(),self.cfg,digest)


try:
    import pysam
except ImportError:
    pysam = None
SAMTOOLS = os.environ.get('MVA_SAMTOOLS')


@unittest.skipUnless(pysam is not None and SAMTOOLS and Path(SAMTOOLS).is_file(), 'native pysam/samtools required')
class NativeControls(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.ref = self.root / 'ref.fa'
        self.ref.write_text('>1\n' + 'ACGT' * 15 + '\n')
        pysam.faidx(str(self.ref))
        self.header = dict(HD=dict(VN='1.6', SO='coordinate'), SQ=[dict(SN='1', LN=60)], RG=[dict(ID='RG1', SM='FICTIONAL')])
        self.bam = self.root / 'original.bam'
        with pysam.AlignmentFile(str(self.bam), 'wb', header=self.header) as stream:
            for name, start, flag, qualities in [('same', 0, 0, [30]*12), ('same', 0, 0, [30]*12), ('excluded', 1, 1024, None), ('low', 7, 0, [0]*5), ('outside', 40, 0, [30]*5)]:
                r = pysam.AlignedSegment(stream.header)
                r.query_name = name; r.reference_id = 0; r.reference_start = start
                r.mapping_quality = 60; r.flag = flag; r.cigartuples = [(0, 12 if name == 'same' else 5)]
                r.query_sequence = 'A' * r.cigartuples[0][1]; r.query_qualities = qualities
                r.set_tag('RG', 'RG1'); r.set_tag('XI', 300); r.set_tag('XB', array('i', [1, 2]))
                stream.write(r)
        pysam.index(str(self.bam))
        self.design = fixture(self.root)
        self.config = dict(sample='FICTIONAL', expected_read_groups=1)
        self.bed = self.root / 'regions.bed'
        self.bed.write_text(''.join(f"1\t{w['start']}\t{w['end']}\n" for w in self.design['windows']))
        self.selected = self.root / 'selected.bam'
        subprocess.run([SAMTOOLS, 'view', '--no-PG', '-b', '-M', '-L', str(self.bed), '-o', str(self.selected), str(self.bam)], check=True, capture_output=True)

    def test_original_BAM_complete_payload_and_identical_record_multiplicity(self):
        result = audit.compare_extraction(self.bam, self.selected, self.ref, Path(str(self.bam)+'.bai'), self.design, self.config)
        self.assertEqual(result['records'], 4)
        self.assertEqual(result['unique_payloads'], 3)
        self.assertEqual(result['groups'], {'RG1': 4})

    def test_original_CRAM_payload_selection(self):
        cram = self.root / 'original.cram'
        subprocess.run([SAMTOOLS, 'view', '--no-PG', '-T', str(self.ref), '-O', 'CRAM,version=3.0,store_md=1,store_nm=1', '-o', str(cram), str(self.bam)], check=True, capture_output=True)
        subprocess.run([SAMTOOLS, 'index', str(cram)], check=True, capture_output=True)
        selected = self.root / 'cram-selected.bam'
        subprocess.run([SAMTOOLS, 'view', '--no-PG', '-T', str(self.ref), '-b', '-M', '-L', str(self.bed), '-o', str(selected), str(cram)], check=True, capture_output=True)
        result = audit.compare_extraction(cram, selected, self.ref, Path(str(cram)+'.crai'), self.design, self.config)
        self.assertEqual(result['records'], 4)

    def test_selected_record_omission_refused(self):
        wrong = self.root / 'wrong.bam'
        with pysam.AlignmentFile(str(self.selected)) as inp, pysam.AlignmentFile(str(wrong), 'wb', header=inp.header) as out:
            next(inp)
            for read in inp: out.write(read)
        with self.assertRaisesRegex(ValueError, 'payload or multiplicity'):
            audit.compare_extraction(self.bam, wrong, self.ref, Path(str(self.bam)+'.bai'), self.design, self.config)

    def test_selected_tag_mutation_refused(self):
        wrong = self.root / 'wrong.bam'
        with pysam.AlignmentFile(str(self.selected)) as inp, pysam.AlignmentFile(str(wrong), 'wb', header=inp.header) as out:
            for read in inp:
                read.set_tag('XI', 301); out.write(read)
        with self.assertRaisesRegex(ValueError, 'payload or multiplicity'):
            audit.compare_extraction(self.bam, wrong, self.ref, Path(str(self.bam)+'.bai'), self.design, self.config)

    def test_selected_unknown_BQ_refused_before_reporting(self):
        wrong = self.root / 'wrong.bam'
        with pysam.AlignmentFile(str(self.selected)) as inp, pysam.AlignmentFile(str(wrong), 'wb', header=inp.header) as out:
            for read in inp:
                if read.flag == 0: read.query_qualities = None
                out.write(read)
        with self.assertRaisesRegex(ValueError, 'unknown BQ'):
            audit.compare_extraction(self.bam, wrong, self.ref, Path(str(self.bam)+'.bai'), self.design, self.config)

    def test_pinned_same_engine_recount_has_known_depths(self):
        result = subprocess.run([SAMTOOLS, 'depth', '-aa', '-b', str(self.bed), '-q', '20', '-Q', '20', '-G', '0xF04', '-s', str(self.selected)], check=True, capture_output=True, text=True)
        depths = audit.positions_from_lines(result.stdout.splitlines(), self.design['windows'])
        self.assertEqual(set(depths.values()), {2})

    def test_full_fictional_producer_and_independent_audit(self):
        # Only frozen source code is read from the project; all alignment,
        # annotation, configuration and output data are fictional and temporary.
        project = Path(os.environ['MVA_PROJECT_ROOT'])
        frozen = project / 'methods/coding-depth-20260920/qualification-v2'
        method = self.root / 'methods/coding-depth-20260920'
        qualified = method / 'qualification-v2'
        scripts = qualified / 'workflow/scripts'
        scripts.mkdir(parents=True)
        source_names = ['measure_coding_depth.py', 'coding_depth_targets.py', 'survey_chromosome_depth.py', 'prepare_gatk_compatible_cram.py']
        for name in source_names:
            shutil.copyfile(frozen / 'workflow/scripts' / name, scripts / name)
        cfg = dict(schema_version=1, assembly='GRCh38', sample='FICTIONAL', expected_read_groups=1,
            genes=['ONE', 'TWO'], depth_modes=audit.MODES, source_sha256={}, inputs={})
        for name in source_names:
            cfg['source_sha256'][name] = hashlib.sha256((scripts / name).read_bytes()).hexdigest()
        paths = dict(alignment=self.bam, alignment_index=Path(str(self.bam)+'.bai'),
            reference=self.ref, reference_index=Path(str(self.ref)+'.fai'),
            gff=self.root / 'annotation.gff3.gz', samtools=Path(SAMTOOLS))
        for name, path in paths.items():
            cfg['inputs'][name] = dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        config_path = self.root / audit.CONFIG_REL
        config_path.parent.mkdir(parents=True)
        config_path.write_text(json.dumps(cfg))
        config_sha = hashlib.sha256(config_path.read_bytes()).hexdigest()
        controls = dict(status='complete', exit_code=0, scope='Fictional fixture qualification binding only',
            source_sha256={'workflow/scripts/'+n: cfg['source_sha256'][n] for n in source_names})
        controls_path = qualified / 'controls.json'
        controls_path.write_text(json.dumps(controls))
        controls_sha = hashlib.sha256(controls_path.read_bytes()).hexdigest()
        output = self.root / audit.RUN_REL
        command = [sys.executable, '-B', str(scripts / 'measure_coding_depth.py'), '--config', str(config_path),
            '--expected-config-sha256', config_sha, '--output', str(output)]
        subprocess.run(command, check=True, capture_output=True)
        with pysam.FastaFile(str(self.ref)) as fasta:
            sequences=[fasta.fetch(w['chrom'],w['start'],w['end']).upper() for w in self.design['windows']]
        _,_,public_sha=public_fixture(self.root,cfg,self.design,sequences)
        pins=audit.Pins()
        launch_sha=launch_fixture(self.root,cfg,pins.add(output/'status.json'),pins.add(output/'summary.json'))
        own = self.root / 'audit'; own.mkdir()
        with patch.object(audit, 'CONFIG_SHA', config_sha), patch.object(audit, 'QUALIFICATION_SHA', controls_sha), patch.object(audit,'PUBLIC_STATUS_SHA',public_sha), patch.object(audit,'EXPECTED_GENES',2):
            report = audit.audit(self.root, own, SOURCE, launch_sha)
        self.assertEqual(report['reference_only_target_reconstruction']['panel_union_bases'], 8)
        self.assertEqual(report['reference_only_target_reconstruction']['gene_union_bases_sum'], 10)
        self.assertEqual(report['extraction']['records'], 4)
        self.assertEqual(report['depth_modes']['MQ20_BQ20']['groups']['panel_union']['span']['histogram'], {'2': 8})
        self.assertEqual(len(report['same_engine_technical_recounts']), 3)


if __name__ == '__main__':
    unittest.main()
