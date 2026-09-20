#!/usr/bin/env python3
"""Fictional overlay controls, including a bounded completed-input fixture."""
import argparse, copy, hashlib, importlib.util, json, os, tempfile, unittest
from pathlib import Path
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('overlay',HERE/'normalized-overlay.py')
M=importlib.util.module_from_spec(spec);spec.loader.exec_module(M)

def fixture():
    reference='ACGTACGTACGTGG';contexts={}
    def context(caller,ordinal,pos,ref,alts,gt,pl):
        identity=caller+':'+('a' if caller=='callerA' else 'b')*64+':'+str(ordinal)
        row={'source_record_id':identity,'caller':caller,'chrom':'1','pos':pos,'ref':ref,'alts':alts,
            'FORMAT':{'GT':gt,'PL':pl,'AD':[7]*(len(alts)+1),'DP':14},'GT_phased':False,'FILTER':[],
            'QUAL':0.0,'rendered_record':'FICTIONAL native record with full payload'}
        contexts[identity]=row;return {'source_record_id':identity,'source_alt_index':1}
    raw=context('callerB',1,1,'ACGTACGTA',['ACGTA'],[0,1],[9,0,8])
    disjoint=context('callerB',2,6,'C',['T'],[0,1],[7,0,6])
    overlap=context('callerB',3,3,'G',['A'],[1,1],[8,4,0])
    first=context('callerA',1,1,'ACG',['A'],[0,1],[5,0,6])
    second=context('callerA',2,4,'TAC',['T'],[0,1],[5,0,6])
    star=context('callerA',3,2,'C',['*','G'],[1,2],[9,8,7,6,0,5]);star['source_alt_index']=2
    selection={'target_id':'FICTIONAL','chrom':'1','normalized_deletion':{'chrom':'1','pos':1,'ref':'ACGTA','alt':'A'},
        'normalized_key':'1:1:ACGTA:A','roles':{'raw_deletion':raw,'disjoint_substitution':disjoint,'overlapping_substitution':overlap,
            'comparison_deletions':[first,second],'spanning_call':star}}
    return reference,selection,contexts

def old_classes(reference,selection,contexts):
    a=M.overlay_lineages(reference,0,selection,contexts)
    raw=M.native_edit(contexts,selection['roles']['raw_deletion'])
    comparison=[M.native_edit(contexts,r) for r in selection['roles']['comparison_deletions']]
    return [dict(id='H0',sequence=reference,source_combinations=[dict(caller='callerA',edits=[],reference_only=True)]),
        dict(id='H1',sequence=a['lineages'][0]['sequence'],source_combinations=[dict(caller='callerB',edits=[raw],reference_only=False)]),
        dict(id='H2',sequence=a['lineages'][1]['sequence'],source_combinations=[dict(caller='callerA',edits=comparison,reference_only=False)])]

def arm(sequences,known):
    obs=[{'observation_id':i+1,'RG':'FICTIONAL_RG','QNAME':'FICTIONAL_pair','status':'accepted' if seq is not None else 'rejected','sequence':seq,
        'strand':'forward' if i==0 else 'reverse','mapping_quality':60,'base_qualities':[30]*len(seq or '')} for i,seq in enumerate(sequences)]
    accepted=[o for o in obs if o['status']=='accepted'];unique=sorted({o['sequence'] for o in accepted})
    conflict=len(unique)>1;nonconflicting=len(unique)==1;matched=nonconflicting and unique[0] in known
    status='mate_or_repeated_observation_conflict' if conflict else ('matched' if matched else 'unmatched') if nonconflicting else 'no_eligible_spanning_observation'
    fragment={'RG':'FICTIONAL_RG','QNAME':'FICTIONAL_pair','observation_ids':[o['observation_id'] for o in obs],
        'accepted_observation_ids':[o['observation_id'] for o in accepted],'sequences':unique,'status':status}
    return {'observations':obs,'fragments':[fragment],
        'counts':{'raw_observations':len(obs),'fragments_with_identity':1,'accepted_observations':len(accepted),'conflicting_fragments':int(conflict),
            'nonconflicting_spanning_fragments':int(nonconflicting),'matched_fragments':int(matched),'unmatched_fragments':int(nonconflicting and not matched)},
        'rejection_counts':{},'read_group_counts':{'FICTIONAL_RG':{'raw_observations':len(obs)}},'accepted_observation_strand_counts':{},
        'window':{'start0':0,'end0':14},'thresholds':{'min_MQ':0,'min_raw_BQ':20,'excluded_flags':0xF04,'reject_MQ_255':True,'BAQ_applied':False}}

class Controls(unittest.TestCase):
    def setUp(self):
        self.reference,self.selection,self.contexts=fixture();self.old=old_classes(self.reference,self.selection,self.contexts)
        self.overlay=M.overlay_lineages(self.reference,0,self.selection,self.contexts)
        self.classes=M.union_classes(self.old,self.overlay['lineages'],'FICTIONAL')
    def test_long_REF_normalized_literal_equality(self):
        self.assertEqual(self.overlay['lineages'][0]['sequence'],'ACGTACGTGG')
        self.assertTrue(self.overlay['raw_normalized_equal'])
    def test_disjoint_substitution_equals_two_other_deletions(self):
        self.assertEqual(self.overlay['lineages'][1]['sequence'],'ATGTACGTGG')
        self.assertTrue(self.overlay['joint_sequence_equal'])
    def test_all_three_overlap_refusals_retained(self):
        self.assertEqual(len(self.overlay['rejected_compositions']),3)
        self.assertTrue(all(r['reason']=='overlapping REF edits' for r in self.overlay['rejected_compositions']))
    def test_full_star_GT_ALT2_six_PL_context_unchanged(self):
        before=copy.deepcopy(self.contexts);out=M.overlay_lineages(self.reference,0,self.selection,self.contexts)
        self.assertEqual(out['preserved_spanning_call'],before[self.selection['roles']['spanning_call']['source_record_id']])
        self.assertEqual(out['preserved_spanning_call']['FORMAT']['GT'],[1,2]);self.assertEqual(out['preserved_spanning_call']['FORMAT']['PL'],[9,8,7,6,0,5])
        self.assertEqual(out['spanning_concrete_ALT_edit']['source_alt_index'],2);self.assertTrue(out['spanning_concrete_ALT_edit']['called']);self.assertEqual(before,self.contexts)
    def test_STAR_index_or_PL_recode_refused(self):
        for kind in ('index','GT','PL'):
            s,c=copy.deepcopy(self.selection),copy.deepcopy(self.contexts);ctx=c[s['roles']['spanning_call']['source_record_id']]
            if kind=='index':s['roles']['spanning_call']['source_alt_index']=1
            elif kind=='GT':ctx['FORMAT']['GT']=[0,1]
            else:ctx['FORMAT']['PL']=[9,0,5]
            with self.subTest(kind=kind),self.assertRaises(ValueError):M.overlay_lineages(self.reference,0,s,c)
    def test_equivalent_classes_preserve_both_lineages(self):
        self.assertEqual(len(self.classes),len(self.old));joint=self.classes['ATGTACGTGG']
        self.assertEqual(joint['original_lineages'][0]['caller'],'callerA');self.assertEqual(joint['derived_lineages'][0]['caller'],'callerB')
    def test_two_observations_one_fragment_no_double_count(self):
        a=arm(['ATGTACGTGG']*2,self.classes);r=M.recount(a,self.classes)
        self.assertEqual(r['counts']['matched_fragments'],1);self.assertEqual(r['counts']['accepted_observations'],2)
        support=r['class_support'][self.classes['ATGTACGTGG']['id']];self.assertEqual(support,{'observations':2,'fragments':1})
        self.assertEqual(r['original_observation_evidence'],a)
    def test_conflicting_mates_remain_unassigned(self):
        r=M.recount(arm(['ATGTACGTGG','ACGTACGTGG'],self.classes),self.classes)
        self.assertEqual(r['counts']['conflicting_fragments'],1);self.assertEqual(r['counts']['matched_fragments'],0)
    def test_unmatched_and_no_span_not_reference(self):
        self.assertEqual(M.recount(arm(['CCCC'],self.classes),self.classes)['counts']['unmatched_fragments'],1)
        r=M.recount(arm([None,None],self.classes),self.classes);self.assertEqual(r['fragment_assignments'][0]['status'],'no_eligible_span')
        self.assertEqual(r['counts']['accepted_observations'],0)
    def test_missing_fragment_member_refused(self):
        a=arm(['ATGTACGTGG']*2,self.classes);a['fragments'][0]['observation_ids']=[1]
        with self.assertRaises(ValueError):M.recount(a,self.classes)
    def test_duplicate_observation_refused(self):
        a=arm(['ATGTACGTGG']*2,self.classes);a['observations'][1]['observation_id']=1
        with self.assertRaises(ValueError):M.recount(a,self.classes)
    def test_source_identity_missing_or_ALT_index_refused(self):
        for patch in ({'source_record_id':'absent'},{'source_alt_index':9}):
            with self.assertRaises(ValueError):M.native_edit(self.contexts,self.selection['roles']['raw_deletion']|patch)
    def test_comparison_caller_mixing_refused(self):
        c=copy.deepcopy(self.contexts);c[self.selection['roles']['comparison_deletions'][1]['source_record_id']]['caller']='callerB'
        with self.assertRaises(ValueError):M.overlay_lineages(self.reference,0,self.selection,c)
    def test_REF_mismatch_and_window_crossing_refused(self):
        for edits in ([{'pos':1,'ref':'T','alt':'A'}],[{'pos':14,'ref':'GG','alt':'G'}]):
            with self.assertRaises(ValueError):M.apply_edits(self.reference,0,edits)
    def test_wrong_normalization_link_refused(self):
        s=copy.deepcopy(self.selection);s['normalized_deletion']['ref']='ACG';s['normalized_key']='1:1:ACG:A'
        with self.assertRaises(ValueError):M.overlay_lineages(self.reference,0,s,self.contexts)
    def test_new_sequence_silently_added_refused(self):
        with self.assertRaises(ValueError):M.union_classes(self.old,[{'sequence':'CCCC','lineage':{}}],'bad')
    def test_source_qualification_refused(self):
        q={'status':'complete','tests_run':20,'skipped':0,'failures':0,'errors':0,'source_sha256':{'normalized-overlay.py':'a'*64}}
        M.qualified(q,'a'*64)
        for patch in ({'tests_run':0},{'skipped':1},{'errors':1},{'failures':1},{'source_sha256':{'normalized-overlay.py':'b'*64}}):
            with self.assertRaises(ValueError):M.qualified(q|patch,'a'*64)
    def test_output_scope_and_existing_refused(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);native=root/'results/adamts1-joint-review-20260920/native-v1';native.mkdir(parents=True)
            M.new_output(root,native,native.parent/'overlay-v1')
            for out in (native,native/'nested',root/'outside',native.parent):
                with self.assertRaises(ValueError):M.new_output(root,native,out)
    def test_exact_native_pysam_STAR_GT_PL(self):
        import pysam
        scripts=Path(os.environ['MVA_QUALIFIED_SCRIPTS']);s=importlib.util.spec_from_file_location('qualified_records',scripts/'native_coding_records.py')
        records=importlib.util.module_from_spec(s);s.loader.exec_module(records)
        h=pysam.VariantHeader();h.contigs.add('1',length=100);h.add_sample('FICTIONAL')
        h.formats.add('GT',1,'String','Genotype');h.formats.add('PL','G','Integer','Likelihoods')
        r=h.new_record(contig='1',start=1,alleles=('C','*','G'));r.samples['FICTIONAL']['GT']=(1,2);r.samples['FICTIONAL']['PL']=(9,8,7,6,0,5)
        c=records.source_context(r,'FICTIONAL','callerA','a'*64,3)
        self.assertEqual(c['alts'],['*','G']);self.assertEqual(c['FORMAT']['GT'],[1,2]);self.assertEqual(c['FORMAT']['PL'],[9,8,7,6,0,5])
        contexts=copy.deepcopy(self.contexts);contexts[c['source_record_id']]=c
        self.assertEqual(M.overlay_lineages(self.reference,0,self.selection,contexts)['preserved_spanning_call'],c)
    def completed_fixture(self,root):
        native=root/'results/adamts1-joint-review-20260920/native-v1';native.mkdir(parents=True)
        config={'mapping_quality_thresholds':[0],'base_quality_thresholds':[20]};cfg=root/'fictional-config.json';M.save(cfg,config)
        source=M.fingerprint(HERE/'normalized-overlay.py');ctl=root/'controls.json'
        M.save(ctl,{'status':'complete','tests_run':20,'skipped':0,'failures':0,'errors':0,'source_sha256':{'normalized-overlay.py':source['sha256']}})
        pre={'chrom':'1','start0':0,'end0':14,'flank_bp':3,'reference_sequence':self.reference,
            'native_family':{'haplotypes':[{'sequence':r['sequence'],'source_combinations':r['source_combinations']} for r in self.old],'rejected_combinations':[]},'omitted_native_alleles':[{'alt':'*'}]}
        selection=root/'selection.json';M.save(selection,{'status':'complete','mode':'reference_only_preparation','verified_files':{},'native_config_body':config,
            'native_config_pin':M.fingerprint(cfg),'overlay_spec':self.selection,'contexts':self.contexts,'windows':[pre]})
        w={k:pre[k] for k in ('chrom','start0','end0','flank_bp','reference_sequence')};w.update(hypothesis_ids=[r['id'] for r in self.old],
            enumeration_metadata={'rejected_combinations':[]},omitted_overlapping_native_alleles=[{'alt':'*'}],
            arms=[{'minimum_mapping_quality':0,'minimum_base_quality':20,'evidence':arm(['ATGTACGTGG']*2,self.classes)}])
        M.save(native/'evidence.json',{'config':config,'hypotheses':{r['id']:r for r in self.old},'windows':[w]})
        M.save(native/'regional-native-context.json',{'targets':{'FICTIONAL':[{'context':c} for c in self.contexts.values()]}})
        M.save(native/'summary.json',{'status':'complete','inputs':{'config':M.fingerprint(cfg)},'outputs':{n:M.fingerprint(native/n) for n in ('evidence.json','regional-native-context.json')}})
        M.save(native/'status.json',{'status':'complete','summary_sha256':M.fingerprint(native/'summary.json')['sha256']})
        verification=root/'verification.json';M.save(verification,{'status':'complete','verified_files':{n:M.fingerprint(native/n) for n in ('summary.json','status.json','evidence.json','regional-native-context.json')}})
        return SimpleNamespace(root=root,native_run=native,output=native.parent/'overlay-v1',selection=selection,expected_selection_sha256=M.fingerprint(selection)['sha256'],
            controls=ctl,expected_controls_sha256=M.fingerprint(ctl)['sha256'],native_verification=verification,expected_native_verification_sha256=M.fingerprint(verification)['sha256'])
    def test_fictional_end_to_end_no_new_classes_or_fragments(self):
        with tempfile.TemporaryDirectory() as d:
            args=self.completed_fixture(Path(d))
            self.assertNotIn('status',json.loads((args.native_run/'evidence.json').read_text()))
            M.run(args);out=json.loads((args.output/'evidence.json').read_text())
            self.assertEqual(out['status'],'complete');self.assertEqual(out['native_contexts'],self.contexts)
            self.assertEqual(out['windows'][0]['original_class_count'],out['windows'][0]['union_class_count'])
            self.assertEqual(out['windows'][0]['arms'][0]['counts']['matched_fragments'],1)
    def test_stale_independent_input_pin_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            args=self.completed_fixture(Path(d));p=args.native_run/'evidence.json';p.write_text(p.read_text()+' ')
            with self.assertRaises(ValueError):M.run(args)
            self.assertEqual(json.loads((args.output/'status.json').read_text())['status'],'failed')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():raise ValueError('new qualification receipt required')
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Controls))
    names=('normalized-overlay.py','reference-preflight.py','test_normalized_overlay.py','adamts1-native-config.json','overlay-spec.json')
    receipt={'status':'complete' if result.wasSuccessful() else 'failed','completed_at':M.stamp(),'tests_run':result.testsRun,'skipped':len(result.skipped),
        'failures':len(result.failures),'errors':len(result.errors),'source_sha256':{n:M.fingerprint(HERE/n)['sha256'] for n in names},
        'scope':'Fictional controls only; no original alignment read, case extraction or native case run.'}
    args.output.parent.mkdir(parents=True,exist_ok=True);M.save(args.output,receipt)
    raise SystemExit(0 if result.wasSuccessful() else 1)
