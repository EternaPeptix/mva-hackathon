#!/usr/bin/env python3
"""Independent fictional controls; no case data or producer imports."""
import argparse, copy, hashlib, importlib.util, json, tempfile, unittest
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('audit',HERE/'verify-overlay.py')
M=importlib.util.module_from_spec(spec);spec.loader.exec_module(M)
REF='ACGTACGTACGTGG';DEL='ACGTACGTGG';JOINT='ATGTACGTGG'

def fixture():
    contexts={}
    def context(caller,ordinal,pos,ref,alts,gt,pl):
        identity=caller+':'+('a' if caller=='callerA' else 'b')*64+':'+str(ordinal)
        c=dict(source_record_id=identity,caller=caller,chrom='1',pos=pos,ref=ref,alts=alts,FORMAT=dict(GT=gt,PL=pl,AD=[7]*len(alts+[ref]),DP=14),GT_phased=False,FILTER=[],QUAL=0.0,rendered_record='FICTIONAL full original record')
        contexts[identity]=c
        return dict(source_record_id=identity,source_alt_index=1)
    roles={}
    roles['raw_deletion']=context('callerB',1,1,'ACGTACGTA',['ACGTA'],[0,1],[9,0,8])
    roles['disjoint_substitution']=context('callerB',2,6,'C',['T'],[0,1],[7,0,6])
    roles['overlapping_substitution']=context('callerB',3,3,'G',['A'],[1,1],[8,4,0])
    roles['comparison_deletions']=[context('callerA',1,1,'ACG',['A'],[0,1],[5,0,6]),context('callerA',2,4,'TAC',['T'],[0,1],[5,0,6])]
    roles['spanning_call']=context('callerA',3,2,'C',['*','G'],[1,2],[9,8,7,6,0,5]);roles['spanning_call']['source_alt_index']=2
    spec=dict(target_id='FICTIONAL',chrom='1',roles=roles,normalized_key='1:1:ACGTA:A',normalized_deletion=dict(chrom='1',pos=1,ref='ACGTA',alt='A'))
    def edit(role):
        c=contexts[role['source_record_id']];idx=role['source_alt_index']
        return dict(chrom=c['chrom'],pos=c['pos'],ref=c['ref'],alt=c['alts'][idx-1],caller=c['caller'],source_record_id=c['source_record_id'],source_alt_index=idx,called=idx in c['FORMAT']['GT'])
    raw,snv,overlap=[edit(roles[x]) for x in ['raw_deletion','disjoint_substitution','overlapping_substitution']]
    comps=list(map(edit,roles['comparison_deletions']));normal=spec['normalized_deletion'];star=edit(roles['spanning_call'])
    common=dict(caller='callerB',normalization_key=spec['normalized_key'],original_deletion_edit=raw,normalized_deletion=normal,raw_normalized_sequence_equal=True,interpretation='Constructed normalized-representation lineage; not native genotype or biological phase.')
    derived=[dict(sequence=DEL,lineage=common|dict(label='normalized_equivalent_deletion',substitution_edits=[])),dict(sequence=JOINT,lineage=common|dict(label='normalized_deletion_plus_disjoint_substitution',substitution_edits=[snv],exact_sequence_comparison_lineage=comps))]
    proof=dict(lineages=derived,rejected_compositions=[dict(label=label,edits=edits,reason='overlapping REF edits') for label,edits in [('raw_deletion_plus_disjoint_substitution',[raw,snv]),('normalized_deletion_plus_overlapping_substitution',[normal,overlap]),('normalized_deletion_plus_both_substitutions',[normal,overlap,snv])]],preserved_spanning_call=contexts[star['source_record_id']],spanning_concrete_ALT_edit=star,raw_normalized_equal=True,joint_sequence_equal=True)
    old=[dict(id='H0',sequence=REF,source_combinations=[dict(caller='callerA',edits=[],reference_only=True)]),dict(id='H1',sequence=DEL,source_combinations=[dict(caller='callerB',edits=[raw],reference_only=False)]),dict(id='H2',sequence=JOINT,source_combinations=[dict(caller='callerA',edits=comps,reference_only=False)])]
    classes={}
    for h in old:
        sequence=h['sequence'];classes[sequence]=dict(sequence=sequence,original_hypothesis_ids=[h['id']],original_lineages=h['source_combinations'],derived_lineages=[d['lineage'] for d in derived if d['sequence']==sequence],id='ADAMTS1_F3_U%04d'%sorted([REF,DEL,JOINT]).index(sequence),sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),sequence_length=len(sequence))
    return spec,contexts,proof,old,classes

def original_arm(sequences=(JOINT,JOINT),names=None):
    names=names or [('RG1','pair')]*len(sequences)
    rows=[dict(observation_id=i+1,RG=identity[0],QNAME=identity[1],status='accepted' if seq is not None else 'rejected',sequence=seq,mapping_quality=60,base_qualities=None if seq is None else [30]*len(seq),original_FLAG=99 if i==0 else 147) for i,(seq,identity) in enumerate(zip(sequences,names))]
    fragments=[];counts=dict(raw_observations=len(rows),fragments_with_identity=0,accepted_observations=sum(s is not None for s in sequences),conflicting_fragments=0,nonconflicting_spanning_fragments=0,matched_fragments=0,unmatched_fragments=0)
    for identity in sorted(set(names)):
        r=[v for v in rows if (v['RG'],v['QNAME'])==identity];a=[v for v in r if v['status']=='accepted'];seqs=sorted({v['sequence'] for v in a})
        status='no_eligible_spanning_observation'
        if len(seqs)>1:status='mate_or_repeated_observation_conflict';counts['conflicting_fragments']+=1
        elif seqs:
            counts['nonconflicting_spanning_fragments']+=1;status='matched' if seqs[0] in [REF,DEL,JOINT] else 'unmatched';counts[status+'_fragments']+=1
        fragments.append(dict(RG=identity[0],QNAME=identity[1],observation_ids=[v['observation_id'] for v in r],accepted_observation_ids=[v['observation_id'] for v in a],sequences=seqs,status=status))
    counts['fragments_with_identity']=len(fragments)
    return dict(observations=rows,fragments=fragments,counts=counts,window=dict(start0=0,end0=14),thresholds=dict(min_MQ=0,min_raw_BQ=20,excluded_flags=0xF04,reject_MQ_255=True,BAQ_applied=False),rejection_counts={},raw_context='FICTIONAL kept verbatim')

def content():
    spec,contexts,proof,old,classes=fixture();arm=original_arm();cfg=dict(window_flanks_bp=[3],mapping_quality_thresholds=[0],base_quality_thresholds=[20])
    pre=dict(chrom='1',start0=0,end0=14,flank_bp=3,reference_sequence=REF,native_family=dict(haplotypes=old,rejected_combinations=[]),omitted_native_alleles=[dict(alt='*')])
    w={k:pre[k] for k in ['chrom','start0','end0','flank_bp','reference_sequence']};w.update(hypothesis_ids=[h['id'] for h in old],enumeration_metadata=dict(rejected_combinations=[]),omitted_overlapping_native_alleles=[dict(alt='*')],arms=[dict(minimum_mapping_quality=0,minimum_base_quality=20,evidence=arm)])
    native=dict(config=cfg,hypotheses={h['id']:h for h in old},windows=[w]);context=dict(targets={'FICTIONAL':[dict(context=c) for c in contexts.values()]});selection=dict(overlay_spec=spec,native_config_body=cfg,contexts=contexts,windows=[pre])
    counts=dict(raw_observations=2,raw_fragments=1,accepted_observations=2,conflicting_fragments=0,nonconflicting_spanning_fragments=1,matched_fragments=1,unmatched_fragments=0,fragments_in_derived_lineage_classes=1)
    support={v['id']:dict(observations=2 if seq==JOINT else 0,fragments=1 if seq==JOINT else 0) for seq,v in classes.items()}
    assignment=[dict(RG='RG1',QNAME='pair',observation_ids=[1,2],accepted_observation_ids=[1,2],class_id=classes[JOINT]['id'],status='matched')]
    a=dict(minimum_mapping_quality=0,minimum_base_quality=20,counts=counts,class_support=support,fragment_assignments=assignment,original_observation_evidence=arm)
    ow=dict(chrom='1',start0=0,end0=14,flank_bp=3,original_class_count=3,union_class_count=3,classes=list(classes.values()),overlay=proof,arms=[a]);overlay=dict(native_contexts=contexts,windows=[ow]);summary=dict(windows=[dict(flank_bp=3,original_class_count=3,union_class_count=3,arms=[dict(MQ=0,raw_BQ=20,counts=counts,class_support=support)])])
    return tuple(json.loads(json.dumps(v)) for v in (native,context,selection,overlay,summary))

class Controls(unittest.TestCase):
    def test_01_literal_raw_normalization_and_joint_equality(self):
        s,c,p,old,classes=fixture();self.assertEqual(M.overlay_proof(REF,0,s,c),p);self.assertEqual(p['lineages'][0]['sequence'],DEL);self.assertEqual(p['lineages'][1]['sequence'],JOINT)
    def test_02_full_literal_fixture_actual_config_key(self):
        result=M.audit_content(*content(),expected_classes=3);self.assertEqual(result['arms'],1);self.assertEqual(result['new_sequence_classes'],0)
    def test_03_wrong_REF_refused(self):
        with self.assertRaises(ValueError):M.compose(REF,0,[dict(pos=1,ref='T',alt='A')])
    def test_04_boundary_and_original_overlap_refused(self):
        for edits in [[dict(pos=14,ref='GG',alt='G')],[dict(pos=1,ref='ACG',alt='A'),dict(pos=3,ref='G',alt='A')]]:
            with self.assertRaises(ValueError):M.compose(REF,0,edits)
    def test_05_star_GT_PL_not_rewritten(self):
        s,c,p,_,_=fixture();before=copy.deepcopy(c);M.overlay_proof(REF,0,s,c);self.assertEqual(c,before)
        for gt,pl in [([0,1],[9,8,7,6,0,5]),([1,2],[9,0,5])]:
            changed=copy.deepcopy(c);changed[s['roles']['spanning_call']['source_record_id']]['FORMAT'].update(GT=gt,PL=pl)
            with self.assertRaises(ValueError):M.overlay_proof(REF,0,s,changed)
    def test_06_wrong_normalized_sequence_refused(self):
        s,c,_,_,_=fixture();s['normalized_deletion']['ref']='ACG';s['normalized_key']='1:1:ACG:A'
        with self.assertRaises(ValueError):M.overlay_proof(REF,0,s,c)
    def test_07_wrong_original_ALT_index_refused(self):
        s,c,_,_,_=fixture();s['roles']['raw_deletion']['source_alt_index']=2
        with self.assertRaises(ValueError):M.overlay_proof(REF,0,s,c)
    def test_08_original_source_payload_mutation_refused(self):
        values=content();values[3]['native_contexts'][next(iter(values[3]['native_contexts']))]['FORMAT']['PL'][0]=100
        with self.assertRaises(ValueError):M.audit_content(*values,expected_classes=3)
    def test_09_original_class_omission_refused(self):
        values=content();values[3]['windows'][0]['classes'].pop()
        with self.assertRaises(ValueError):M.audit_content(*values,expected_classes=3)
    def test_10_original_lineage_mutation_refused(self):
        values=content();values[3]['windows'][0]['classes'][1]['original_lineages'][0]['caller']='wrong'
        with self.assertRaises(ValueError):M.audit_content(*values,expected_classes=3)
    def test_11_derived_lineage_mutation_refused(self):
        values=content();values[3]['windows'][0]['classes'][1]['derived_lineages'][0]['original_deletion_edit']['source_alt_index']=9
        with self.assertRaises(ValueError):M.audit_content(*values,expected_classes=3)
    def test_12_repeated_mates_count_one_fragment(self):
        *_,classes=fixture();result=M.count_arm(original_arm(),classes);self.assertEqual(result['counts']['matched_fragments'],1);self.assertEqual(result['counts']['accepted_observations'],2)
    def test_13_RG_distinguishes_same_QNAME(self):
        *_,classes=fixture();r=M.count_arm(original_arm(names=[('RG1','pair'),('RG2','pair')]),classes);self.assertEqual(r['counts']['matched_fragments'],2)
    def test_14_conflicting_mates_unassigned(self):
        *_,classes=fixture();r=M.count_arm(original_arm([JOINT,DEL]),classes);self.assertEqual(r['counts']['conflicting_fragments'],1);self.assertEqual(r['counts']['matched_fragments'],0);self.assertEqual(sum(v['observations'] for v in r['class_support'].values()),2)
    def test_15_unmatched_and_rejected_are_not_reference(self):
        *_,classes=fixture()
        for seq,status in [('CCCC','unmatched'),(None,'no_eligible_span')]:
            r=M.count_arm(original_arm([seq]),classes);self.assertEqual(r['fragment_assignments'][0]['status'],status);self.assertEqual(r['counts']['matched_fragments'],0)
    def test_16_duplicated_observation_refused(self):
        *_,classes=fixture();a=original_arm();a['observations'][1]['observation_id']=1
        with self.assertRaises(ValueError):M.count_arm(a,classes)
    def test_17_missing_fragment_membership_refused(self):
        *_,classes=fixture();a=original_arm();a['fragments'][0]['accepted_observation_ids']=[1]
        with self.assertRaises(ValueError):M.count_arm(a,classes)
    def test_18_observations_preserved_verbatim(self):
        values=content();values[3]['windows'][0]['arms'][0]['original_observation_evidence']['observations'][0]['mapping_quality']=0
        with self.assertRaises(ValueError):M.audit_content(*values,expected_classes=3)
    def test_19_summary_counter_mutation_refused(self):
        values=content();values[4]['windows'][0]['arms'][0]['counts']['matched_fragments']=2
        with self.assertRaises(ValueError):M.audit_content(*values,expected_classes=3)
    def test_20_missing_arm_and_window_refused(self):
        for what in ['arms','windows']:
            values=content()
            if what=='arms':values[3]['windows'][0]['arms']=[]
            else:values[3]['windows']=[]
            with self.assertRaises(ValueError):M.audit_content(*values,expected_classes=3)
    def test_21_json_types_duplicatekeys_and_nan(self):
        for text in ['{"x":1,"x":2}','{"x":NaN}']:
            with self.assertRaises(ValueError):M.parse(text)
        with self.assertRaises(ValueError):M.same(True,1,'bool is not int')
    def test_22_fresh_hash_and_read_access_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'source.json';p.write_text('{}');pins=M.Pins();pin=pins.add(p)
            with self.assertRaises(ValueError):M.Pins().add(p,dict(sha256='0'*64))
            p.write_text('{ }')
            with self.assertRaises(ValueError):pins.stable()
            with self.assertRaises(ValueError):M.Pins().add(Path(d)/'forbidden.cram')
    def test_23_new_sequence_and_wrong_class_count_refused(self):
        s,c,p,old,classes=fixture()
        with self.assertRaises(ValueError):M.sequence_classes(old,[dict(sequence='CCCC',lineage={})],'x',REF,0,3)
        with self.assertRaises(ValueError):M.sequence_classes(old,[],'x',REF,0,9)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise ValueError('preserve existing qualification')
    r=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Controls))
    receipt=dict(status='complete' if r.wasSuccessful() else 'failed',completed_at=M.stamp(),tests_run=r.testsRun,errors=len(r.errors),failures=len(r.failures),skipped=len(r.skipped),source_sha256={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ['verify-overlay.py','test_verify_overlay.py']},scope='Fictional literal sequence and fragment controls; no producer imports or case reads.')
    a.output.parent.mkdir(parents=True,exist_ok=True);M.write(a.output,receipt)
    raise SystemExit(0 if r.wasSuccessful() else 1)
