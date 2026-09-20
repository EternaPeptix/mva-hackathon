#!/usr/bin/env python3
"""Predeclared normalized lineages over pinned, completed native observations.

No read extraction, alignment, genotype caller, observer, or allele annotation is
run here. Original observations and native context remain unchanged evidence.
"""
import argparse, copy, hashlib, json, re, time
from collections import Counter, defaultdict
from pathlib import Path

def require(ok, message):
    if not ok: raise ValueError(message)
def canonical(value): return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)
def same(a, b, message): require(canonical(a)==canonical(b), message)
def stamp(): return time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
def save(path, value): path.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')
def seq_sha(sequence): return hashlib.sha256(sequence.encode()).hexdigest()

def fingerprint(path, expected=None):
    path=Path(path).resolve();before=path.stat();digest=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''): digest.update(block)
    after=path.stat();same([before.st_size,before.st_mtime_ns],[after.st_size,after.st_mtime_ns],'input changed while hashing')
    pin={'path':str(path),'bytes':after.st_size,'mtime_ns':after.st_mtime_ns,'sha256':digest.hexdigest()}
    if expected:
        for field in ('bytes','sha256'): same(pin[field],expected[field],'input '+field+' differs')
    return pin

def apply_edits(reference,start0,edits):
    require(type(start0) is int and start0>=0 and re.fullmatch('[ACGT]+',reference or '') is not None,'invalid reference window')
    ordered=sorted(edits,key=lambda row:row['pos']);cursor=0;pieces=[]
    for edit in ordered:
        require(type(edit['pos']) is int and edit['pos']>0,'invalid edit position')
        require(all(re.fullmatch('[ACGT]+',edit[k] or '') is not None for k in ('ref','alt')),'nonconcrete edit')
        offset=edit['pos']-1-start0
        require(0<=offset and offset+len(edit['ref'])<=len(reference),'edit crosses window')
        require(reference[offset:offset+len(edit['ref'])]==edit['ref'],'exact REF mismatch')
        require(offset>=cursor,'overlapping REF edits')
        pieces.extend((reference[cursor:offset],edit['alt']));cursor=offset+len(edit['ref'])
    pieces.append(reference[cursor:]);return ''.join(pieces)

def native_edit(contexts,role):
    identity=role['source_record_id'];index=role['source_alt_index']
    require(identity in contexts,'required native source identity missing')
    context=contexts[identity]
    require(type(index) is int and 1<=index<=len(context['alts']),'original ALT index invalid')
    require(context['source_record_id']==identity,'native identity mismatch')
    return {'chrom':context['chrom'],'pos':context['pos'],'ref':context['ref'],'alt':context['alts'][index-1],
        'caller':context['caller'],'source_record_id':identity,'source_alt_index':index,
        'called':index in (context['FORMAT'].get('GT') or [])}

def overlay_lineages(reference,start0,spec,contexts):
    before=canonical(contexts);roles=spec['roles']
    raw=native_edit(contexts,roles['raw_deletion']);disjoint=native_edit(contexts,roles['disjoint_substitution'])
    overlap=native_edit(contexts,roles['overlapping_substitution']);comparison=[native_edit(contexts,r) for r in roles['comparison_deletions']]
    star=native_edit(contexts,roles['spanning_call']);star_context=contexts[star['source_record_id']]
    require(star['source_alt_index']==2 and star_context['alts']==['*',star['alt']], 'spanning-star ALT order differs')
    require(re.fullmatch('[ACGT]+',star['alt'] or '') is not None,'spanning concrete ALT is nonconcrete')
    require(star_context['FORMAT'].get('GT')==[1,2] and len(star_context['FORMAT'].get('PL') or [])==6,'spanning-star GT/PL shape differs')
    require(raw['caller']==disjoint['caller']==overlap['caller'],'derived source callers differ')
    require(len(comparison)==2 and len({e['source_record_id'] for e in comparison})==2 and len({e['caller'] for e in comparison})==1,'comparison source lineage differs')
    require(len({raw['source_record_id'],disjoint['source_record_id'],overlap['source_record_id']})==3,'derived native record identities collide')
    require(all(e['chrom']==spec['chrom'] for e in [raw,disjoint,overlap,star]+comparison),'context chromosome differs')
    normalized=copy.deepcopy(spec['normalized_deletion'])
    require(normalized['chrom']==spec['chrom'],'normalized chromosome differs')
    key=':'.join(map(str,(normalized['chrom'],normalized['pos'],normalized['ref'],normalized['alt'])))
    require(key==spec['normalized_key'],'normalized allele key differs')
    raw_sequence=apply_edits(reference,start0,[raw]);normalized_sequence=apply_edits(reference,start0,[normalized])
    same(raw_sequence,normalized_sequence,'raw/normalized full sequence differs')
    joint=apply_edits(reference,start0,[normalized,disjoint]);comparison_sequence=apply_edits(reference,start0,comparison)
    same(joint,comparison_sequence,'joint full sequence equality differs')
    rejected=[]
    for label,edits in [('raw_deletion_plus_disjoint_substitution',[raw,disjoint]),
                        ('normalized_deletion_plus_overlapping_substitution',[normalized,overlap]),
                        ('normalized_deletion_plus_both_substitutions',[normalized,overlap,disjoint])]:
        try: apply_edits(reference,start0,edits)
        except ValueError as error:
            require(str(error)=='overlapping REF edits','unexpected overlap refusal')
            rejected.append({'label':label,'edits':copy.deepcopy(edits),'reason':str(error)})
        else: raise ValueError('expected overlapping REF refusal missing')
    base={'caller':raw['caller'],'normalization_key':key,'original_deletion_edit':raw,
        'normalized_deletion':normalized,'raw_normalized_sequence_equal':True,
        'interpretation':'Constructed normalized-representation lineage; not native genotype or biological phase.'}
    lineages=[{'sequence':normalized_sequence,'lineage':dict(copy.deepcopy(base),label='normalized_equivalent_deletion',substitution_edits=[])},
        {'sequence':joint,'lineage':dict(copy.deepcopy(base),label='normalized_deletion_plus_disjoint_substitution',substitution_edits=[disjoint],
            exact_sequence_comparison_lineage=comparison)}]
    same(canonical(contexts),before,'native context mutated')
    return {'lineages':lineages,'rejected_compositions':rejected,'preserved_spanning_call':copy.deepcopy(star_context),
        'spanning_concrete_ALT_edit':star,'raw_normalized_equal':True,'joint_sequence_equal':True}

def union_classes(original,lineages,prefix):
    classes={};ids=set()
    for row in original:
        require(row['id'] not in ids,'duplicate original hypothesis id');ids.add(row['id'])
        value=classes.setdefault(row['sequence'],{'sequence':row['sequence'],'original_hypothesis_ids':[],'original_lineages':[],'derived_lineages':[]})
        value['original_hypothesis_ids'].append(row['id']);value['original_lineages'].extend(copy.deepcopy(row['source_combinations']))
    for row in lineages:
        require(row['sequence'] in classes,'expected derived equivalence absent from native class set')
        classes[row['sequence']]['derived_lineages'].append(copy.deepcopy(row['lineage']))
    for index,sequence in enumerate(sorted(classes)):
        classes[sequence].update(id=prefix+'_U%04d'%index,sequence_sha256=seq_sha(sequence),sequence_length=len(sequence))
    require(len(classes)==len(original),'original panel contains duplicate sequence classes')
    return classes

def recount(arm,classes):
    observations=arm['observations'];fragments=arm['fragments'];by_id={o['observation_id']:o for o in observations}
    require(len(by_id)==len(observations),'duplicate observation id')
    expected=defaultdict(list)
    for observation in observations:
        require(observation['status'] in ('accepted','rejected'),'invalid observation status')
        identity=(observation['RG'],observation['QNAME'])
        if all(isinstance(x,str) and x for x in identity): expected[identity].append(observation['observation_id'])
    require(len(fragments)==len(expected),'fragment identity count differs')
    counts=Counter(raw_observations=len(observations),raw_fragments=len(fragments));obs_counts=Counter();frag_counts=Counter();assignments=[];seen=set()
    for observation in observations:
        if observation['status']=='accepted':
            counts['accepted_observations']+=1
            if observation['sequence'] in classes: obs_counts[classes[observation['sequence']]['id']]+=1
    for fragment in fragments:
        identity=(fragment['RG'],fragment['QNAME']);require(identity not in seen,'duplicate fragment');seen.add(identity)
        same(fragment['observation_ids'],expected.get(identity),'fragment record membership/order differs')
        accepted=[by_id[i] for i in fragment['observation_ids'] if by_id[i]['status']=='accepted']
        same([o['observation_id'] for o in accepted],fragment['accepted_observation_ids'],'accepted fragment membership differs')
        sequences=sorted({o['sequence'] for o in accepted});same(sequences,fragment['sequences'],'fragment sequences differ')
        assignment={'RG':identity[0],'QNAME':identity[1],'observation_ids':copy.deepcopy(fragment['observation_ids']),
            'accepted_observation_ids':copy.deepcopy(fragment['accepted_observation_ids']),'class_id':None,'status':'no_eligible_span'}
        if len(sequences)>1:
            require(fragment['status']=='mate_or_repeated_observation_conflict','conflict status differs')
            counts['conflicting_fragments']+=1;assignment['status']='conflict'
        elif sequences:
            counts['nonconflicting_spanning_fragments']+=1;sequence=sequences[0]
            require(fragment['status']==('matched' if sequence in classes else 'unmatched'),'original class assignment differs')
            if sequence in classes:
                value=classes[sequence];counts['matched_fragments']+=1;frag_counts[value['id']]+=1
                assignment.update(status='matched',class_id=value['id'])
                if value['derived_lineages']: counts['fragments_in_derived_lineage_classes']+=1
            else: counts['unmatched_fragments']+=1;assignment['status']='unmatched'
        else: require(fragment['status']=='no_eligible_spanning_observation','no-span status differs')
        assignments.append(assignment)
    for field in ('accepted_observations','conflicting_fragments','nonconflicting_spanning_fragments','matched_fragments','unmatched_fragments'):
        counts.setdefault(field,0);same(counts[field],arm['counts'][field],'original '+field+' differs')
    same(counts['raw_observations'],arm['counts']['raw_observations'],'raw observation denominator differs')
    same(counts['raw_fragments'],arm['counts']['fragments_with_identity'],'raw fragment denominator differs')
    counts.setdefault('fragments_in_derived_lineage_classes',0)
    return {'counts':dict(counts),'class_support':{v['id']:{'observations':obs_counts[v['id']],'fragments':frag_counts[v['id']]} for v in classes.values()},
        'fragment_assignments':assignments,'original_observation_evidence':copy.deepcopy(arm)}

def new_output(root,native_run,output):
    root,native_run,output=map(lambda p:Path(p).resolve(),(root,native_run,output))
    allowed=root/'results/adamts1-joint-review-20260920'
    require(output.is_relative_to(allowed) and output!=allowed,'output outside new owned result directory')
    require(not output.exists() and not output.is_relative_to(native_run) and not native_run.is_relative_to(output),'new separate overlay output required')
    return output

def qualified(controls,source_sha):
    require(controls['status']=='complete' and type(controls['tests_run']) is int and controls['tests_run']>=12,'overlay controls incomplete')
    require(all(type(controls[k]) is int and controls[k]==0 for k in ('skipped','failures','errors')),'overlay controls not wholly passing')
    require(controls['source_sha256']['normalized-overlay.py']==source_sha,'qualified overlay source differs')

def run(args):
    root=args.root.resolve();native=args.native_run.resolve();output=new_output(root,native,args.output)
    output.mkdir(parents=True);state={'status':'running','started_at':stamp()};save(output/'status.json',state);pins={}
    def load(path,expected=None):
        path=Path(path).resolve();pins[str(path)]=fingerprint(path,expected);return json.loads(path.read_text())
    try:
        source=fingerprint(__file__);pins[source['path']]=source
        controls=load(args.controls);require(pins[str(args.controls.resolve())]['sha256']==args.expected_controls_sha256,'controls receipt pin differs');qualified(controls,source['sha256'])
        selection=load(args.selection);require(pins[str(args.selection.resolve())]['sha256']==args.expected_selection_sha256,'selection pin differs')
        require(selection['status']=='complete' and selection['mode']=='reference_only_preparation','preparation not complete')
        for pin in selection['verified_files'].values(): pins[pin['path']]=fingerprint(pin['path'],pin)
        config=selection['native_config_body'];spec=selection['overlay_spec'];contexts=selection['contexts']
        verification=load(args.native_verification);require(pins[str(args.native_verification.resolve())]['sha256']==args.expected_native_verification_sha256,'native verification pin differs')
        require(verification['status']=='complete','native verification incomplete')
        verified={str(Path(v['path']).resolve()):v for v in verification['verified_files'].values()}
        for filename in ('summary.json','status.json','evidence.json','regional-native-context.json'):
            path=native/filename;require(str(path) in verified,'native artifact absent from independent verification');load(path,verified[str(path)])
        summary=json.loads((native/'summary.json').read_text());status=json.loads((native/'status.json').read_text());evidence=json.loads((native/'evidence.json').read_text())
        require(summary['status']==status['status']=='complete','native case incomplete')
        require(status['summary_sha256']==pins[str(native/'summary.json')]['sha256'],'native summary/status binding differs')
        for name in ('evidence.json','regional-native-context.json'):
            same(summary['outputs'][name]['sha256'],pins[str(native/name)]['sha256'],'native output binding differs')
        same(summary['inputs']['config']['sha256'],selection['native_config_pin']['sha256'],'native configuration pin differs')
        same(evidence['config'],config,'native configuration payload differs')
        saved=json.loads((native/'regional-native-context.json').read_text())['targets'][spec['target_id']]
        native_contexts={row['context']['source_record_id']:row['context'] for row in saved}
        require(len(native_contexts)==len(saved),'duplicate native context');same(native_contexts,contexts,'native context changed from preparation')
        preflight={w['flank_bp']:w for w in selection['windows']};require(len(evidence['windows'])==len(preflight),'window membership differs')
        windows=[]
        for window in evidence['windows']:
            pre=preflight[window['flank_bp']]
            for key in ('chrom','start0','end0','reference_sequence'): same(window[key],pre[key],'window '+key+' differs')
            old=[evidence['hypotheses'][i] for i in window['hypothesis_ids']]
            same({h['sequence']:h['source_combinations'] for h in old},{h['sequence']:h['source_combinations'] for h in pre['native_family']['haplotypes']},'native sequence lineages differ')
            same(window['enumeration_metadata']['rejected_combinations'],pre['native_family']['rejected_combinations'],'native overlap refusals differ')
            same(window['omitted_overlapping_native_alleles'],pre['omitted_native_alleles'],'nonconcrete native omissions differ')
            overlay=overlay_lineages(window['reference_sequence'],window['start0'],spec,contexts)
            classes=union_classes(old,overlay['lineages'],'ADAMTS1_F'+str(window['flank_bp']))
            expected_arms=[(mq,bq) for mq in config['mapping_quality_thresholds'] for bq in config['base_quality_thresholds']]
            same([(a['minimum_mapping_quality'],a['minimum_base_quality']) for a in window['arms']],expected_arms,'quality arm membership/order differs')
            arms=[]
            for arm in window['arms']:
                observed=arm['evidence'];mq,bq=arm['minimum_mapping_quality'],arm['minimum_base_quality']
                same(observed['window'],{'start0':window['start0'],'end0':window['end0']},'observation window differs')
                same(observed['thresholds'],{'min_MQ':mq,'min_raw_BQ':bq,'excluded_flags':0xF04,'reject_MQ_255':True,'BAQ_applied':False},'observer threshold policy differs')
                arms.append({'minimum_mapping_quality':mq,'minimum_base_quality':bq,**recount(observed,classes)})
            windows.append({'chrom':window['chrom'],'start0':window['start0'],'end0':window['end0'],'flank_bp':window['flank_bp'],
                'original_class_count':len(old),'union_class_count':len(classes),'classes':list(classes.values()),'overlay':overlay,'arms':arms})
        for pin in list(pins.values()): same(fingerprint(pin['path']),pin,'pinned input changed')
        result={'status':'complete','completed_at':stamp(),'verified_files':pins,'native_contexts':copy.deepcopy(contexts),'windows':windows,
            'limitations':['Predeclared normalized lineages reuse exactly completed original observations and quality gates; no new extraction, alignment, genotyper or annotation.',
                'Both derived sequences already exist in the original panel. Lineages are added without new sequence classes or doubled observation/fragment support.',
                'Native GT/ALT/FORMAT including spanning star are preserved. Sequence identity is not biological phase, genotype, locus origin, function or a clinical finding.']}
        save(output/'evidence.json',result)
        bounded={'status':'complete','completed_at':result['completed_at'],'evidence':fingerprint(output/'evidence.json'),'verified_files':pins,
            'windows':[{'flank_bp':w['flank_bp'],'original_class_count':w['original_class_count'],'union_class_count':w['union_class_count'],
                'arms':[{'MQ':a['minimum_mapping_quality'],'raw_BQ':a['minimum_base_quality'],'counts':a['counts'],'class_support':a['class_support']} for a in w['arms']]} for w in windows],
            'limitations':result['limitations']}
        save(output/'summary.json',bounded);state.update(status='complete',completed_at=stamp(),summary_sha256=fingerprint(output/'summary.json')['sha256']);save(output/'status.json',state)
    except BaseException as error:
        state.update(status='failed',failed_at=stamp(),error_type=type(error).__name__,error=str(error));save(output/'status.json',state);raise

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,required=True);parser.add_argument('--native-run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--selection',type=Path,required=True);parser.add_argument('--expected-selection-sha256',required=True)
    parser.add_argument('--controls',type=Path,required=True);parser.add_argument('--expected-controls-sha256',required=True)
    parser.add_argument('--native-verification',type=Path,required=True);parser.add_argument('--expected-native-verification-sha256',required=True);run(parser.parse_args())
