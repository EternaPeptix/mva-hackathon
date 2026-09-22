#!/usr/bin/env python3
"""Public clinical-drug census. Python 3.9+, standard library only; see README.md."""
import argparse, csv, gzip, hashlib, io, json, re, time, urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CHEMBL = 'https://www.ebi.ac.uk/chembl/api/data/'
LICENSES = {
    'ChEMBL': {'license': 'CC-BY-SA-3.0', 'url': 'https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/general-questions'},
    'DrugCentral': {'license': 'CC-BY-SA-4.0', 'url': 'https://drugcentral.org/privacy'},
    'UniProt': {'license': 'CC-BY-4.0', 'url': 'https://www.uniprot.org/help/license'},
}
MOLECULE_FIELDS = 'molecule_chembl_id,pref_name,max_phase,withdrawn_flag,first_approval,molecule_hierarchy,molecule_synonyms,molecule_type,molecule_structures,cross_references,veterinary,therapeutic_flag,atc_classifications'

def stamp(): return datetime.now(timezone.utc).isoformat()
def js(x): return json.dumps(x, ensure_ascii=False, sort_keys=True, separators=(',', ':'))
def write_json(path, data): path.write_text(json.dumps(data, indent=2, ensure_ascii=False, sort_keys=True)+'\n')
def norm(x): return ' '.join((x or '').casefold().split())
def sha(b): return hashlib.sha256(b).hexdigest()

class Fetch:
    def __init__(self, cache, offline=False):
        self.raw = cache; self.raw.mkdir(parents=True, exist_ok=True)
        self.offline = offline; self.used = {}
    def get(self, name, url):
        p=self.raw/(name+'.gz'); m=self.raw/(name+'.meta.json')
        if p.exists() and m.exists():
            meta=json.loads(m.read_text())
            if meta['url'] != url: raise ValueError(f'Cache URL mismatch: {name}')
            b=gzip.decompress(p.read_bytes())
            if sha(b) != meta['sha256']: raise ValueError(f'Cache checksum mismatch: {name}')
        else:
            if self.offline: raise RuntimeError(f'Missing cached source: {name}')
            for attempt in range(5):
                try:
                    req=urllib.request.Request(url,headers={'User-Agent':'MVA-public-drug-census/1.0','Accept-Encoding':'identity'})
                    with urllib.request.urlopen(req,timeout=150) as response:
                        b=response.read()
                        meta={'url':url,'resolved_url':response.url,'retrieved_utc':stamp(),
                              'last_modified':response.headers.get('Last-Modified'),
                              'etag':response.headers.get('ETag'),'sha256':sha(b),'bytes':len(b)}
                    break
                except Exception:
                    if attempt==4: raise
                    time.sleep(2**attempt)
            p.write_bytes(gzip.compress(b,mtime=0)); write_json(m,meta)
        self.used[name]=meta
        return b
    def json(self,name,url): return json.loads(self.get(name,url))
    def pages(self,name,endpoint,key,params=None):
        from urllib.parse import urlencode
        params=dict(params or {}); params['limit']=1000
        def page(offset):
            u=CHEMBL+endpoint+'.json?'+urlencode(dict(params,offset=offset))
            return self.json(f'{name}-{offset:06d}.json',u)
        first=page(0); total=first['page_meta']['total_count']
        offsets=range(1000,total,1000)
        with ThreadPoolExecutor(max_workers=4) as pool: rest=list(pool.map(page,offsets))
        pages=[first]+rest
        assert all(p['page_meta']['total_count']==total for p in pages), 'Source changed during pagination'
        rows=[r for p in pages for r in p[key]]
        assert len(rows)==total,(name,len(rows),total)
        print(f'{name}: {len(rows)} / {total}',flush=True)
        return rows,total

def write_tsv(path,rows,fields):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader()
        for row in rows:
            w.writerow({k:js(v) if isinstance(v,(dict,list)) else '' if v is None else v for k,v in row.items()})

def dc_data(fetch):
    sdf=gzip.decompress(fetch.get('drugcentral-structures-2023.sdf.gz','https://unmtid-dbs.net/download/DrugCentral/2023/structures.molV3.sdf.gz')).decode()
    records={}
    for block in sdf.split('$$$$'):
        if not block.strip(): continue
        props={m.group(1):m.group(2).strip() for m in re.finditer(r'>\s*<([^>]+)>[^\n]*\n(.*?)(?=\n\n|\Z)',block,re.S)}
        assert props.get('ID'), 'SDF record lacks ID'
        sid=props['ID']; assert sid not in records
        records[sid]={'drugcentral_id':sid,'name':props.get('PREFERRED_NAME',''), 'cas_rn':props.get('CAS_RN',''),
                      'aliases':props.get('SYNONYMS','').splitlines(),'inchikey':'','approval_agencies':[],'sources':['2023_SDF']}
    smiles=fetch.get('drugcentral-structures-2021.tsv','https://unmtid-dbs.net/download/DrugCentral/2021_09_01/structures.smiles.tsv').decode()
    smiles_rows=list(csv.DictReader(io.StringIO(smiles),delimiter='\t'))
    for row in smiles_rows:
        sid=row['ID']
        if sid not in records:
            records[sid]={'drugcentral_id':sid,'name':row['INN'],'cas_rn':row['CAS_RN'],'aliases':[row['INN']], 'approval_agencies':[], 'sources':[]}
        records[sid]['inchikey']=row['InChIKey'];records[sid]['sources'].append('2021_SMILES')
    approval_counts={}
    for agency in ['FDA','EMA','PMDA']:
        text=fetch.get('drugcentral-'+agency+'.csv',f'https://drugcentral.org/static/{agency}_Approved.csv').decode()
        rows=list(csv.reader(io.StringIO(text)));approval_counts[agency]=len(rows)
        for sid,name in rows:
            if sid not in records:
                records[sid]={'drugcentral_id':sid,'name':name,'cas_rn':'','aliases':[name],'inchikey':'','approval_agencies':[],'sources':[]}
            records[sid]['approval_agencies'].append(agency);records[sid]['sources'].append(agency+'_CSV')
    interactions=list(csv.DictReader(io.StringIO(gzip.decompress(fetch.get('drugcentral-interactions-2021.tsv.gz','https://unmtid-dbs.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz')).decode()),delimiter='\t'))
    # Target export has a small number of IDs absent from the newer structure export.
    # Preserve these as status-unverified records rather than silently dropping them.
    for row in interactions:
        sid=row['STRUCT_ID']
        if sid not in records:
            records[sid]={'drugcentral_id':sid,'name':row['DRUG_NAME'],'cas_rn':'','aliases':[row['DRUG_NAME']],'inchikey':'','approval_agencies':[],'sources':['2021_target_interaction_only']}
    return records,interactions,approval_counts,len(smiles_rows)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--offline',action='store_true');p.add_argument('--out',type=Path,default=ROOT);p.add_argument('--cache',type=Path,default=Path.home()/'.cache'/'mva-track2-drug-census');args=p.parse_args()
    out=args.out;out.mkdir(parents=True,exist_ok=True);fetch=Fetch(args.cache,args.offline)
    before=fetch.json('chembl-status-before.json',CHEMBL+'status.json')
    clinical,nclinical=fetch.pages('chembl-clinical','molecule','molecules',{'max_phase__gte':1,'only':MOLECULE_FIELDS,'order_by':'molecule_chembl_id'})
    withdrawn,nwithdrawn=fetch.pages('chembl-withdrawn','molecule','molecules',{'withdrawn_flag':'true','only':MOLECULE_FIELDS,'order_by':'molecule_chembl_id'})
    molecules={r['molecule_chembl_id']:r for r in clinical}
    assert len(molecules)==nclinical
    for r in withdrawn: molecules.setdefault(r['molecule_chembl_id'],r)
    mechanisms,nmechanisms=fetch.pages('chembl-mechanisms','mechanism','mechanisms',{'order_by':'mec_id'})
    targets,ntargets=fetch.pages('chembl-targets','target','targets',{'only':'target_chembl_id,pref_name,organism,target_type,target_components','order_by':'target_chembl_id'})
    assert len({r['target_chembl_id'] for r in targets})==ntargets
    dc,dc_interactions,approval_counts,nsmiles=dc_data(fetch)
    uni=fetch.json('uniprot-O60566.json','https://rest.uniprot.org/uniprotkb/O60566.json')
    assert uni['primaryAccession']=='O60566' and uni['genes'][0]['geneName']['value']=='BUB1B' and uni['organism']['taxonId']==9606
    ensembl_refs=[x for x in uni.get('uniProtKBCrossReferences',[]) if x['database']=='Ensembl']
    assert 'ENSG00000156970' in js(ensembl_refs)
    identity={'gene':'BUB1B','uniprot':'O60566','ensembl':'ENSG00000156970','taxon':9606,'audit':uni['entryAudit'],'ensembl_cross_references':ensembl_refs}
    write_json(out/'gene-identity.json',identity)
    # Exact human interacting-protein neighborhood is annotation, never rescue evidence.
    neighbors={}
    for comment in uni['comments']:
        if comment['commentType']=='INTERACTION':
            for row in comment['interactions']:
                if not row.get('organismDiffer'):
                    it=row['interactantTwo'];neighbors[it['uniProtKBAccession'].split('-')[0]]=it.get('geneName','')
    write_json(out/'bub1b-uniprot-interactors.json',neighbors)
    tindex={t['target_chembl_id']:t for t in targets}
    accessions=lambda t: sorted({c['accession'] for c in t.get('target_components',[]) if c.get('accession')})
    btids=[t['target_chembl_id'] for t in targets if 'O60566' in accessions(t) and t['organism']=='Homo sapiens']
    activities=[]
    if btids:
        activities,_=fetch.pages('chembl-bub1b-activities','activity','activities',{'target_chembl_id__in':','.join(sorted(btids)),'order_by':'activity_id'})
    after=fetch.json('chembl-status-after.json',CHEMBL+'status.json?after=census')
    assert before['chembl_db_version']==after['chembl_db_version'],'ChEMBL release changed during retrieval'
    # Parent hierarchy joins salts/forms but retains each source molecule explicitly.
    groups={}; molecule_to_group={}; aliases=set()
    def blank(cid,name):
        return {'catalog_id':cid,'preferred_name':name,'chembl_ids':[],'drugcentral_ids':[],
                'max_phase':None,'withdrawn_flag':False,'first_approval_years':[], 'molecule_types':[],
                'veterinary_flags':[],'drugcentral_approval_agencies':[],'inclusion_basis':[],
                'inchikeys':[], 'atc_codes':[], 'chembl_mechanism_count':0,'drugcentral_interaction_count':0,
                'drugcentral_moa_count':0,'bub1b_activity_count':0}
    for mid,m in sorted(molecules.items()):
        gid=(m.get('molecule_hierarchy') or {}).get('parent_chembl_id') or mid
        molecule_to_group[mid]=gid
        g=groups.setdefault(gid,blank(gid,m.get('pref_name') or mid))
        if gid==mid and m.get('pref_name'):g['preferred_name']=m['pref_name']
        g['chembl_ids'].append(mid)
        phase=float(m['max_phase']) if m.get('max_phase') is not None else None
        if phase is not None:g['max_phase']=max(g['max_phase'] or 0,phase)
        if phase is not None and phase>=1:g['inclusion_basis'].append('ChEMBL_max_phase_ge_1')
        if m.get('withdrawn_flag'):g['withdrawn_flag']=True;g['inclusion_basis'].append('ChEMBL_withdrawn_flag')
        if m.get('first_approval'):g['first_approval_years'].append(m['first_approval'])
        if m.get('molecule_type'):g['molecule_types'].append(m['molecule_type'])
        if m.get('veterinary') is not None:g['veterinary_flags'].append(m['veterinary'])
        ik=(m.get('molecule_structures') or {}).get('standard_inchi_key')
        if ik:g['inchikeys'].append(ik)
        g['atc_codes']+=m.get('atc_classifications') or []
        if m.get('pref_name'):aliases.add((gid,'ChEMBL',mid,'preferred',m['pref_name']))
        for a in m.get('molecule_synonyms') or []:
            name=a.get('molecule_synonym') or a.get('synonyms')
            if name:aliases.add((gid,'ChEMBL',mid,a.get('syn_type',''),name))
    # Conservative DC crosswalk: unique exact full InChIKey, otherwise unique exact preferred name;
    # reject preferred-name joins with conflicting known structure keys. No fuzzy brand/salt merging.
    ik_index=defaultdict(set);name_index=defaultdict(set)
    for gid,g in groups.items():
        for ik in g['inchikeys']:ik_index[ik].add(gid)
        name_index[norm(g['preferred_name'])].add(gid)
    dc_to_group={};crosswalk=[];supplement=[]
    for sid,d in sorted(dc.items(),key=lambda x:int(x[0])):
        match='unmatched';candidates=ik_index.get(d['inchikey'],set()) if d.get('inchikey') else set()
        if len(candidates)==1:match='unique_exact_full_inchikey';gid=next(iter(candidates))
        else:
            nc=name_index.get(norm(d['name']),set())
            if len(nc)==1:
                ng=next(iter(nc));known=groups[ng]['inchikeys']
                if not d.get('inchikey') or not known or d['inchikey'] in known:
                    gid=ng;match='unique_exact_preferred_name'
                else:match='name_structure_conflict'
            elif len(candidates)>1 or len(nc)>1:match='ambiguous'
        if match in ('unmatched','ambiguous','name_structure_conflict'):
            gid='DC'+sid
            if d['approval_agencies']:groups[gid]=blank(gid,d['name'])
        included=gid in groups
        dc_to_group[sid]=gid
        crosswalk.append({'drugcentral_id':sid,'catalog_id':gid if included else '', 'name':d['name'],'match_basis':match,'included':included,'inchikey':d.get('inchikey','')})
        if not included:
            supplement.append(dict(d,exclusion_reason='No ChEMBL clinical/withdrawal match and no FDA/EMA/PMDA approval-list membership; status unverified'))
            continue
        g=groups[gid];g['drugcentral_ids'].append(sid);g['drugcentral_approval_agencies']+=d['approval_agencies']
        if d['approval_agencies']:g['inclusion_basis'].append('DrugCentral_regulator_approval_list')
        if d.get('inchikey'):g['inchikeys'].append(d['inchikey'])
        for name in set(d['aliases']+[d['name']]):
            if name:aliases.add((gid,'DrugCentral',sid,'name_or_synonym',name))
    mrows=[];screen=[]
    for m in mechanisms:
        mid=m['molecule_chembl_id'];gid=molecule_to_group.get(mid) or molecule_to_group.get(m.get('parent_molecule_chembl_id'))
        if not gid:continue
        t=tindex.get(m.get('target_chembl_id'),{});acc=accessions(t)
        row=dict(m,catalog_id=gid,drug_name=groups[gid]['preferred_name'],target_name=t.get('pref_name'),target_organism=t.get('organism'),target_type=t.get('target_type'),target_accessions=acc)
        mrows.append(row);groups[gid]['chembl_mechanism_count']+=1
        if t.get('organism')=='Homo sapiens':
            overlap=set(acc)&set(neighbors)
            if 'O60566' in acc or overlap:
                screen.append({'catalog_id':gid,'drug_name':groups[gid]['preferred_name'],'source':'ChEMBL','record_id':m['mec_id'],'evidence_kind':'curated_mechanism','target_name':t.get('pref_name'),'target_accessions':acc,'action_type':m.get('action_type'),'bub1b_direct_target':'O60566' in acc,'bub1b_neighbor_targets':sorted(overlap),'references':m.get('mechanism_refs'),'interpretation':'Target annotation only; direction, exposure and functional benefit unestablished'})
    dcrows=[]
    for idx,r in enumerate(dc_interactions,1):
        gid=dc_to_group.get(r['STRUCT_ID'])
        if not gid:continue
        moa=r.get('MOA')=='1';acc=r['ACCESSION'].split('|')
        row=dict(r,catalog_id=gid,source_record_number=idx,catalog_scope='clinical_or_approved' if gid in groups else 'status_unverified_supplement');dcrows.append(row)
        if gid in groups:
            groups[gid]['drugcentral_interaction_count']+=1
            if moa:groups[gid]['drugcentral_moa_count']+=1
        overlap=set(acc)&set(neighbors)
        if r['ORGANISM']=='Homo sapiens' and ('O60566' in acc or overlap):
            screen.append({'catalog_id':gid,'drug_name':groups[gid]['preferred_name'] if gid in groups else dc[r['STRUCT_ID']]['name'],'source':'DrugCentral','record_id':idx,'evidence_kind':'curated_mechanism' if moa else 'interaction_or_assay','target_name':r['TARGET_NAME'],'target_accessions':acc,'action_type':r['ACTION_TYPE'],'bub1b_direct_target':'O60566' in acc,'bub1b_neighbor_targets':sorted(overlap),'references':[r['ACT_SOURCE_URL'],r['MOA_SOURCE_URL']],'interpretation':'Target annotation only; direction, exposure and functional benefit unestablished'})
    arows=[]
    for a in activities:
        gid=molecule_to_group.get(a['molecule_chembl_id'])
        if gid:
            groups[gid]['bub1b_activity_count']+=1;arows.append(dict(a,catalog_id=gid))
    for g in groups.values():
        for key,value in list(g.items()):
            if isinstance(value,list):g[key]=sorted(set(value))
        g['status_summary']='withdrawal_flagged' if g['withdrawn_flag'] else 'approval_record' if g['max_phase']==4 or g['drugcentral_approval_agencies'] else 'clinical_development_record'
        g['annotation_status']='has_target_or_mechanism_annotation' if g['chembl_mechanism_count'] or g['drugcentral_interaction_count'] or g['bub1b_activity_count'] else 'unannotated_in_retrieved_mechanism_sources'
    grows=sorted(groups.values(),key=lambda r:r['catalog_id'])
    write_tsv(out/'catalog.tsv',grows,list(grows[0]))
    write_tsv(out/'aliases.tsv',[dict(zip(['catalog_id','source','source_id','alias_type','alias'],r)) for r in sorted(aliases)],['catalog_id','source','source_id','alias_type','alias'])
    write_tsv(out/'drugcentral-crosswalk.tsv',crosswalk,list(crosswalk[0]))
    write_tsv(out/'drugcentral-status-unverified.tsv',supplement,list(supplement[0]) if supplement else ['drugcentral_id'])
    write_tsv(out/'chembl-mechanisms.tsv',mrows,list(mrows[0]) if mrows else ['catalog_id'])
    write_tsv(out/'drugcentral-interactions.tsv',dcrows,list(dcrows[0]) if dcrows else ['catalog_id'])
    write_tsv(out/'bub1b-target-neighborhood.tsv',screen,list(screen[0]) if screen else ['catalog_id'])
    write_tsv(out/'bub1b-all-target-activities.tsv',activities,list(activities[0]) if activities else ['activity_id'])
    write_tsv(out/'bub1b-clinical-drug-activities.tsv',arows,list(arows[0]) if arows else ['catalog_id'])
    referenced_targets={m['target_chembl_id'] for m in mrows}|set(btids)
    target_rows=[{k:t.get(k) for k in ['target_chembl_id','pref_name','organism','target_type']} | {'target_accessions':accessions(t),'gene_symbols':sorted({x['component_synonym'] for c in t.get('target_components',[]) for x in c.get('target_component_synonyms',[]) if x.get('syn_type')=='GENE_SYMBOL'})} for t in targets if t['target_chembl_id'] in referenced_targets]
    write_tsv(out/'chembl-targets.tsv',target_rows,list(target_rows[0]))
    # Keep all source-level identifiers/statuses, including clinical salts collapsed in catalog.tsv.
    source_rows=[{k:v for k,v in m.items() if k not in ['molecule_structures','molecule_synonyms']} for m in molecules.values()]
    write_tsv(out/'chembl-clinical-records.tsv',sorted(source_rows,key=lambda r:r['molecule_chembl_id']),list(source_rows[0]))
    coverage={'built_utc':stamp(),'chembl_release':before,'source_scope':'ChEMBL max_phase>=1 OR withdrawn_flag; union DrugCentral FDA/EMA/PMDA approval lists',
      'source_counts':{'chembl_phase_ge_1':nclinical,'chembl_withdrawn_flag':nwithdrawn,'chembl_unique_clinical_or_withdrawn_records':len(molecules),'chembl_mechanisms_all':nmechanisms,'chembl_targets_all':ntargets,'drugcentral_record_union':len(dc),'drugcentral_interaction_only_ids':sum(d['sources']==['2021_target_interaction_only'] for d in dc.values()),'drugcentral_smiles_2021':nsmiles,'drugcentral_approval_list_rows':approval_counts,'drugcentral_unique_approved_list_ids':sum(bool(d['approval_agencies']) for d in dc.values()),'drugcentral_interactions_all':len(dc_interactions),'drugcentral_interactions_screened':len(dcrows),'drugcentral_interactions_status_unverified':sum(r['catalog_scope']=='status_unverified_supplement' for r in dcrows),'bub1b_chembl_activity_records_all':len(activities)},
      'catalog_counts':{'parent_or_unmatched_dc_entities':len(grows),'chembl_parent_groups':sum(bool(g['chembl_ids']) for g in grows),'drugcentral_only_approved_entities':sum(not g['chembl_ids'] for g in grows),'drugcentral_status_unverified_excluded':len(supplement),'aliases':len(aliases),'status_summary':dict(Counter(g['status_summary'] for g in grows)),'maximum_phase':dict(Counter(str(g['max_phase']) for g in grows)),'withdrawal_flagged_entities':sum(g['withdrawn_flag'] for g in grows),'veterinary_flagged_entities':sum(1 in g['veterinary_flags'] for g in grows)},
      'coverage':{'chembl_mechanism_entities':sum(g['chembl_mechanism_count']>0 for g in grows),'drugcentral_interaction_entities':sum(g['drugcentral_interaction_count']>0 for g in grows),'drugcentral_moa_entities':sum(g['drugcentral_moa_count']>0 for g in grows),'any_target_or_mechanism_annotated_entities':sum(g['annotation_status']=='has_target_or_mechanism_annotation' for g in grows),'unannotated_entities':sum(g['annotation_status']=='unannotated_in_retrieved_mechanism_sources' for g in grows),'bub1b_annotated_targets':btids,'direct_bub1b_mechanism_or_interaction_rows':sum(r['bub1b_direct_target'] for r in screen),'bub1b_neighborhood_rows':len(screen),'bub1b_neighborhood_entities':len({r['catalog_id'] for r in screen}),'bub1b_clinical_drug_activity_rows':len(arows),'bub1b_clinical_drug_activity_entities':len({r['catalog_id'] for r in arows})},
      'crosswalk':dict(Counter(r['match_basis'] for r in crosswalk)),
      'limitations':['Source-bounded census, not every drug ever developed or every medicine worldwide.','Maximum development phase is historical, not proof of current authorization, development activity, jurisdiction, formulation, pediatric safety, or efficacy.','Withdrawn flag does not mean withdrawal everywhere; no false-to-currently-marketed inference.','DrugCentral 2023 SDF / 2021 SMILES and interactions have different dates; approval CSV content version is unreported, HTTP metadata retained.','Status-unverified DrugCentral entries outside clinical/withdrawal/approval filters are retained separately, not promoted to medicines by chemical identity alone.','Parent-group counts retain molecular/formulation heterogeneity and unresolved cross-source duplicates; not a definitive count of distinct active ingredients.','Target annotations, assay records and interactions do not establish pharmacological benefit or mutant rescue.','Missing annotations are UNKNOWN, not evidence of inefficacy.','No whole-transcriptome perturbation, all-drug docking, clinical intervention, or experiment performed.']}
    write_json(out/'coverage.json',coverage)
    write_json(out/'provenance.json',{'licenses':LICENSES,'sources':dict(sorted(fetch.used.items())),'gene_identity':identity,'source_version_check':{'before':before,'after':after}})
    print(json.dumps(coverage,indent=2),flush=True)

if __name__=='__main__':main()
