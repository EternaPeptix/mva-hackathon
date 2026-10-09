#!/usr/bin/env python3
"""Replay NIH DSLD ingredient groups and all-alias discovery; no efficacy or reading credit.

Inputs are immutable. Outputs are staged, inputs rehashed, then an immutable run
is published through derived/current (an atomically replaced relative symlink).
Consumers must resolve current once and validate summary.json output hashes.
Only --out is written. No network, fuzzy names, chemical equivalence or raw
abstract output. A source synonym relation is preserved, never authenticated.
"""
import argparse
import csv
import hashlib
import importlib.util
import io
import json
import os
import re
import shutil
import string
import sys
import tempfile
import unicodedata
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import parse_qs, urlparse

sys.dont_write_bytecode = True
WORK = Path('/Users/jeweled/mva-hackathon/track2-drug-census-work')
DEFAULT_SOURCE = WORK / 'supplement-ingredients-20260923'
DEFAULT_CATALOG = Path('/Users/jeweled/mva-hackathon/artifacts/source-revalidation-20260922/repo/workflow/track2-drug-census')
MANIFEST_SHA = '28409b93b597b4fc1633a47646c188fb530b1250b4a479f0324952cdaf911ad3'
CORPUS_SHA = 'e5b8b017bb68af38818f1cc8c1bee47adb16afe0367981104a6b17886d2954d9'
PARTITION_TOTALS = dict(zip(list(string.ascii_uppercase)+['Other'], [502,484,838,222,233,198,358,217,149,56,71,295,345,119,123,606,17,174,581,284,26,254,105,18,36,28,128]))
HIT_KEYS = {'_id','_index','_type','_score','_source'}
GROUP_KEYS = {'groupId','groupName','category','synonyms','factsheets','nutrientInfo'}
WORD = re.compile(r'[^\W_]+', re.UNICODE)
HYPHENS = set('-‐‑‒–—−')
GENE = re.compile(r'\b(?:bub1b|bubr1|bub r1|mosaic variegated aneuploidy)\b')
# Curated ambiguity flags, not exclusions or a claim that the list is exhaustive.
COMMON = set('''acid alcohol antioxidant antioxidants amino amino acids ash base bile bio black blood blue bone brain calcium carbon cell cellulose chloride cholesterol chromium collagen copper core culture enzyme enzymes essential extract extracts fat fatty fiber fibre folate glucose gold green growth gum hormone hormones human immune iron lactate lead lipid lipids magnesium manganese matrix milk mineral minerals natural nitrogen oil oils oxygen peptide peptides phosphorus plant plants potassium powder protein proteins red resin root salt selenium serum silicon silver sodium starch sulfur support tissue total vitamin vitamins water white yeast zinc'''.split()) | {'amino acid','amino acids','fatty acid','fatty acids','green tea','essential oil','trace minerals','dietary fiber','dietary fibre','casein kinase'}
# Added after bounded inspection of the highest-frequency lexical hits; never exclude them.
COMMON |= set('complex same threonine phosphatase muscle skin lysine tyrosine insulin estrogen bladder cultures gamma fibers none paraffin cytochrome galactosidase glutathione mono adenine alanine cysteine camp leucine bacteria fruit adenosine hydrogen casein biotin glycine arginine aspartate glutamate serine histidine isoleucine methionine phenylalanine proline tryptophan valine lung liver prostate date'.split()) | {'stem cells','growth factors','transforming growth factor','other proteins','nicotinamide adenine dinucleotide'}
METADATA_FIELDS = ['id','source','pmid','pmcid','doi','title','pubTypeList','publicationStatus','commentCorrectionList','versionList','versionNumber','firstPublicationDate']


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def sha(path):
    return digest(path.read_bytes())


def jdump(value):
    return json.dumps(value, ensure_ascii=False, separators=(',', ':'))


def checked(path, expected, snapshots):
    path = path.resolve()
    data = path.read_bytes()
    actual = digest(data)
    require(actual == expected, 'Input hash mismatch: '+str(path))
    snapshots[str(path)] = actual
    return data


def capture(path, snapshots):
    path = path.resolve();data = path.read_bytes()
    snapshots[str(path)] = digest(data)
    return data


def recheck(snapshots):
    for path, expected in snapshots.items():
        require(sha(Path(path)) == expected, 'Input changed before publish: '+path)


def parse_page(payload, expected_total, partition):
    require(isinstance(payload, dict) and set(payload)=={'total','max_score','hits'}, 'Response schema differs: '+partition)
    require(payload['total']=={'value':expected_total,'relation':'eq'}, 'Total not exact/expected: '+partition)
    require(isinstance(payload['hits'], list), 'hits is not a list')
    rows=[]
    for index, hit in enumerate(payload['hits']):
        require(isinstance(hit, dict) and set(hit)==HIT_KEYS, 'Hit schema differs')
        source=hit['_source']
        require(isinstance(source,dict) and set(source)==GROUP_KEYS, 'Group schema differs')
        require(isinstance(hit['_id'],str) and bool(hit['_id']) and hit['_id']==source['groupId'], 'ID/groupId mismatch')
        require(isinstance(source['groupName'],str), 'Nonstring group name')
        require(isinstance(source['synonyms'],list) and all(isinstance(x,str) for x in source['synonyms']), 'Nonstring synonym cell')
        require(isinstance(source['category'],list) and all(x is None or isinstance(x,str) for x in source['category']), 'Unexpected category type')
        require(isinstance(source['factsheets'],list) and all(isinstance(x,dict) and set(x)=={'link','name'} and all(isinstance(v,str) for v in x.values()) for x in source['factsheets']), 'Factsheet schema differs')
        require(isinstance(source['nutrientInfo'],list), 'NutrientInfo is not a list')
        require(isinstance(hit['_index'],str) and isinstance(hit['_type'],str) and isinstance(hit['_score'],(float,int)), 'Hit metadata type differs')
        rows.append((index,hit))
    return rows


def require_full_page(page, expected_total, offset, page_size, partition):
    require(len(page)==min(page_size,expected_total-offset), 'Partial page: '+partition)


def register_group(ids, hit):
    gid=hit['_id']
    require(gid not in ids, 'Duplicate/conflicting group ID: '+gid)
    ids[gid]=hit


def load_source(source_dir, snapshots):
    source_dir=source_dir.resolve()
    manifest=json.loads(checked(source_dir/'retrieval-manifest.json',MANIFEST_SHA,snapshots))
    require(manifest['schema_version']==1, 'Manifest schema version differs')
    parts=manifest['partitions']
    require([p['partition'] for p in parts]==list(PARTITION_TOTALS), 'Expected all 27 ordered partitions')
    require(len(manifest['version_checks'])==2, 'Need before and after version checks')
    versions=[]
    for v in manifest['version_checks']:
        data=checked(source_dir/v['file'],v['sha256'],snapshots)
        require(v['status']==200 and len(data)==v['bytes'], 'Version response not successful/exact')
        versions.append(json.loads(data))
    require(versions[0]==versions[1]==manifest['api_version'], 'API versions changed')
    require(versions[0]['version']=='9.5.0' and versions[0]['versionTimeStamp']=='May 2026', 'Unexpected API snapshot')
    require(versions[0]['esIndexProcessed']=='dsldnxt_labels-2025-09-25n and dsldnxt_ing_groups-2025-09-25n', 'Processed index differs')
    rows=[];ids={};rawfiles=[];part_summary=[]
    for part in parts:
        letter=part['partition'];expected=PARTITION_TOTALS[letter];collected=[];offset=0
        require(part['declared_total']==part['retrieved_rows']==expected, 'Manifest total differs: '+letter)
        for response in part['responses']:
            file=(source_dir/response['file']).resolve()
            require(file.is_relative_to(source_dir/'raw'), 'Response path escapes raw source')
            require(response['status']==200, 'Failed group response in manifest')
            url=urlparse(response['url']);q=parse_qs(url.query)
            require(url.scheme=='https' and url.netloc=='api.ods.od.nih.gov' and url.path.rstrip('/')=='/dsld/v9/ingredient-groups', 'Unexpected group endpoint')
            require(q.get('method')==['by_letter'] and q.get('term')==[letter] and q.get('from')==[str(offset)], 'Pagination/partition URL differs')
            page_size=int(q['size'][0]);require(page_size>0, 'Invalid page size')
            data=checked(file,response['sha256'],snapshots);require(len(data)==response['bytes'], 'Response byte count differs')
            payload=json.loads(data);page=parse_page(payload,expected,letter)
            require_full_page(page,expected,offset,page_size,letter)
            for local,hit in page:
                gid=hit['_id'];register_group(ids,hit)
                rows.append({'partition':letter,'source_file':response['file'],'source_hit_index':local,'partition_row_index':offset+local,'source_url':response['url'],'hit':hit})
                collected.append(gid)
            offset+=len(page);rawfiles.append(file)
        require(offset==expected and collected==part['source_ids'], 'Source IDs/order/count differs: '+letter)
        part_summary.append({'partition':letter,'declared_exact_total':expected,'retrieved_rows':offset,'responses':len(part['responses'])})
    require(set((source_dir/'raw').glob('groups-*.json'))==set(rawfiles), 'Missing/unmanifested group file')
    require(len(rows)==len(ids)==manifest['sum_partition_rows']==manifest['distinct_group_ids']==6467, 'Group denominator differs')
    require(manifest['cross_partition_duplicate_ids']=={}, 'Manifest reports duplicate groups')
    require(sum(len(r['hit']['_source']['synonyms']) for r in rows)==78788, 'Synonym-cell denominator differs')
    require(all(r['hit']['_source']['nutrientInfo']==[] for r in rows), 'NutrientInfo changed from this snapshot')
    return rows,manifest,part_summary


def load_clean(catalog_dir,snapshots):
    path=catalog_dir/'broad_hub_screen.py';capture(path,snapshots)
    spec=importlib.util.spec_from_file_location('dsld_corrected_broad_screen',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    clean=module.clean
    require('GSK923295' in clean('A <0.1 GSK923295 > 2'), 'Unsafe legacy text stripping detected')
    require(clean('&lt;i&gt;test&lt;/i&gt;').strip()=='test','HTML decode/tag policy differs')
    return clean


def normal_text(text,clean):
    return unicodedata.normalize('NFKC',clean(text)).casefold()


def norm(text,clean):
    return ' '.join(WORD.findall(normal_text(text,clean)))


def ambiguity(key,groups):
    flags=[]
    if not key:flags.append('empty_after_normalization_not_searchable')
    if key and len(key.replace(' ',''))<=3:flags.append('short_alias_3_or_fewer_characters')
    if key and all(t.isdigit() for t in key.split()):flags.append('numeric_only_alias')
    if key in COMMON:flags.append('manual_common_or_class_term')
    if len(groups)>1:flags.append('alias_shared_across_DSLD_groups')
    return flags


def make_trie(keys):
    trie={}
    for key in keys:
        if not key:continue
        node=trie
        for token in key.split():node=node.setdefault(token,{})
        node[None]=key
    return trie


def scan(text,trie,clean):
    text=normal_text(text,clean);spans=list(WORD.finditer(text));tokens=[m.group() for m in spans];found=defaultdict(list)
    for start in range(len(tokens)):
        node=trie;end=start
        while end<len(tokens) and tokens[end] in node:
            node=node[tokens[end]];end+=1
            if None in node:
                a=spans[start].start();b=spans[end-1].end()
                hyphen=(a>0 and text[a-1] in HYPHENS) or (b<len(text) and text[b] in HYPHENS)
                found[node[None]].append({'token_start':start,'token_end_exclusive':end,'adjacent_hyphen':hyphen})
    return dict(found)


def read_tsv(path,snapshots):
    text=capture(path,snapshots).decode();reader=csv.DictReader(io.StringIO(text),delimiter='\t',strict=True);rows=list(reader)
    require(all(None not in x and all(v is not None for v in x.values()) for x in rows),'Malformed TSV: '+str(path))
    return rows


def target_names(catalog_dir,clean,snapshots,expected_appraisals):
    targets=defaultdict(set);counts={}
    def add(name,layer,id,field):
        key=norm(name,clean)
        if key:targets[key].add((layer,id,field,name))
    catalog=read_tsv(catalog_dir/'catalog.tsv',snapshots);counts['main_catalog_rows']=len(catalog)
    require(len(catalog)==14882,'Main catalog count differs');catids={x['catalog_id'] for x in catalog};require(len(catids)==len(catalog),'Duplicate catalog ID')
    for row in catalog:add(row['preferred_name'],'main_catalog',row['catalog_id'],'preferred_name')
    aliases=read_tsv(catalog_dir/'aliases.tsv',snapshots);counts['main_catalog_alias_cells']=len(aliases)
    for row in aliases:
        require(row['catalog_id'] in catids,'Orphan catalog alias')
        add(row['alias'],'main_catalog',row['catalog_id'],'alias:'+row['source']+':'+row['source_id'])
    supp=read_tsv(catalog_dir/'drugcentral-status-unverified.tsv',snapshots);counts['status_unverified_rows']=len(supp);require(len(supp)==475,'Supplement count differs')
    for row in supp:
        add(row['name'],'status_unverified',row['drugcentral_id'],'name')
        for alias in json.loads(row['aliases']):add(alias,'status_unverified',row['drugcentral_id'],'alias')
    food=read_tsv(catalog_dir/'food-substances/substances.tsv',snapshots);counts['FDA_food_source_rows']=len(food);require(len(food)==3971,'FDA source count differs')
    for row in food:
        add(row['substance_as_reported'],'FDA_food',row['source_row_id'],'substance_as_reported')
        add(row['retrieval_name'],'FDA_food',row['source_row_id'],'explicit_regulatory_suffix_removed_retrieval_name')
    decisions=json.loads(capture(catalog_dir/'decisions.json',snapshots))['decisions'];counts['appraisals']=len(decisions)
    require(len(decisions)==expected_appraisals and len({x['decision_id'] for x in decisions})==expected_appraisals,'Expected '+str(expected_appraisals)+' unique appraisals')
    for row in decisions:
        for name in row['match_aliases']:add(name,'curated_appraisal',row['decision_id'],'match_aliases')
    return targets,counts


def write_tsv(path,rows,fields):
    with path.open('w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fields,delimiter='\t',lineterminator='\n',extrasaction='raise');w.writeheader();w.writerows(rows)


def write_json(path,value):
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')


def publish(stage,out,snapshots,summary):
    recheck(snapshots)
    files={f.name:sha(f) for f in sorted(stage.iterdir()) if f.is_file()}
    summary['output_sha256']=files
    write_json(stage/'summary.json',summary)
    recheck(snapshots)
    run_id=digest((stage/'summary.json').read_bytes())[:20]
    runs=out/'runs';runs.mkdir(exist_ok=True)
    dest=runs/run_id
    if dest.exists():
        require(all(sha(dest/f.name)==sha(f) for f in stage.iterdir()),'Existing output run conflicts')
    else:os.replace(stage,dest)
    # Immutable run is complete before the sole public pointer is replaced.
    pointer=out/('.current-'+run_id)
    if pointer.is_symlink():pointer.unlink()
    require(not pointer.exists(),'Unexpected publication path')
    pointer.symlink_to(Path('runs')/run_id,target_is_directory=True)
    os.replace(pointer,out/'current')
    return dest


def build(args):
    snapshots={};capture(Path(__file__),snapshots)
    source=args.source_dir.resolve();catalog_dir=args.catalog_dir.resolve();out=args.out.resolve()
    groups,manifest,partitions=load_source(source,snapshots)
    clean=load_clean(catalog_dir,snapshots)
    # Source documentation is provenance, not proof of clinical availability.
    documentation={}
    for name in ['api-doc.html','api-guide.html','retrieval.json']:
        documentation[name]=digest(capture(source/name,snapshots))
    corpus=json.loads(checked(args.corpus,CORPUS_SHA,snapshots));records=corpus['results']
    require(len(records)==7362 and corpus['retrievedCount']==7362,'Corpus count differs')
    require(len({(x['source'],x['id']) for x in records})==7362,'Duplicate corpus source/ID')
    require(sum(bool(x.get('abstractText')) for x in records)==7191,'Saved abstract availability differs')
    targets,target_counts=target_names(catalog_dir,clean,snapshots,args.expected_appraisals)
    cells=[];normgroups=defaultdict(lambda:defaultdict(list));groupnames={};rawdup=0
    for row in groups:
        g=row['hit']['_source'];gid=g['groupId'];groupnames[gid]=g['groupName']
        seen=Counter();seen_norm=Counter()
        for role,index,value in [('group_name',-1,g['groupName'])]+[('synonym',i,s) for i,s in enumerate(g['synonyms'])]:
            key=norm(value,clean);cell_id='dsld:'+gid+(':'+role if index<0 else ':synonym:'+str(index))
            cells.append({'alias_cell_id':cell_id,'group_id':gid,'role':role,'synonym_index_zero_based':index if index>=0 else '', 'value_as_reported':value,'normalized_alias':key,'prior_identical_cells_in_group':seen[value],'prior_normalized_cells_in_group':seen_norm[key]})
            rawdup+=bool(seen[value]);seen[value]+=1;seen_norm[key]+=1;normgroups[key][gid].append(cell_id)
    require(len(cells)==85255,'Name+synonym cells differ')
    flags={key:ambiguity(key,gids) for key,gids in normgroups.items()}
    for cell in cells:cell['ambiguity_flags']=jdump(flags[cell['normalized_alias']])
    trie=make_trie(normgroups);mentions=[];source_rows=[];group_records=defaultdict(set);group_unflagged=defaultdict(set);alias_records=defaultdict(set);alias_occurrences=Counter();context_flags=Counter()
    for rec in records:
        sid=rec['source']+':'+str(rec['id']);title=rec.get('title') or '';abstract=rec.get('abstractText') or ''
        require(isinstance(title,str) and isinstance(abstract,str),'Nonstrings in source title/abstract')
        by_field={field:scan(value,trie,clean) for field,value in [('title',title),('abstract',abstract)]}
        keys=sorted(set(by_field['title'])|set(by_field['abstract']))
        direct_gene=bool(GENE.search(norm(title+' '+abstract,clean)))
        seen_groups=set()
        for key in keys:
            positions={field:hits[key] for field,hits in by_field.items() if key in hits}
            local_flags=set(flags[key])
            if any(x['adjacent_hyphen'] for seq in positions.values() for x in seq):local_flags.add('adjacent_hyphen_may_be_derivative_or_compound_phrase')
            alias_records[key].add(sid);alias_occurrences[key]+=sum(map(len,positions.values()))
            for gid,cellids in sorted(normgroups[key].items()):
                group_records[gid].add(sid);seen_groups.add(gid)
                if not local_flags:group_unflagged[gid].add(sid)
                for flag in local_flags:context_flags[flag]+=1
                mentions.append({'group_id':gid,'group_name_as_reported':groupnames[gid],'normalized_alias':key,'source_alias_cell_ids':jdump(cellids),'source_key':sid,'source':rec['source'],'id':str(rec['id']),'pmid':rec.get('pmid') or '', 'matched_fields':jdump(list(positions)),'token_positions':jdump(positions),'occurrences':sum(map(len,positions.values())),'direct_gene_in_title_or_abstract':direct_gene,'ambiguity_flags':jdump(sorted(local_flags)),'interpretation':'LEXICAL_DISCOVERY_CONTEXT_UNREVIEWED'})
        source_rows.append({'source_key':sid,'source':rec['source'],'id':str(rec['id']),'pmid':rec.get('pmid') or '', 'title_as_reported':title,'abstract_present':bool(abstract),'abstract_characters':len(abstract),'abstract_sha256':digest(abstract.encode()) if abstract else '', 'matched_group_count':len(seen_groups),'matched_normalized_alias_count':len(keys),'direct_gene_in_title_or_abstract':direct_gene,'source_metadata':jdump({k:rec[k] for k in METADATA_FIELDS if k in rec}),'reading_credit':'NONE_FROM_AUTOMATED_SCREEN'})
    crossrefs=[];group_targets=defaultdict(lambda:defaultdict(set));cross_layers=Counter()
    for key,gids in sorted(normgroups.items()):
        for gid,cellids in sorted(gids.items()):
            for layer,tid,field,name in sorted(targets.get(key,())):
                group_targets[gid][layer].add(tid);cross_layers[layer]+=1
                crossrefs.append({'group_id':gid,'normalized_alias':key,'source_alias_cell_ids':jdump(cellids),'target_layer':layer,'target_id':tid,'target_name_field':field,'target_name_as_reported':name,'identity_state':'EXACT_NORMALIZED_NAME_ONLY_NOT_CHEMICAL_EQUIVALENCE','ambiguity_flags':jdump(flags[key])})
    group_rows=[];categories=Counter();empty_syn=[];null_categories=0;case_categories=defaultdict(set)
    raw_group_names=Counter(r['hit']['_source']['groupName'] for r in groups)
    for r in groups:
        hit=r['hit'];g=hit['_source'];gid=g['groupId']
        for cat in g['category']:
            categories[jdump(cat)]+=1
            if cat is None:null_categories+=1
            else:case_categories[cat.casefold()].add(cat)
        if not g['synonyms']:empty_syn.append(gid)
        group_rows.append({'group_id':gid,'group_name_as_reported':g['groupName'],'source_partition':r['partition'],'source_file':r['source_file'],'source_hit_index_zero_based':r['source_hit_index'],'partition_row_index_zero_based':r['partition_row_index'],'source_url':r['source_url'],'source_index':hit['_index'],'source_type':hit['_type'],'source_score':hit['_score'],'category_as_reported':jdump(g['category']),'synonym_cell_count':len(g['synonyms']),'factsheets_as_reported':jdump(g['factsheets']),'nutrient_info_as_reported':jdump(g['nutrientInfo']),'main_catalog_name_candidates':jdump(sorted(group_targets[gid]['main_catalog'])),'status_unverified_name_candidates':jdump(sorted(group_targets[gid]['status_unverified'])),'FDA_food_name_candidates':jdump(sorted(group_targets[gid]['FDA_food'])),'appraisal_name_candidates':jdump(sorted(group_targets[gid]['curated_appraisal'])),'literature_source_record_count':len(group_records[gid]),'unflagged_lexical_source_record_count':len(group_unflagged[gid]),'identity_flags':jdump(['raw_group_name_shared_across_source_IDs'] if raw_group_names[g['groupName']]>1 else []),'evidence_state':'UNKNOWN_NOT_ASSESSED_BY_INVENTORY'})
    alias_rows=[]
    for key,gids in sorted(normgroups.items()):
        alias_rows.append({'normalized_alias':key,'group_ids':jdump(sorted(gids)),'group_count':len(gids),'source_cell_count':sum(map(len,gids.values())),'ambiguity_flags':jdump(flags[key]),'literature_source_record_count':len(alias_records[key]),'literal_occurrences':alias_occurrences[key]})
    policy={'normalization':'Corrected broad_hub_screen.clean (iterated HTML entity decoding and allowlisted tags only), Unicode NFKC, casefold, Unicode alphanumeric tokens excluding underscore; join tokens by single spaces. Punctuation including hyphens is a boundary, not part of tokens.','matching':'All overlapping normalized token sequences; fields searched separately; no substring within a token. All normalized aliases, including flagged terms, remain searchable except empty normalization.','hyphen_limit':'Component aliases may match within hyphenated chemical names; adjoining hyphen is explicitly flagged. No derivative-to-parent efficacy attribution is allowed.','name_crossrefs':'Full normalized-name equality only against catalog preferred/alias names, unverified-supplement names/aliases, FDA reported/retrieval names and curated match_aliases. FDA retrieval names only remove documented explicit regulatory suffixes. No fuzzy, CAS, structure, moiety, dose, mixture or product equivalence.','grouping_limit':'DSLD source-defined synonym relations are not verified chemical equivalence. Repeated cells retain indices; pattern indexing does not erase source cells.','manual_common_terms':sorted(COMMON),'other_flags':['short <=3 normalized characters','numeric-only','alias shared across groups','contextual adjacent hyphen'],'flag_limit':'Flags are triage aids, not complete ambiguity detection, intervention review, efficacy or exclusion decisions.','positions':'Zero-based token positions, end-exclusive, in the normalized decoded individual field. No raw abstract or snippet text written.','credit':'No abstract-reading, primary-paper, biological-validation or drug-appraisal credit from this screen.'}
    matched_source_ids={m['source_key'] for m in mentions};unflagged=[m for m in mentions if m['ambiguity_flags']=='[]']
    anomaly={'distinct_raw_group_names':len(raw_group_names),'raw_group_names_shared_by_multiple_IDs':sum(n>1 for n in raw_group_names.values()),'group_IDs_in_raw_duplicate_name_sets':sum(n for n in raw_group_names.values() if n>1),'duplicate_synonym_cells_within_group':sum(len(r['hit']['_source']['synonyms'])-len(set(r['hit']['_source']['synonyms'])) for r in groups),'category_cell_counts_preserving_JSON_type_and_case':dict(categories),'JSON_null_category_cells':null_categories,'literal_string_null_category_cells':categories.get('"null"',0),'category_case_variants':{k:sorted(v) for k,v in case_categories.items() if len(v)>1},'groups_with_empty_synonym_list':len(empty_syn),'empty_synonym_cells':[c['alias_cell_id'] for c in cells if c['role']=='synonym' and not c['value_as_reported'].strip()],'empty_normalized_alias_cells':sum(not c['normalized_alias'] for c in cells),'raw_duplicate_cells_including_group_name_repetitions':rawdup,'within_group_normalization_collision_cells':sum(bool(c['prior_normalized_cells_in_group']) for c in cells),'all_nutrient_info_empty':True,'documentation_version':'9.4.0 in saved API HTML; live before/after API version9.5.0','guide_visible_extraction':'Existing api-guide-visible.txt is empty; it was not treated as reviewed guide prose.','snapshot_limit':'27 API partitions are separately retrieved, not a transactional database snapshot; version equality does not prove no in-place changes.'}
    summary={'schema_version':1,'source':'NIH Office of Dietary Supplements DSLD ingredient groups','source_api_url':'https://api.ods.od.nih.gov/dsld/v9/ingredient-groups/','official_documentation_url':'https://dsldapi.od.nih.gov/','official_guide_url':'https://dsld.od.nih.gov/api-guide','retrieval_date':manifest['retrieval_date'],'api_version':manifest['api_version'],'partitions':partitions,'counts':{'source_partitions':27,'source_groups':6467,'distinct_source_group_ids':6467,'synonym_cells':78788,'group_name_cells':6467,'all_source_alias_cells':len(cells),'distinct_nonempty_normalized_aliases':sum(bool(k) for k in normgroups),'normalized_aliases_shared_across_groups':sum(len(x)>1 for x in normgroups.values()),'corpus_source_records':7362,'corpus_nonempty_saved_abstracts':7191,'corpus_without_saved_abstract':171,'groups_with_any_lexical_match':len([g for g in group_records.values() if g]),'groups_without_lexical_match':6467-len([g for g in group_records.values() if g]),'matched_normalized_aliases':sum(bool(x) for x in alias_records.values()),'matched_source_records':len(matched_source_ids),'group_alias_record_mentions':len(mentions),'unflagged_group_alias_record_mentions':len(unflagged),'unflagged_matched_source_records':len({m['source_key'] for m in unflagged}),'groups_with_unflagged_lexical_match':sum(bool(x) for x in group_unflagged.values()),'direct_gene_group_alias_record_mentions':sum(m['direct_gene_in_title_or_abstract'] for m in mentions),'unflagged_direct_gene_mentions':sum(m['direct_gene_in_title_or_abstract'] for m in unflagged),'name_crossref_rows':len(crossrefs),'groups_with_name_crossref_by_layer':{layer:sum(bool(group_targets[g][layer]) for g in groupnames) for layer in ['main_catalog','status_unverified','FDA_food','curated_appraisal']},'crossref_rows_by_layer':dict(cross_layers),'target_inventory':target_counts,'automatic_manual_abstract_reads':0,'automatic_primary_reads':0,'automatic_new_appraisals':0},'mention_ambiguity_flag_counts':dict(context_flags),'metadata_anomalies':anomaly,'documentation_sha256':documentation,'input_sha256':snapshots,'policy':'screen-policy.json','limitations':['Ingredient groups are not unique molecules, product labels, currently available products or verified medicines.','Source synonym grouping is not independently verified chemical identity.','Co-mentions can be endogenous, reagents, derivatives, harms, genetic contexts or background. No intervention/efficacy role is adjudicated here.','Unflagged means no listed lexical flag, not specific or beneficial by manual review.','Absence of a hit is unknown, not absence of biological activity.','Catalog/FDA/appraisal names are cross-references, not regulatory status, source equivalence or transfer of previous efficacy.','No current product, dose, supplier, availability or label inventory was retrieved.'],'publication':'All output files in one immutable run; current is switched atomically only after input hash recheck. Consumers should resolve current once and verify hashes; historical runs remain intact.'}
    require(not out.is_relative_to(source/'raw'),'Output cannot be inside raw sources')
    out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.stage-',dir=out) as tmp:
        stage=Path(tmp)
        for name,rows in [('groups.tsv',group_rows),('alias-cells.tsv',cells),('normalized-aliases.tsv',alias_rows),('mentions.tsv',mentions),('name-crossrefs.tsv',crossrefs),('source-records.tsv',source_rows)]:
            require(bool(rows),'Unexpected empty output: '+name);write_tsv(stage/name,rows,list(rows[0]))
        write_json(stage/'screen-policy.json',policy);write_json(stage/'metadata-anomalies.json',anomaly)
        dest=publish(stage,out,snapshots,summary)
    print(json.dumps({'published':str(dest),'current':str(out/'current'),'counts':summary['counts']},indent=2))
    return summary


def self_test(catalog_dir,out,source_dir):
    results=[];snapshots={};clean=load_clean(catalog_dir,snapshots)
    def test(name,fn):
        fn();results.append({'name':name,'passed':True})
    def equal(actual,expected):require(actual==expected,repr((actual,expected)))
    test('corrected_clean_preserves_inequalities_and_names',lambda:equal(norm('P <0.01 with GSK923295 > 2 &lt;i&gt;niacin&lt;/i&gt;',clean),'p 0 01 with gsk923295 2 niacin'))
    test('case_NFKC_whitespace_and_hyphen_equivalence',lambda:equal(norm('ＦＯＬＩＣ\tACID; N-acetyl',clean),'folic acid n acetyl'))
    def overlap():
        hits=scan('alpha lipoic acid and lipoic acid',make_trie(['alpha lipoic acid','lipoic acid','acid']),clean)
        equal({k:len(v) for k,v in hits.items()},{'alpha lipoic acid':1,'lipoic acid':2,'acid':2})
    test('all_nested_overlapping_aliases_retained',overlap)
    test('no_substring_inside_word',lambda:equal(scan('niacinamide is distinct from niacin',make_trie(['niacin']),clean)['niacin'][0]['token_start'],4))
    def hyphen():
        hits=scan('N-acetyl-L-cysteine and cysteine',make_trie(['l cysteine','cysteine']),clean)
        equal([x['adjacent_hyphen'] for x in hits['cysteine']],[True,False]);equal(hits['l cysteine'][0]['adjacent_hyphen'],True)
    test('hyphen_component_matches_flagged',hyphen)
    test('common_short_multigroup_not_suppressed',lambda:equal(ambiguity('cat',{'1','2'}),['short_alias_3_or_fewer_characters','alias_shared_across_DSLD_groups']))
    test('manual_common_term_flag',lambda:equal(ambiguity('protein',{'1'}),['manual_common_or_class_term']))
    def fixture_failure():
        try:parse_page({'total':{'value':2,'relation':'gte'},'max_score':1,'hits':[]},2,'A')
        except ValueError:return
        raise AssertionError('inexact totals accepted')
    test('inexact_source_total_rejected',fixture_failure)
    def partial_page():
        try:require_full_page([1],2,0,1000,'A')
        except ValueError:return
        raise AssertionError('partial exact-total page accepted')
    test('partial_page_with_exact_declared_total_rejected',partial_page)
    def source_cell_preservation():
        source={'groupId':'1','groupName':'Folic acid','category':[None,'null','Botanical','botanical'],'synonyms':['folic acid','folic acid','folic-acid',''],'factsheets':[],'nutrientInfo':[]}
        hit={'_id':'1','_index':'fixture','_type':'_doc','_score':1,'_source':source}
        result=parse_page({'total':{'value':1,'relation':'eq'},'max_score':1,'hits':[hit]},1,'A')
        equal(result[0][1]['_source'],source)
        equal(result[0][1]['_source']['synonyms'],['folic acid','folic acid','folic-acid',''])
    test('duplicate_empty_cells_and_category_null_type_case_preserved',source_cell_preservation)
    def conflict_id():
        ids={};register_group(ids,{'_id':'1','name':'one'})
        try:register_group(ids,{'_id':'1','name':'different'})
        except ValueError:return
        raise AssertionError('conflicting duplicate ID accepted')
    test('conflicting_duplicate_group_ID_rejected',conflict_id)
    def mismatch_id():
        hit={'_id':'1','_index':'fixture','_type':'_doc','_score':1,'_source':{'groupId':'2','groupName':'name','category':[],'synonyms':[],'factsheets':[],'nutrientInfo':[]}}
        try:parse_page({'total':{'value':1,'relation':'eq'},'max_score':1,'hits':[hit]},1,'A')
        except ValueError:return
        raise AssertionError('ID/groupId mismatch accepted')
    test('source_ID_and_groupID_mismatch_rejected',mismatch_id)
    out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.fixtures-',dir=out) as tmp:
        root=Path(tmp);input_file=root/'source';input_file.write_bytes(b'original');ss={str(input_file):sha(input_file)}
        def tamper():
            input_file.write_bytes(b'tampered')
            try:checked(input_file,ss[str(input_file)],{})
            except ValueError:return
            raise AssertionError('tampered source accepted')
        test('tampered_input_hash_rejected',tamper)
        input_file.write_bytes(b'original');dest=root/'published';dest.mkdir();(dest/'current').symlink_to('existing');(dest/'existing').mkdir();(dest/'existing'/'sentinel').write_text('keep')
        def late_change():
            stage=root/'stage';stage.mkdir();(stage/'output.tsv').write_text('staged')
            input_file.write_bytes(b'late mutation')
            try:publish(stage,dest,ss,{'test':True})
            except ValueError:
                equal(os.readlink(dest/'current'),'existing');equal((dest/'current'/'sentinel').read_text(),'keep');return
            raise AssertionError('late mutation published')
        test('changed_input_before_publish_preserves_prior_outputs',late_change)
        def missing_partition():
            fixture=root/'partial-source';fixture.mkdir()
            shutil.copytree(source_dir/'raw',fixture/'raw')
            shutil.copy2(source_dir/'retrieval-manifest.json',fixture/'retrieval-manifest.json')
            (fixture/'raw/groups-Z-00000.json').unlink()
            try:load_source(fixture,{})
            except (ValueError,FileNotFoundError):
                equal(os.readlink(dest/'current'),'existing');equal((dest/'current'/'sentinel').read_text(),'keep');return
            raise AssertionError('missing partition accepted')
        test('missing_partition_rejected_prior_output_preserved',missing_partition)
    return {'scope':'Deterministic data-handling fixtures only, no literature review or biological validation.','tests':results,'passed':len(results),'failed':0}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-dir',type=Path,default=DEFAULT_SOURCE)
    ap.add_argument('--catalog-dir',type=Path,default=DEFAULT_CATALOG)
    ap.add_argument('--corpus',type=Path,default=WORK/'literature-gene-unrestricted.json')
    ap.add_argument('--out',type=Path,default=DEFAULT_SOURCE/'derived')
    ap.add_argument('--expected-appraisals',type=int,default=168,help='Explicit ledger snapshot denominator; use 170 only after a separately authorized ledger update')
    ap.add_argument('--self-test',action='store_true')
    args=ap.parse_args()
    if args.self_test:
        print(json.dumps(self_test(args.catalog_dir.resolve(),args.out.resolve(),args.source_dir.resolve()),indent=2))
    else:build(args)

if __name__=='__main__':
    main()
