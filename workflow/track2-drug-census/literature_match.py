#!/usr/bin/env python3
"""Conservative title/abstract alias screen; co-mentions are not efficacy evidence."""
import argparse,csv,hashlib,html,json,re,unicodedata
from collections import defaultdict,Counter
from pathlib import Path
from catalog import write_tsv,write_json,stamp
ROOT=Path(__file__).resolve().parent
# These are chemistry/biology words, endogenous substances, assay terms, or ambiguous
# ordinary words. Keep their hits separate; do not silently classify them as drug experiments.
# Single-token entries and multi-token terms are explicit to avoid splitting compounds.
LOW_SPECIFICITY=set('water oxygen nitrogen carbon hydrogen calcium sodium potassium magnesium iron copper zinc lead mercury glucose dextrose fructose sucrose insulin albumin thyroxine cortisol adrenaline epinephrine norepinephrine dopamine serotonin melatonin glutamate glutamine glycine alanine leucine arginine histidine cysteine proline valine tryptophan histamine adenosine adenine guanine cytosine thymidine thymine uracil uridine inosine cytidine creatine creatinine cholesterol choline taurine urea ethanol methanol acetone formaldehyde fluoride chloride sulfate phosphate glycerol mannitol sorbitol starch cellulose heparin thrombin calcitonin erythropoietin interferon somatostatin estradiol progesterone testosterone vasopressin oxytocin saline serum plasma biotin riboflavin thiamine'.split())|{'hydrogen peroxide','dimethyl sulfoxide','nitric oxide','carbon dioxide','carbon monoxide','folic acid','ascorbic acid','vitamin a','vitamin b','vitamin c','vitamin d','vitamin e','retinoic acid','all trans retinoic acid','retinol','serine','l serine','threonine','l threonine','lysine','l lysine','tyrosine','l tyrosine','l cysteine','thyroid','estrogen','amino acids','annexin a5','galactosidase','glutathione','l glutamine','l arginine','l histidine','l leucine','l proline','l valine','l tryptophan','l tyrosine','l alanine','pyrimidine','triphosphate','arsenic','paraffin','alcohol','platinum','lactic acid','folate','fibrinogen','transferrin','ferritin','sars cov'}
def normalize(x):
    x=unicodedata.normalize('NFKC',html.unescape(re.sub('<[^>]+>',' ',x or ''))).casefold()
    return ' '.join(re.findall(r'[^\W_]+',x,flags=re.UNICODE))
def readtsv(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--corpus',type=Path,required=True);ap.add_argument('--title-abstract-corpus',type=Path,required=True);ap.add_argument('--out',type=Path,default=ROOT);a=ap.parse_args();out=a.out
    drugs=readtsv(out/'catalog.tsv');clinical_count=len(drugs);aliases=readtsv(out/'aliases.tsv');index=defaultdict(list)
    supplements=readtsv(out/'drugcentral-status-unverified.tsv')
    for d in supplements:
        gid='DC'+d['drugcentral_id'];drugs.append({'catalog_id':gid,'preferred_name':d['name'],'status_summary':'STATUS_UNVERIFIED_SUPPLEMENT','annotation_status':'supplemental_source_record_status_unverified'})
        for name in set(json.loads(d['aliases'])+[d['name']]):
            aliases.append({'catalog_id':gid,'source':'DrugCentral','source_id':d['drugcentral_id'],'alias_type':'preferred' if name==d['name'] else 'name_or_synonym','alias':name})
    for row in aliases:index[normalize(row['alias'])].append(row)
    trie={};policies=[]
    for key,rows in sorted(index.items()):
        reasons=[];ids={r['catalog_id'] for r in rows};letters=sum(c.isalpha() for c in key)
        if len(key.replace(' ',''))<5 or letters<3:reasons.append('short_or_code_ambiguity')
        if key in LOW_SPECIFICITY:reasons.append('common_biological_chemical_or_assay_term')
        if len(ids)>1:reasons.append('multiple_catalog_entities')
        named_types={r['alias_type'] for r in rows}
        # Brand/untyped synonyms can be ordinary words; preserve hits for manual review.
        if not (named_types-{'OTHER','TRADE_NAME','name_or_synonym'}):reasons.append('trade_name_or_untyped_synonym_only')
        policy='ambiguous_manual_review' if reasons else 'specific_alias_comention_only'
        policies.append({'normalized_alias':key,'catalog_ids':sorted(ids),'policy':policy,'reasons':reasons,'original_aliases':sorted({r['alias'] for r in rows})})
        if not key or letters<2 or len(key.replace(' ',''))<3:continue
        node=trie
        for token in key.split():node=node.setdefault(token,{})
        node['']=policies[-1]
    raw=a.corpus.read_bytes();corpus=json.loads(raw);ta=json.loads(a.title_abstract_corpus.read_bytes())
    ta_ids={(r['source'],r['id']) for r in ta['results']};by_drug=defaultdict(set);amb_by_drug=defaultdict(set);mentions=[];ambiguous=[]
    names={r['catalog_id']:r['preferred_name'] for r in drugs};preferred={r['catalog_id']:normalize(r['preferred_name']) for r in drugs}
    abstract_present=0;direct_gene_count=0
    for record in corpus['results']:
        rid=(record['source'],record['id']);title=normalize(record.get('title',''));abstract=normalize(record.get('abstractText',''))
        abstract_present+=bool(abstract);text=title+' '+abstract;tokens=text.split();found={};direct_gene=bool(re.search(r'\b(bub1b|bubr1|bub r1|mosaic variegated aneuploidy)\b',text));direct_gene_count+=direct_gene
        for i in range(len(tokens)):
            node=trie;j=i
            while j<len(tokens) and tokens[j] in node:
                node=node[tokens[j]];j+=1
                if '' in node:found[node['']['normalized_alias']]=node['']
        hits=defaultdict(list);ambhits=defaultdict(list)
        for key,policy in found.items():
            for gid in policy['catalog_ids']:
                (ambhits if policy['policy']=='ambiguous_manual_review' else hits)[gid].append(key)
        for hitdict,dest,bucket in [(hits,mentions,by_drug),(ambhits,ambiguous,amb_by_drug)]:
            for gid,keys in hitdict.items():
                bucket[gid].add(rid)
                dest.append({'catalog_id':gid,'drug_name':names[gid],'source':record['source'],'record_id':record['id'],'pmid':record.get('pmid',''),'pmcid':record.get('pmcid',''),'doi':record.get('doi',''),'title':html.unescape(re.sub('<[^>]+>','',record.get('title',''))),'matched_normalized_aliases':sorted(keys),'preferred_name_matched':preferred[gid] in keys,'direct_gene_name_in_title_or_abstract':direct_gene,'member_of_gene_title_abstract_query':rid in ta_ids,'abstract_available':bool(abstract),'evidence_state':'CO_MENTION_ONLY_NOT_EFFICACY_OR_RESCUE','source_url':'https://europepmc.org/article/'+record['source']+'/'+record['id']})
    rows=[]
    for d in drugs:
        gid=d['catalog_id'];rows.append({'catalog_id':gid,'preferred_name':d['preferred_name'],'status_summary':d['status_summary'],'specific_alias_paper_count':len(by_drug[gid]),'ambiguous_alias_paper_count':len(amb_by_drug[gid]),'literature_evidence_state':'CO_MENTION_REQUIRES_MANUAL_REVIEW' if by_drug[gid] else 'UNKNOWN_NO_SPECIFIC_ALIAS_TITLE_ABSTRACT_MATCH','annotation_status':d['annotation_status']})
    write_tsv(out/'literature-alias-policy.tsv',policies,list(policies[0]));write_tsv(out/'literature-mentions.tsv',mentions,list(mentions[0]) if mentions else ['catalog_id']);write_tsv(out/'literature-ambiguous-mentions.tsv',ambiguous,list(ambiguous[0]) if ambiguous else ['catalog_id']);write_tsv(out/'literature-drug-coverage.tsv',rows,list(rows[0]))
    summary={'built_utc':stamp(),'corpus_query':corpus['query'],'corpus_hit_count':corpus['hitCount'],'corpus_records_scanned':len(corpus['results']),'records_with_abstracts':abstract_present,'records_with_gene_literal_in_title_abstract':direct_gene_count,'gene_title_abstract_query_count':ta['hitCount'],'corpus_sha256':hashlib.sha256(raw).hexdigest(),'title_abstract_corpus_sha256':hashlib.sha256(a.title_abstract_corpus.read_bytes()).hexdigest(),'clinical_or_approved_catalog_entities_screened':clinical_count,'status_unverified_supplement_entities_screened':len(supplements),'catalog_entities_screened':len(drugs),'normalized_aliases':len(policies),'alias_policy_counts':dict(Counter(x['policy'] for x in policies)),'entities_with_specific_alias_comention':sum(bool(by_drug[d['catalog_id']]) for d in drugs),'entities_without_specific_alias_comention':sum(not by_drug[d['catalog_id']] for d in drugs),'clinical_catalog_entities_with_comention':sum(bool(by_drug[d['catalog_id']]) for d in drugs[:clinical_count]),'status_unverified_entities_with_comention':sum(bool(by_drug[d['catalog_id']]) for d in drugs[clinical_count:]),'specific_alias_comention_rows':len(mentions),'specific_alias_direct_gene_comention_rows':sum(r['direct_gene_name_in_title_or_abstract'] for r in mentions),'ambiguous_comention_rows':len(ambiguous),'interpretation':'Conservative lexical retrieval screen only. A match can refer to a control, assay reagent, harmful exposure, other disorder, endogenous molecule, cited background, or drug combination. No match is UNKNOWN, not inefficacy. Broad-query records can match the gene only in full text/references; the explicit direct-gene flag distinguishes title/abstract co-mentions. Text normalization permits punctuation/hyphen differences; ambiguous aliases are separated. Exact studies and intervention direction require manual appraisal.','reproduce_command':'python3 literature_match.py --corpus /path/literature-gene-unrestricted.json --title-abstract-corpus /path/literature-gene-title-abstract.json'}
    write_json(out/'literature-coverage.json',summary)
    if (out/'decisions.json').exists():join_decisions(out,drugs,index)
    print(json.dumps(summary,indent=2))
def join_decisions(out, drugs, index):
    payload=json.loads((out/'decisions.json').read_text()); matches=[];candidate_decisions=defaultdict(list)
    for dec in payload['decisions']:
        hits={}
        for alias in dec['match_aliases']:
            ids=sorted({r['catalog_id'] for r in index.get(normalize(alias),[])})
            if ids:hits[alias]=ids
        candidates=sorted({gid for ids in hits.values() for gid in ids})
        collision=any(len(ids)>1 for ids in hits.values()) or len(candidates)>1
        for gid in candidates:candidate_decisions[gid].append((dec,collision))
        matches.append({'decision_id':dec['decision_id'],'route':dec['route'],'catalog_ids':candidates,'matched_aliases':hits,'join_status':'AMBIGUOUS_MULTIPLE_ENTITIES' if collision else 'EXACT_ALIAS_SINGLE_ENTITY' if candidates else 'NO_CATALOG_MATCH','disposition':dec['disposition']})
    rows=[]
    for drug in drugs:
        vals=candidate_decisions[drug['catalog_id']];amb=any(v[1] for v in vals) or len(vals)>1
        rows.append({'catalog_id':drug['catalog_id'],'preferred_name':drug['preferred_name'],'catalog_status':drug['status_summary'],'decision_ids':[v[0]['decision_id'] for v in vals],'appraisal_join_status':'AMBIGUOUS_REQUIRES_IDENTITY_REVIEW' if amb else 'MANUAL_DECISION_EXACT_ALIAS_JOIN' if vals else 'UNKNOWN_NOT_IN_CURATED_APPRAISALS','manual_dispositions':[v[0]['disposition'] for v in vals] if not amb else [],'bub1b_variant_efficacy':'NOT_ESTABLISHED_BY_CENSUS','rationale_source':'decisions.json; source-level appraisal required'})
    write_tsv(out/'catalog-appraisals.tsv',rows,list(rows[0]));write_tsv(out/'decision-alias-joins.tsv',matches,list(matches[0]))

if __name__=='__main__':main()
