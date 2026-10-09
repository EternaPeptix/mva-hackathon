#!/usr/bin/env python3
"""Read-only exact alias cross-references and complete saved title/abstract screen."""
import argparse,collections,csv,datetime,hashlib,json,re,unicodedata
from html.parser import HTMLParser
from pathlib import Path

def norm(x):return ' '.join(unicodedata.normalize('NFKC',x or '').casefold().split())
class TextOnly(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.parts=[]
    def handle_data(self,data):self.parts.append(data)
    def handle_starttag(self,tag,attrs):
        if tag in {'p','br','div','h1','h2','h3','h4','h5','h6','li'}:self.parts.append(' ')
    def handle_endtag(self,tag):
        if tag in {'p','div','h1','h2','h3','h4','h5','h6','li'}:self.parts.append(' ')
def plain_text(value):
    parser=TextOnly();parser.feed(value or '');parser.close()
    return ''.join(parser.parts)
def textnorm(x):
    x=plain_text(x)
    x=re.sub('[\u2010-\u2015\u2212]','-',x)
    return norm(x)
def readtsv(p):
    with p.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write_tsv(p,rows,fields):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v,ensure_ascii=False,separators=(',',':')) if isinstance(v,(list,dict,bool)) else v for k,v in r.items()})
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--catalog-dir',type=Path,required=True);p.add_argument('--corpus',type=Path,required=True);p.add_argument('--aliases',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
    d=a.catalog_dir;out=a.output_dir;out.mkdir(parents=True,exist_ok=True)
    idx={k:collections.defaultdict(set) for k in ['main','status_unverified_supplement','fda2026','broad_preferred_names']}
    for r in readtsv(d/'catalog.tsv'):idx['main'][norm(r['preferred_name'])].add(r['catalog_id'])
    for r in readtsv(d/'aliases.tsv'):idx['main'][norm(r['alias'])].add(r['catalog_id'])
    for r in readtsv(d/'drugcentral-status-unverified.tsv'):
        for v in [r['name']]+json.loads(r['aliases']):idx['status_unverified_supplement'][norm(v)].add('DC:'+r['drugcentral_id'])
    for r in json.loads((d/'approval-overrides.json').read_text())['records']:
        for v in [r['brand'],r['active_ingredient']]+r['search_aliases']:idx['fda2026'][norm(v)].add(r['brand'])
    for r in readtsv(d/'broad-hub-screen.tsv'):idx['broad_preferred_names'][norm(r['preferred_source_name'])].add(r['hub_row_id'])
    corpus=json.loads(a.corpus.read_bytes());records=corpus['results'];assert corpus['hitCount']==corpus['retrievedCount']==len(records),'Incomplete saved corpus'
    corpus_text=[(r,textnorm(r.get('title',''))+' '+textnorm(r.get('abstractText',''))) for r in records]
    aliases=json.loads(a.aliases.read_text());matches=[];mentions=[]
    for alias in aliases:
        key=norm(alias['alias']);row={**alias,'normalized_alias':key,'catalog_cross_references':{layer:sorted(index[key]) for layer,index in idx.items() if key in index}}
        pat=re.compile(r'(?<!\w)'+re.escape(textnorm(alias['alias']))+r'(?!\w)');hits=[]
        for r,text in corpus_text:
            if pat.search(text):
                hit={'medicine':alias['medicine'],'matched_alias':alias['alias'],'alias_kind':alias['alias_kind'],'source':r['source'],'record_id':r['id'],'pmid':r.get('pmid',''),'doi':r.get('doi',''),'title':plain_text(r.get('title','')),'abstract_available':bool(r.get('abstractText')),'screen_state':'LEXICAL_CANDIDATE_ONLY_NOT_EFFICACY'};mentions.append(hit);hits.append({'source':r['source'],'record_id':r['id']})
        row['saved_corpus_record_matches']=hits;row['identity_merge_performed']=False;matches.append(row)
    fields=list(matches[0]);write_tsv(out/'ema-identity-match-rows.tsv',matches,fields)
    mention_fields=['medicine','matched_alias','alias_kind','source','record_id','pmid','doi','title','abstract_available','screen_state'];write_tsv(out/'ema-identity-mentions.tsv',mentions,mention_fields)
    (out/'ema-identity-match-rows.json').write_text(json.dumps(matches,ensure_ascii=False,indent=2)+'\n')
    meds={m:{'aliases_screened':sum(r['medicine']==m for r in matches),'catalog_match_rows':sum(r['medicine']==m and bool(r['catalog_cross_references']) for r in matches),'saved_corpus_distinct_record_matches':len({(r['source'],r['record_id']) for r in mentions if r['medicine']==m})} for m in sorted({r['medicine'] for r in aliases})}
    summary={'schema_version':1,'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'medicines':meds,'catalog_normalization':'NFKC, casefold, whitespace collapse; explicit official names and separately labeled typography variants only. No fuzzy, salt, combination or platform merge.','corpus_normalization':'HTML text extraction/unescape, Unicode dash-to-hyphen, NFKC/casefold/whitespace collapse, word-boundary literal alias matching. Titles and saved abstracts only; not full texts.','corpus_query':corpus['query'],'corpus_records':len(records),'corpus_records_with_abstract':sum(bool(r.get('abstractText')) for r in records),'corpus_sha256':sha(a.corpus),'alias_input_sha256':sha(a.aliases),'catalog_input_sha256':{f:sha(d/f) for f in ['catalog.tsv','aliases.tsv','drugcentral-status-unverified.tsv','approval-overrides.json','broad-hub-screen.tsv']},'limits':'Zero matches do not demonstrate null efficacy, prove a drug absent under every possible alias, or certify worldwide uniqueness. Historical 180-product mismatch snapshot and 15,357 entity denominator unchanged.'}
    (out/'ema-identity-screen-summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
