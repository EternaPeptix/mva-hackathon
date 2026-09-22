#!/usr/bin/env python3
"""Reconcile fixed next-256 coverage, not biological efficacy or actual reading."""
import argparse,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sha=lambda b:hashlib.sha256(b).hexdigest()
key=lambda r:(r['source'],r.get('record_id',r.get('id')))
def check(root=ROOT,corpus=None,gene_corpus=None):
    read=lambda n:json.loads((root/n).read_text())
    s=read('outer-unfiltered-next256-selection.json');v=read('outer-unfiltered-next256-review.json')
    assert sha((root/'outer-unfiltered-next256-selection.json').read_bytes())==v['selection_sha256']
    assert s['source_sha256']==v['source_sha256']
    previous=[]
    for name,digest in s['prior_review_sha256'].items():
        assert sha((root/name).read_bytes())==digest,name
        previous+=read(name)['records']
    prev={key(r) for r in previous};gene={key(r) for r in read('abstract-review-notes.json')['records']}
    assert len(previous)==len(prev)==1393 and len(gene)==1460
    selected=s['records'];reviewed=v['records'];seen={key(r) for r in selected}
    assert len(selected)==len(reviewed)==len(seen)==s['selected_count']==v['reviewed_count']==256
    assert not seen&(prev|gene)
    for i,(a,b) in enumerate(zip(selected,reviewed),961):
        assert all(a[k]==b[k] for k in a),key(a)
        assert a['unfiltered_index']==i
        assert a['order_key']==sha((s['seed']+'\0'+a['source']+'\0'+a['record_id']).encode())
        assert b['review_depth']=='complete_saved_abstract'
        assert b['note'] and b['uncertainty'] and b['model'] and b['intervention_type']
        assert isinstance(b['requires_primary_followup'],bool)
        assert isinstance(b['administered_agents'],list)
        assert all(x.get('name') and x.get('role') for x in b['administered_agents'])
        assert not {'abstractText','abstract','body_text'}&b.keys()
    assert [a['order_key'] for a in selected]==sorted(a['order_key'] for a in selected)
    assert sum(a['requires_primary_followup'] for a in reviewed)==v['followup_count']
    assert s['remaining_outer_records']==4509 and s['remaining_available_abstracts']==4366 and s['remaining_missing_abstracts']==143
    if (corpus is None)!=(gene_corpus is None):raise ValueError('Supply both source retrievals')
    if corpus is not None:
        assert sha(corpus.read_bytes())==s['source_sha256']
        assert sha(gene_corpus.read_bytes())==s['gene_source_sha256']
        raw=json.loads(corpus.read_text())['results'];inside=json.loads(gene_corpus.read_text())['results']
        assert {key(a) for a in inside}==gene
        available=[a for a in raw if a.get('abstractText') and key(a) not in prev|gene]
        assert len(available)==4366
        order=lambda a:sha((s['seed']+'\0'+a['source']+'\0'+a['id']).encode())
        ordered=sorted(available,key=order)[:256]
        for a,b in zip(ordered,reviewed):
            assert key(a)==key(b) and a['title']==b['title']
            assert sha(a['abstractText'].encode())==b['abstract_sha256']
    return dict(result='PASS',additional_saved_abstract_entries=256,combined_outer_source_records=1649,original_saved_outer_abstracts_unread=4110,followup_flags=v['followup_count'],source_selection_replayed=corpus is not None)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT);p.add_argument('--corpus',type=Path);p.add_argument('--gene-corpus',type=Path);a=p.parse_args();print(json.dumps(check(a.root,a.corpus,a.gene_corpus),indent=2))
