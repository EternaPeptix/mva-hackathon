#!/usr/bin/env python3
"""Reconcile the final outer-abstract selection, not biological efficacy or actual reading.

This selection takes every remaining saved outer abstract, so the selection rule
is "all of them in hash order" rather than a fixed head count. Passing this check
means the 3,598 reviewed records are exactly the records that remained, replayed
from the original retrievals; it does not mean any drug works.
"""
import argparse,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sha=lambda b:hashlib.sha256(b).hexdigest()
key=lambda r:(r['source'],r.get('record_id',r.get('id')))
PRIOR=['outer-corpus-review.json','outer-next120-review.json','outer-ranked-completion-review.json',
       'outer-unfiltered-review.json','outer-unfiltered-next576-review.json',
       'outer-unfiltered-next256-review.json','outer-unfiltered-next512-review.json']
def check(root=ROOT,corpus=None,gene_corpus=None):
    read=lambda n:json.loads((root/n).read_text())
    s=read('outer-unfiltered-final-selection.json');v=read('outer-unfiltered-final-review.json')
    assert sha((root/'outer-unfiltered-final-selection.json').read_bytes())==v['selection_sha256']
    assert s['source_sha256']==v['source_sha256']
    assert set(s['prior_review_sha256'])==set(PRIOR)
    previous=[]
    for name,digest in s['prior_review_sha256'].items():
        assert sha((root/name).read_bytes())==digest,name
        previous+=read(name)['records']
    prev={key(r) for r in previous};gene={key(r) for r in read('abstract-review-notes.json')['records']}
    assert len(previous)==len(prev)==2161 and len(gene)==1460
    selected=s['records'];reviewed=v['records'];seen={key(r) for r in selected}
    assert len(selected)==len(reviewed)==len(seen)==s['selected_count']==v['reviewed_count']==3598
    assert not seen&(prev|gene)
    for i,(a,b) in enumerate(zip(selected,reviewed),1729):
        assert all(a[k]==b[k] for k in a),key(a)
        assert a['unfiltered_index']==i
        assert a['order_key']==sha((s['seed']+'\0'+a['source']+'\0'+a['record_id']).encode())
        assert b['review_depth']=='complete_saved_abstract'
        assert b['note'] and b['uncertainty'] and b['model'] and b['intervention_type']
        assert b['outcome_direction'] and b['source_abstract_quality']
        assert isinstance(b['requires_primary_followup'],bool)
        assert isinstance(b['administered_agents'],list)
        assert all(x.get('name') and x.get('role') for x in b['administered_agents'])
        assert not {'abstractText','abstract','body_text'}&b.keys()
    assert [a['order_key'] for a in selected]==sorted(a['order_key'] for a in selected)
    # the queue resumes strictly after the last record of the preceding tranche
    last=max(r['order_key'] for r in read('outer-unfiltered-next512-review.json')['records'])
    assert selected[0]['order_key']>last
    assert sum(a['requires_primary_followup'] for a in reviewed)==v['followup_count']
    assert sum(c['reviewed_count'] for c in v['chunks'])==3598 and len(v['chunks'])==57
    assert sum(c['followup_count'] for c in v['chunks'])==v['followup_count']
    assert s['remaining_outer_records']==3741 and s['remaining_available_abstracts']==3598
    assert s['remaining_missing_abstracts']==143
    if (corpus is None)!=(gene_corpus is None):raise ValueError('Supply both source retrievals')
    replayed=False
    if corpus is not None:
        assert sha(corpus.read_bytes())==s['source_sha256']
        assert sha(gene_corpus.read_bytes())==s['gene_source_sha256']
        raw=json.loads(corpus.read_text())['results'];inside=json.loads(gene_corpus.read_text())['results']
        assert {key(a) for a in inside}==gene
        remaining=[a for a in raw if key(a) not in prev|gene]
        available=[a for a in remaining if a.get('abstractText')]
        assert len(remaining)==3741 and len(available)==3598
        order=lambda a:sha((s['seed']+'\0'+a['source']+'\0'+a['id']).encode())
        ordered=sorted(available,key=order)
        assert len(ordered)==len(reviewed)
        for a,b in zip(ordered,reviewed):
            assert key(a)==key(b) and a['title']==b['title']
            assert sha(a['abstractText'].encode())==b['abstract_sha256']
        replayed=True
    return dict(result='PASS',additional_saved_abstract_entries=3598,
                combined_outer_source_records=5759,original_saved_outer_abstracts_unread=0,
                outer_records_without_saved_abstract=143,
                followup_flags=v['followup_count'],source_selection_replayed=replayed)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT)
    p.add_argument('--corpus',type=Path);p.add_argument('--gene-corpus',type=Path)
    a=p.parse_args();print(json.dumps(check(a.root,a.corpus,a.gene_corpus),indent=2))
