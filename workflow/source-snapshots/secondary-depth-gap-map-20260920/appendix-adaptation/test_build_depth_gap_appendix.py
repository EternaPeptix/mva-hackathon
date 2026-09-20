"""Fictional controls for CLI/count adaptation; no case inputs are opened."""
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SOURCE=Path(__file__).with_name('build_depth_gap_appendix.py')
MODES=('MQ20_BQ20','MQ60_BQ20')


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


class AppendixAdaptationControls(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.base=Path(self.temp.name).resolve()
        self.input=self.base/'fictional-input';self.input.mkdir();self.out=self.base/'fictional-output'
        metrics={MODES[0]:dict(CDS_union_bases=4,below20_bases=2,zero_depth_bases=2),
                 MODES[1]:dict(CDS_union_bases=4,below20_bases=2,zero_depth_bases=0)}
        rows=[dict(chrom='1',start=i,end=i+1,gene_ids=['gene:fictional'],
                   transcript_ids=['transcript:one','transcript:two'],depth_modes=list(MODES),categories=['zero','1-19']) for i in (0,1)]
        with gzip.open(self.input/'gap-intervals.jsonl.gz','wt') as f:
            for row in rows:f.write(json.dumps(row)+'\n')
        summary=dict(status='complete',gap_intervals=2,paired_categories=[],MQ60_depth_greater_than_MQ20_bases=0,
                     panel_union=metrics,genes=[dict(gene_id='gene:fictional',name='FICTIONAL',modes=metrics)],
                     transcripts=[dict(transcript_id=t,gene_id='gene:fictional',modes=metrics,cds_rows=[]) for t in rows[0]['transcript_ids']])
        (self.input/'summary.json').write_text(json.dumps(summary))
        self.write_status()
    def tearDown(self):self.temp.cleanup()
    def write_status(self):
        (self.input/'status.json').write_text(json.dumps(dict(status='complete',outputs={'gap-intervals.jsonl.gz':{'sha256':sha(self.input/'gap-intervals.jsonl.gz')}})))
    def run_builder(self,digest=None):
        return subprocess.run([sys.executable,str(SOURCE),'--input-dir',str(self.input),'--output-dir',str(self.out),
            '--expected-summary-sha256',digest or sha(self.input/'summary.json'),'--result-link','fictional-result.md'],capture_output=True,text=True)
    def test_dynamic_counts_and_adjacent_merge(self):
        r=self.run_builder();self.assertEqual(r.returncode,0,r.stderr)
        text=(self.out/'GAP-LOCATIONS.md').read_text()
        self.assertIn('[All 2 transcript summaries',text);self.assertIn('[Original 2 membership/category rows]',text)
        j=json.loads((self.out/'gap-locations.json').read_text());rows=j['genes'][0]['intervals'][MODES[0]]['below20']
        self.assertEqual([(x['start'],x['end'],x['atomic_row_numbers']) for x in rows],[(0,2,[1,2])])
    def test_input_output_links_and_source_preservation(self):
        before={p.name:sha(p) for p in self.input.iterdir()};r=self.run_builder();self.assertEqual(r.returncode,0,r.stderr)
        self.assertIn('](../fictional-input/gap-intervals.jsonl.gz)',(self.out/'GAP-LOCATIONS.md').read_text())
        self.assertIn('](fictional-result.md)',(self.out/'GAP-LOCATIONS.md').read_text())
        self.assertEqual(before,{p.name:sha(p) for p in self.input.iterdir()})
    def test_wrong_summary_pin_refused(self):
        self.assertNotEqual(self.run_builder('0'*64).returncode,0);self.assertFalse(self.out.exists())
    def test_wrong_interval_pin_refused(self):
        (self.input/'status.json').write_text(json.dumps({'outputs':{'gap-intervals.jsonl.gz':{'sha256':'0'*64}}}))
        self.assertNotEqual(self.run_builder().returncode,0);self.assertFalse(self.out.exists())
    def test_existing_output_preserved(self):
        self.out.mkdir();p=self.out/'keep.txt';p.write_text('keep')
        self.assertNotEqual(self.run_builder().returncode,0);self.assertEqual(p.read_text(),'keep')
    def test_empty_gap_inputs_render_zero_rows(self):
        with gzip.open(self.input/'gap-intervals.jsonl.gz','wt') as f:pass
        p=self.input/'summary.json';x=json.loads(p.read_text());x['gap_intervals']=0
        for entry in [dict(modes=x['panel_union'])]+x['genes']+x['transcripts']:
            for m in entry['modes'].values():m.update(below20_bases=0,zero_depth_bases=0)
        p.write_text(json.dumps(x));self.write_status();r=self.run_builder();self.assertEqual(r.returncode,0,r.stderr)
        self.assertIn('[Original 0 membership/category rows]',(self.out/'GAP-LOCATIONS.md').read_text())


if __name__=='__main__':unittest.main()
