import gzip,tempfile,unittest
from pathlib import Path
import independent_vcf_reader_v2 as v
import verify_saved_residuals as checker
class Controls(unittest.TestCase):
 def test_full_literal_reader_and_selected_original_ordinal(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'fiction.vcf.gz'
   headers=['##fileformat=VCFv4.3']
   for tag in v.ANNOTATIONS:
    kind='Flag' if tag=='TRF' else 'Float' if tag in ('LCR','TRFentropy','TRFdiff','TRFovl','TRFcopies','TRFsim') else 'Integer' if tag in ('TRFstart','TRFend','TRFperiod','TRFscore','RM_score') else 'String'
    headers.append(f'##INFO=<ID={tag},Number={0 if kind=="Flag" else 1},Type={kind},Description="Fictional, annotation">')
   headers.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002')
   rows=['chr20\t1\tID\tA\tACC\t.\tPASS\tTRF;LCR=0;RM_score=0;RM_repeat=.;RM_clsfam=LINE/L1\tGT:BD\t0|1:TP',
         'chr20\t1\tID\tA\tAGG\t.\tPASS\tLCR=.\tGT:BD\t1/0:FN']
   with gzip.open(p,'wt') as f:f.write('\n'.join(headers+rows)+'\n')
   got=v.read_vcf(p);self.assertEqual(got['total_records'],2)
   a,b=got['rows'];self.assertNotEqual(v.identity(a),v.identity(b));self.assertEqual(a['record_literal'],rows[0]);self.assertEqual(a['format_values']['GT'],'0|1')
   self.assertEqual(a['annotations']['LCR']['value'],0.0);self.assertEqual(a['annotations']['RM_score']['value'],0);self.assertIsNone(b['annotations']['LCR']['value'])
   one=v.read_vcf(p,{v.identity(b)});self.assertEqual((len(one['rows']),one['rows'][0]['ordinal'],one['total_records']),(1,2,2))
 def test_four_outcomes_and_missing_zero_aggregation(self):
  def field(present,value):return {'present':present,'value':value,'raw':None}
  rows=[]
  for i,out in enumerate(checker.OUTCOMES):
   a={k:field(False,None) for k in v.ANNOTATIONS};a['LCR']=field(i!=0,None if i in(0,1) else float(i-2));a['RM_clsfam']=field(i==3,'LINE/L1' if i==3 else None)
   a['TRF']=field(i==2,True if i==2 else None)
   rows.append({'type':'INS','size_bin':'50-99','outcome':out,'annotations':a})
  tables=checker.independent_tables(rows);self.assertEqual(tables['total']['outcomes'],dict.fromkeys(checker.OUTCOMES,1));self.assertEqual(tables['total']['shared_miss_fraction'],.25)
  for name in ['class_size','TRF_RM_by_class_size','RM_class_family_by_class_size','REMAP_by_class_size']:self.assertEqual(sum(x['truth_records'] for x in tables[name]),4)
  cells=[r for r in tables['raw_LCR_by_class_size_and_outcome'] if r['type']=='INS' and r['size_bin']=='50-99']
  self.assertEqual([r['median'] for r in cells],[None,None,0.0,1.0]);self.assertEqual([r['not_encoded'] for r in cells],[1,0,0,0]);self.assertEqual([r['explicit_missing'] for r in cells],[0,1,0,0])
 def test_genotype_allele_content_not_phase_claim(self):
  a={'format_values':{'GT':'0|1'}};b={'format_values':{'GT':'1/0'}}
  self.assertEqual(checker.gt(a),checker.gt(b));self.assertNotEqual(a['format_values']['GT'],b['format_values']['GT'])
if __name__=='__main__':unittest.main()
