import unittest
from independent_vcf_reader_v2 import annotation,declaration,identity,info_fields
class ReaderControls(unittest.TestCase):
 def test_quoted_header_commas(self):
  d=declaration('##INFO=<ID=LCR,Number=1,Type=Float,Description="Entropy, not binary">')
  self.assertEqual(d['Description'],'"Entropy, not binary"')
 def test_zero_versus_absent_and_missing(self):
  d={'LCR':{'Number':'1','Type':'Float'}}
  self.assertEqual(annotation({'LCR':'0'},'LCR',d)['value'],0.)
  self.assertFalse(annotation({},'LCR',d)['present'])
  self.assertEqual(annotation({'LCR':'.'},'LCR',d)['value'],None)
  self.assertTrue(annotation({'LCR':'.'},'LCR',d)['present'])
 def test_flag_is_not_numeric(self):
  d={'TRF':{'Number':'0','Type':'Flag'}}
  self.assertEqual(annotation({'TRF':None},'TRF',d)['value'],True)
  with self.assertRaises(ValueError):annotation({'TRF':'0'},'TRF',d)
 def test_same_length_different_ALT(self):
  a=dict(chrom='chrT',pos=1,id='same',ref='A',alts=['ACC']);b={**a,'alts':['AGG']}
  self.assertNotEqual(identity(a),identity(b))
 def test_duplicate_key_refused(self):
  with self.assertRaises(ValueError):info_fields('TRF;TRF')
 def test_nonfinite_refused(self):
  with self.assertRaises(ValueError):annotation({'LCR':'NaN'},'LCR',{'LCR':{'Number':'1','Type':'Float'}})
if __name__=='__main__':unittest.main()
