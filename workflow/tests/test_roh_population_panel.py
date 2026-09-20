"""Public-panel selection semantics: independent frequencies, masks and alleles."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'))
from build_roh_population_panel import build,load_mask,included


def inputs(root,bad_ref=False):
    reference=root/'reference.fa'
    with reference.open('w') as stream:
        for chrom in range(1,23):stream.write(f'>{chrom}\n'+ 'A'*60000+'\n')
    pysam.faidx(str(reference))
    mask=root/'mask.bed'
    with mask.open('w') as stream:
        for chrom in range(1,23):
            stream.write(f'chr{chrom}\t0\t20000\tstrict\nchr{chrom}\t21000\t60000\tstrict\n')
    header=pysam.VariantHeader()
    for chrom in range(1,23):header.contigs.add(str(chrom),length=60000)
    for key in ['AF','AFR_AF','AMR_AF','EAS_AF','EUR_AF','SAS_AF']:
        header.info.add(key,'A','Float','Synthetic external population frequency')
    header.info.add('AN',1,'Integer','Synthetic called allele count')
    header.info.add('NS',1,'Integer','Synthetic sample count')
    source=root/'population.vcf'
    with pysam.VariantFile(str(source),'w',header=header) as writer:
        for chrom in range(1,23):
            for pos,alt in [(100,'G'),(200,'G'),(300,'G'),(10100,'G'),(10100,'C'),(20100,'G'),(30100,'G'),(30200,'G'),(50100,'G')]:
                row=writer.new_record(contig=str(chrom),start=pos-1,alleles=('C' if bad_ref and chrom==1 and pos==100 else 'A',alt))
                row.filter.add('PASS')
                for key in ['AF','AFR_AF','AMR_AF','EAS_AF','EUR_AF','SAS_AF']:row.info[key]=(0.5,)
                if pos==200:row.info['AF']=(0.01,)
                if pos==300:row.info['AFR_AF']=(0.01,)
                row.info['AN']=5008;row.info['NS']=2504
                writer.write(row)
    config=json.loads((Path(__file__).parents[1]/'configs/roh-population-panel.json').read_text())
    config.update(source_bytes=source.stat().st_size,source_publisher_md5=hashlib.md5(source.read_bytes()).hexdigest(),
                  source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
                  mask_sha256=hashlib.sha256(mask.read_bytes()).hexdigest())
    return source,reference,mask,config


class PopulationPanel(unittest.TestCase):
    def test_independent_selection_masks_duplicate_loci_and_preserves_frequencies(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source,reference,mask,config=inputs(root)
            result=build(source,reference,mask,config,root/'output')
            self.assertEqual(result['selected_markers'],66)
            self.assertEqual(result['counts']['source_records'],198)
            self.assertEqual(result['counts']['duplicate_position_records'],44)
            self.assertEqual(result['counts']['outside_passed_mask'],22)
            self.assertEqual(result['counts']['excluded_spacing'],22)
            with pysam.VariantFile(str(root/'output/panel.vcf.gz')) as reader:
                self.assertEqual(list(reader.header.samples),[])
                rows=list(reader.fetch('1'))
                self.assertEqual([r.pos for r in rows],[100,30100,50100])
                self.assertTrue(all(r.info['POP_AF']==0.5 and r.info['POP_AFR_AF']==0.5 for r in rows))
            with pysam.TabixFile(str(root/'output/alleles.tsv.gz')) as reader:
                self.assertEqual(list(reader.fetch('1',99,100)),['1\t100\tA,G'])

    def test_mask_inclusion_uses_zero_based_half_open_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'mask.bed';p.write_text('chr1\t10\t20\tstrict\n')
            intervals,_=load_mask(p,{'1':30})
            self.assertFalse(included(intervals,'1',9))
            self.assertTrue(included(intervals,'1',10))
            self.assertTrue(included(intervals,'1',19))
            self.assertFalse(included(intervals,'1',20))
            p.write_text('chr1\t10\t20\tstrict\nchr1\t19\t25\tstrict\n')
            with self.assertRaisesRegex(ValueError,'overlapping'):load_mask(p,{'1':30})

    def test_selected_ref_mismatch_is_not_silently_excluded(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source,reference,mask,config=inputs(root,bad_ref=True)
            with self.assertRaisesRegex(ValueError,'selected REF'):
                build(source,reference,mask,config,root/'output')

    def test_changed_population_and_mask_sources_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source,reference,mask,config=inputs(root)
            for field in ['source_sha256','source_publisher_md5','mask_sha256']:
                with self.assertRaises(ValueError):
                    build(source,reference,mask,dict(config,**{field:'wrong'}),root/field)


if __name__=='__main__':unittest.main()
