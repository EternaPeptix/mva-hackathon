import importlib.util
from pathlib import Path
import unittest


P = Path(__file__).resolve().parents[1] / 'scripts/audit_bub1b_isoforms.py'
SPEC = importlib.util.spec_from_file_location('bub1b_isoforms', P)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


class IsoformAuditTests(unittest.TestCase):
    def rows(self):
        # Children deliberately precede parents; the readthrough gene is separate.
        return [
            '15\t.\tCDS\t210\t219\t.\t+\t0\tParent=transcript:a',
            '15\t.\texon\t200\t229\t.\t+\t.\tParent=transcript:a',
            '15\t.\tCDS\t110\t119\t.\t+\t0\tParent=transcript:a,transcript:b',
            '15\t.\texon\t100\t129\t.\t+\t.\tParent=transcript:a,transcript:b',
            '15\t.\tmRNA\t100\t229\t.\t+\t.\tID=transcript:a;Parent=gene:g;transcript_id=a;biotype=protein_coding',
            '15\t.\tmRNA\t100\t129\t.\t+\t.\tID=transcript:b;Parent=gene:g;transcript_id=b;biotype=protein_coding',
            '15\t.\tgene\t100\t229\t.\t+\t.\tID=gene:g',
            '15\t.\tmRNA\t100\t229\t.\t+\t.\tID=transcript:readthrough;Parent=gene:other;transcript_id=readthrough;biotype=protein_coding']

    def test_parent_join_and_exon_cds_distinction(self):
        rows = AUDIT.features(self.rows(), start=1, end=500)
        result = AUDIT.summarize(rows, gene='gene:g', targets=(110, 219, 229, 230))
        a, b = result['transcripts']
        self.assertEqual([t['id'] for t in result['transcripts']], ['a', 'b'])
        self.assertEqual(a['loci']['110']['offset_in_annotated_CDS'], 1)
        self.assertEqual(a['loci']['219']['offset_in_annotated_CDS'], 20)
        self.assertTrue(a['loci']['229']['in_exon'])
        self.assertFalse(a['loci']['229']['in_CDS'])
        self.assertFalse(a['loci']['230']['in_exon'])
        self.assertFalse(b['loci']['219']['in_exon'])

    def test_invalid_coding_geometry_is_not_silently_accepted(self):
        rows = self.rows()
        rows[0] = rows[0].replace('210\t219', '210\t239')
        with self.assertRaisesRegex(ValueError, 'CDS outside'):
            AUDIT.summarize(AUDIT.features(rows, start=1, end=500), gene='gene:g')
        rows = self.rows()
        rows[0] = rows[0].replace('\t+\t', '\t-\t')
        with self.assertRaisesRegex(ValueError, 'inconsistent feature strand'):
            AUDIT.summarize(AUDIT.features(rows, start=1, end=500), gene='gene:g')


if __name__ == '__main__':
    unittest.main()
