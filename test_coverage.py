import unittest, tempfile
from pathlib import Path
from plasmid_evidence import blast_import, reference_summaries
class CoverageTests(unittest.TestCase):
    def test_reverse_coordinates_overlap_and_reference_grouping(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'hits.txt'
            p.write_text('Query= q\nLength=100\n>a Species chromosome\nLength=500\n Identities = 60/60 (100%)\nQuery  1  '+ 'A'*60 +'  60\n Identities = 61/61 (100%)\nQuery  100  '+ 'A'*61 +'  40\n>b Other plasmid\nLength=300\n Identities = 10/10 (100%)\nQuery  1  AAAAAAAAAA  10\n')
            hits=blast_import([p],{'q':'A'*100})['q']
            rows=reference_summaries(hits,100)
            self.assertEqual(rows[0]['query_coverage_percent'],100)
            self.assertEqual(rows[0]['alignment_sections'],2)
            self.assertEqual(rows[0]['longest_alignment_bp'],61)
            self.assertEqual(rows[1]['query_coverage_percent'],10)
    def test_missing_coordinates_not_full_coverage(self):
        h=dict(accession='a',title='chromosome',identity=100,qcov=None,alignment_length=100,query_intervals=[])
        self.assertIsNone(reference_summaries([h],100)[0]['query_coverage_percent'])
if __name__=='__main__': unittest.main()
