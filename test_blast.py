import tempfile, unittest
from pathlib import Path
from plasmid_evidence import blast_import, blast_summary
class BlastTests(unittest.TestCase):
    def test_pairwise_wrapped_title_and_exact_identity(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'report.txt'
            p.write_text('Query= NODE_93\n\nLength=2502\n>CP120590.2 Priestia flexa strain PRO116 chromosome,\n complete genome\nLength=4000000\n Identities = 2500/2502 (99%), Gaps = 1/2502 (0%)\n')
            h=blast_import([p],{'NODE_93':'A'*2502})['NODE_93']
            self.assertAlmostEqual(h[0]['identity'],99.92006395)
            self.assertIn('chromosome',blast_summary(h)); self.assertIsNone(h[0]['qcov'])
    def test_tsv_and_conflicting_locations(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'hits.tsv'
            p.write_text('q\ta\t99.9\t100\t100\t0\t200\tSpecies chromosome\nq\tb\t99.8\t100\t100\t0\t199\tSpecies plasmid p1\n')
            h=blast_import([p],{'q':'A'*100})['q']
            self.assertEqual(h[0]['qcov'],100)
            self.assertIn('chromosome, plasmid',blast_summary(h))
    def test_unknown_query_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'report.txt';p.write_text('Query= wrong\n')
            with self.assertRaises(ValueError): blast_import([p],{'q':'A'})
    def test_empty_is_not_negative(self):
        self.assertIn('completion not verified',blast_summary([]))
if __name__=='__main__': unittest.main()
