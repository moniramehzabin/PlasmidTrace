import tempfile, unittest, subprocess, sys, json
from pathlib import Path
from plasmid_evidence import exact_aliases, blast_import
class ReuseTests(unittest.TestCase):
    def test_unique_exact_alias_and_custom_blast(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'alias.fa';p.write_text('>renamed\nACGG\n')
            aliases=exact_aliases(p,{'2':'ACGG'})
            self.assertEqual(aliases,{'renamed':'2'})
            b=Path(d)/'blast.tsv';b.write_text('renamed\tCP1\t100\t4\t4\t100\t1\t4\t20\t23\t0\t8\n')
            with self.assertRaises(ValueError):blast_import([b],{'2':'ACGG'},aliases)
            self.assertEqual(blast_import([b],{'2':'ACGG'},aliases,True)['2'][0]['query_intervals'],[[1,4]])
    def test_ambiguous_alias_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'alias.fa';p.write_text('>renamed\nACGG\n')
            with self.assertRaises(ValueError):exact_aliases(p,{'1':'ACGG','2':'ACGG'})
    def test_reuse_end_to_end(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);u=root/'unicycler';u.mkdir()
            (u/'assembly.fasta').write_text('>1\nAACCGG\n>2\nACGG\n')
            (u/'assembly.gfa').write_text('S\t1\tAACCGG\nS\t2\tACGG\nL\t2\t+\t2\t+\t0M\n')
            script=Path(__file__).with_name('plasmid_evidence.py')
            r=subprocess.run([sys.executable,str(script),'--unicycler-dir',str(u),'--out',str(root/'out')],capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stderr)
            self.assertIn('Same-orientation', (root/'out/report.md').read_text())
            self.assertEqual((root/'out/candidates.fasta').read_text(),'>2\nACGG\n')
if __name__=='__main__':unittest.main()
