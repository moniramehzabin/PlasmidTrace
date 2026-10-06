import unittest, tempfile, subprocess, sys, csv
from pathlib import Path
class GraphContextTests(unittest.TestCase):
    def test_isolated_retained_connected_warned(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);u=root/'u';u.mkdir()
            (u/'assembly.fasta').write_text('>12\nACGG\n>28\nAACC\n>3\nACGT\n>14\nAAAA\n')
            (u/'assembly.gfa').write_text('S\t12\tACGG\nS\t28\tAACC\nS\t3\tACGT\nS\t14\tAAAA\nL\t12\t+\t12\t+\t0M\nL\t28\t+\t28\t+\t0M\nL\t28\t+\t3\t-\t0M\nL\t14\t+\t28\t+\t0M\n')
            r=subprocess.run([sys.executable,str(Path(__file__).with_name('plasmid_evidence.py')),'--unicycler-dir',str(u),'--out',str(root/'out')],capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stderr)
            with (root/'out/evidence.tsv').open() as h: rows={r['contig']:r for r in csv.DictReader(h,delimiter='\t')}
            self.assertEqual(rows['12']['status'],'candidate for review')
            self.assertEqual(rows['28']['status'],'connected graph loop warning')
            self.assertEqual(set(rows['28']['external_graph_neighbors'].split(',')),{'3','14'})
            self.assertEqual((root/'out/candidates.fasta').read_text(),'>12\nACGG\n')
if __name__=='__main__':unittest.main()
