import csv, subprocess, sys, tempfile, unittest
from pathlib import Path
from plasmid_evidence import rrna_intervals, classify

class Tests(unittest.TestCase):
    def test_rrna_overrides_prediction(self):
        self.assertEqual(classify('plasmid',.94,True,True)[0],'rRNA-dominated warning')
    def test_missing_not_negative(self):
        self.assertIn('supplied checks',classify(None,None,False,False)[0])
    def test_overlapping_gff_union(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'r.gff'
            p.write_text('x\tb\trRNA\t1\t60\t.\t+\t.\tName=16S_rRNA\nx\tb\trRNA\t40\t80\t.\t+\t.\tName=23S_rRNA\n')
            self.assertEqual(rrna_intervals(p,{'x':'A'*100})['x'],.8)
    def test_end_to_end_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); assembly=root/'a.fa'; gff=root/'r.gff'; mob=root/'mob.tsv'; out=root/'out'
            n='NODE_129_length_1648_cov_945.307447'
            assembly.write_text('>'+n+'\n'+'A'*1648+'\n>NODE_1_length_12000_cov_70\n'+'ACGT'*3000+'\n>other\nACGT\n')
            gff.write_text(f'{n}\tbarrnap\trRNA\t74\t1623\t0\t-\t.\tName=16S_rRNA\n')
            mob.write_text('contig_id\tmolecule_type\n'+n+'\tplasmid\n')
            cmd=[sys.executable,str(Path(__file__).with_name('plasmid_evidence.py')),'--assembly',str(assembly),'--out',str(out),'--rrna-gff',str(gff),'--mob-report',str(mob)]
            subprocess.run(cmd,check=True,capture_output=True)
            with (out/'evidence.tsv').open() as h: rows=list(csv.DictReader(h,delimiter='\t'))
            self.assertEqual(rows[0]['status'],'rRNA-dominated warning')
            self.assertEqual(rows[2]['coverage_header'],'not assessed')
            self.assertEqual((out/'candidates.fasta').read_text(),'')
            self.assertNotEqual(subprocess.run(cmd,capture_output=True).returncode,0)
if __name__=='__main__': unittest.main()
