import unittest,tempfile,subprocess,sys,json
from pathlib import Path
from annotation_evidence import read_gff,read_genbank,nominates
from finder_adapter import command
from types import SimpleNamespace
class AnnotationTests(unittest.TestCase):
 def test_rules_no_conjugate_false_positive(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'a.gff';p.write_text('x\ttest\tCDS\t1\t30\t.\t+\t0\tproduct=mycothiol%20conjugate%20amidase\n'+'x\ttest\tCDS\t31\t60\t.\t+\t0\tproduct=ParA%20family%20protein\n')
   h,_=read_gff(p,{'x':'A'*100});self.assertEqual(len(h['x']),1);self.assertFalse(nominates(h['x']))
 def test_replication_and_bad_coordinates(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'a.gff';p.write_text('x\ttest\tCDS\t1\t30\t.\t+\t0\tproduct=plasmid%20replication%20protein\n')
   h,_=read_gff(p,{'x':'A'*100});self.assertTrue(nominates(h['x']))
   with self.assertRaises(ValueError):read_gff(p,{'x':'AA'})
 def test_genbank_exact_dna_and_metadata_not_detection(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'a.gb';p.write_text('LOCUS       renamed 30 bp DNA linear\nFEATURES             Location/Qualifiers\n     source          1..30\n                     /plasmid="p1"\n     CDS             1..30\n                     /product="plasmid replication\n                     protein"\nORIGIN\n        1 '+ 'a'*30+'\n//\n')
   h,labels,aliases=read_genbank(p,{'x':'A'*30});self.assertTrue(nominates(h['x']));self.assertEqual(aliases,{'renamed':'x'})
   self.assertEqual(labels['x'],['p1'])
   with self.assertRaises(ValueError):read_genbank(p,{'x':'C'*30})
 def test_pf_module_requests_table(self):
  a=SimpleNamespace(plasmidfinder_python='/env/python',plasmidfinder_executable='unused',plasmidfinder_db=Path('/db'),plasmidfinder_identity=.95,plasmidfinder_coverage=.6)
  c=command(a,Path('/a'),Path('/out'));self.assertEqual(c[:4],['/env/python','-m','plasmidfinder','-x'])
 def test_annotation_only_candidate_end_to_end(self):
  with tempfile.TemporaryDirectory() as d:
   d=Path(d);a=d/'a.fa';a.write_text('>x\n'+'A'*100+'\n');g=d/'a.gff';g.write_text('x\tt\tCDS\t1\t90\t.\t+\t0\tproduct=plasmid%20replication%20protein\n')
   subprocess.run([sys.executable,str(Path(__file__).with_name('plasmid_evidence.py')),'--assembly',str(a),'--annotation-gff',str(g),'--out',str(d/'out')],capture_output=True,check=True)
   self.assertIn('>x',(d/'out/candidates.fasta').read_text())
