import csv, tempfile, subprocess, sys, unittest
from pathlib import Path
from finder_adapter import import_hits

HEADER='Plasmid\tIdentity\tQuery / Template length\tContig\tPosition in contig\tAccession number\n'
class FinderTests(unittest.TestCase):
    def test_marker_thresholds_and_contig_validation(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'results.tsv'
            p.write_text(HEADER+'repA\t99\t70/100\tx\t1..70\tA1\n'+'repB\t90\t100/100\tx\t1..100\tA2\n')
            self.assertEqual(len(import_hits(p,{'x':'A'*100},.95,.6)['x']),1)
            with self.assertRaises(ValueError):import_hits(p,{'y':'A'*100},.95,.6)
    def test_empty_header_table_and_bad_schema(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'results.tsv';p.write_text(HEADER)
            self.assertEqual(import_hits(p,{'x':'AAA'},.95,.6),{'x':[]})
            p.write_text('');
            with self.assertRaises(ValueError):import_hits(p,{'x':'AAA'},.95,.6)
    def test_import_drives_candidates_and_rrna_priority(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);a=d/'a.fa';a.write_text('>x\n'+'A'*100+'\n>y\n'+'A'*100+'\n')
            t=d/'t.tsv';t.write_text(HEADER+'repA\t99\t100/100\tx\t1..100\tA\n'+'repA\t99\t100/100\ty\t1..100\tA\n')
            g=d/'r.gff';g.write_text('y\ttest\trRNA\t1\t100\t.\t+\t.\tName=16S\n')
            subprocess.run([sys.executable,str(Path(__file__).with_name('plasmid_evidence.py')),'--assembly',str(a),'--plasmidfinder-results',str(t),'--rrna-gff',str(g),'--out',str(d/'out')],check=True,capture_output=True)
            with (d/'out/evidence.tsv').open() as handle:
                rows=list(csv.DictReader(handle,delimiter='\t'))
            self.assertEqual(rows[0]['status'],'candidate for review')
            self.assertEqual(rows[1]['status'],'rRNA-dominated warning')
            self.assertEqual(rows[1]['plasmidfinder_markers'],'repA')
