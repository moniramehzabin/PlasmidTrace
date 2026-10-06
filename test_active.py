import unittest, tempfile, subprocess, sys
from pathlib import Path
from reconstruct import cigar_crosses, junction_counts, candidate_matches, spades_launcher
class ActiveTests(unittest.TestCase):
    def test_junction_boundaries(self):
        self.assertTrue(cigar_crosses('100M',450,500,20))
        self.assertFalse(cigar_crosses('100M',481,500,20))
        self.assertFalse(cigar_crosses('40M10D60M',450,500,20))
        self.assertFalse(cigar_crosses('40M10N60M',450,500,20))
        self.assertTrue(cigar_crosses('10S100M',450,500,20))
    def test_read_filters_and_pair_names(self):
        def row(name,flag,cigar='100M',nm=0):return f'{name}\t{flag}\tj\t451\t30\t{cigar}\t*\t0\t0\tAAAA\tIIII\tNM:i:{nm}\n'
        lines=[row('pair',64),row('pair',128),row('other',0,'10S100M'),row('bad',2048),row('duplicate',1024),row('mismatch',0,nm=1)]
        r=junction_counts(lines,'j',500,20,30)
        self.assertEqual(r,dict(crossing_reads=4,distinct_pair_names=3,perfect_unclipped_crossing_reads=2))
    def test_candidate_match_union_and_identity(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'matches';p.write_text('q\ts\t99\t1\t60\t100\t60\t90\nq\ts\t99\t40\t100\t100\t61\t90\nq\tt\t90\t1\t100\t100\t100\t90\n')
            r=candidate_matches(p,{'q':'A'*100},{'s':'A'*100,'t':'A'*100},95,80)
            self.assertEqual(len(r),1);self.assertEqual(r[0]['query_coverage_percent'],100)
    def test_memory_wrapper_overrides_duplicate_flags(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);exe=root/'fake_spades';exe.write_text('#!/usr/bin/env python3\nimport sys\nprint(" ".join(sys.argv[1:]))\n');exe.chmod(0o755)
            wrapper=root/'launcher';spades_launcher(wrapper,exe,7)
            r=subprocess.run([str(wrapper),'-o','out','--memory','5','-m','1024'],capture_output=True,text=True,check=True)
            self.assertEqual(r.stdout.strip(),'-o out -m 7')
            r=subprocess.run([str(wrapper),'--version'],capture_output=True,text=True,check=True)
            self.assertEqual(r.stdout.strip(),'--version')
class VersionOutputTests(unittest.TestCase):
    def test_non_utf8_version_bytes_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            info=subprocess.run([sys.executable,'-c',"import sys;sys.stdout.buffer.write(bytes([0xab]))"],capture_output=True)
            p=Path(d)/'version.txt';p.write_bytes(info.stdout+info.stderr)
            self.assertEqual(p.read_bytes(),bytes([0xab]))
if __name__=='__main__':unittest.main()
