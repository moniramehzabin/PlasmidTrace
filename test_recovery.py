import tempfile, unittest
from pathlib import Path
from recover_graphs import inventory, insert_sizes

class RecoveryTests(unittest.TestCase):
    def test_missing_graphs_refused(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError): inventory(Path(d))
    def test_only_completed_nonempty_graphs(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            (p/'001_spades_graph_k027.gfa').write_text('S\t1\tACGT\n')
            (p/'001_spades_graph_k111.gfa').write_text('S\t1\tACGT\n')
            (p/'001_spades_graph_k119.gfa').touch()
            self.assertEqual([k for k,_ in inventory(p)],[27,111])
    def test_insert_metadata_median(self):
        self.assertEqual(insert_sizes('Insert size = 300, deviation = 50\nInsert size = 320, deviation = 70'),(310,60))
    def test_no_guessed_insert_metadata(self):
        with self.assertRaises(ValueError): insert_sizes('Interrupted')
