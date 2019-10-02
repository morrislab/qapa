import unittest
from io import StringIO

from qapa import fasta

class FastaTestCase(unittest.TestCase):

    def setUp(self):
        self.fasta_file = "python/files/test.fa"

    def test_filter_sequences(self):
        output = StringIO()
        with self.assertLogs("qapa", level="DEBUG") as cm:
            _ = fasta.filter_sequences(self.fasta_file, fout=output)
        self.assertIn("Skipping seq1", cm.output[0])
        self.assertIn("Skipping seq4", cm.output[1])
        results = output.getvalue()
        self.assertIn(">seq2", results)
        self.assertIn(">seq3", results)

if __name__ == '__main__':
    unittest.main()

