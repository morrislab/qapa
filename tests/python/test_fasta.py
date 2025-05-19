from pathlib import Path
import unittest
from io import StringIO

from qapa import fasta


DIR = Path(__file__).parent


class FastaTestCase(unittest.TestCase):

    def setUp(self):
        self.fasta_file = DIR / "files/test.fa"
        self.bed_file = DIR / "files/test_fasta_regions.bed"
        
    def test_get_sequences(self):
        # Get sequences using the test bed file and fasta file
        result = fasta.get_sequences(self.bed_file, self.fasta_file)
        
        with open(result.seqfn, 'r') as f:
            sequence_content = f.read()
            
        # Check that name field of BED file is extracted as the sequence identifier
        self.assertIn(">seq2:0:20:+", sequence_content)
        self.assertIn(">seq3:0:20:-", sequence_content)
        
        # Check that the sequences themselves are extracted correctly
        self.assertIn("AAAACTCTCACACAAGAAAA", sequence_content)
        self.assertIn("GATGGGTCTTTATTTAAGAA", sequence_content)

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

