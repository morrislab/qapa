from pathlib import Path
import unittest
from io import StringIO
import tempfile
import os

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
     
    def test_add_decoys(self):

        extracted_sequences = fasta.get_sequences(self.bed_file, self.fasta_file)
        
        # Create temporary files for decoys.txt and genome-appended fasta
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as decoys_file:
            decoys_path = decoys_file.name
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as output_file:
            output_path = output_file.name
            
        try:
            # Copy extracted sequences to tempfile for use with add_decoys
            with open(extracted_sequences.seqfn, 'r') as src, open(output_path, 'w') as dest:
                dest.write(src.read())

            with open(output_path, 'a') as fout:
                fasta.add_decoys(self.fasta_file, decoys_path, fout=fout)

            # Validate output fasta file (all genome sequence IDs are appended)
            with open(output_path, 'r') as f:
                output_content = f.read()
            
            self.assertIn(">seq2:0:20:+", output_content)
            self.assertIn(">seq3:0:20:-", output_content)
            # genome fasta sequence IDs
            self.assertIn(">seq1", output_content)
            self.assertIn(">seq2", output_content)
            self.assertIn(">seq3", output_content)
            self.assertIn(">seq4", output_content)
            
            # Validate decoys.txt file (genome sequence IDs printed one per line)
            expected_ids = {"seq1", "seq2", "seq3", "seq4"}
            with open(decoys_path, 'r') as f:
                for line in f:
                    stripped_line = line.strip()
                    self.assertEqual(stripped_line.count(" "), 0,
                                     "Each line should contain only one sequence ID")
                    self.assertIn(stripped_line, expected_ids, f"ID {stripped_line} not in expected IDs")
                    expected_ids.remove(stripped_line)
            
            self.assertEqual(len(expected_ids), 0, f"Missing IDs in decoys file: {expected_ids}")

        finally:
            os.unlink(output_path)
            os.unlink(decoys_path)


if __name__ == '__main__':
    unittest.main()

