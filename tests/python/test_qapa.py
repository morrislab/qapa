import unittest
from unittest.mock import Mock, patch, ANY
import argparse
from tempfile import NamedTemporaryFile
from io import StringIO
import pandas as pd

from qapa import qapa

class QapaTestCase(unittest.TestCase):

    @patch('sys.stderr', new_callable=StringIO)
    def test_getoptions_bad_subcommand_error(self, mock_stderr):
        # sub-module must be specified
        with self.assertRaises(SystemExit):
            qapa.getoptions(['bad'])
        self.assertIn('argument subcommand: invalid choice:',
                mock_stderr.getvalue())
        self.assertIn("choose from 'build', 'fasta', 'quant'",
                mock_stderr.getvalue())

    @patch('sys.stderr', new_callable=StringIO)
    def test_getoptions_no_subcommand_error(self, mock_stderr):
        with self.assertRaises(SystemExit):
            qapa.getoptions([])
        self.assertIn("Specify sub-command: 'build', 'fasta', 'quant'",
                mock_stderr.getvalue())

    @patch('qapa.collapse.merge_bed')
    @patch('qapa.extend.main')
    @patch('qapa.annotate.main')
    @patch('qapa.extract.main')
    def test_build(self, mock_extract, mock_anno, mock_extend, mock_collapse):
        mock_collapse.return_value = pd.DataFrame()

        db = NamedTemporaryFile()
        gencode = NamedTemporaryFile()
        polyasite = NamedTemporaryFile()
        genepred = NamedTemporaryFile()
        args = qapa.getoptions(['build', '--db', db.name, '-g', gencode.name,
                                '-p', polyasite.name, genepred.name])
        self.assertIsInstance(args, argparse.Namespace)
        qapa.build(args)
        mock_extract.assert_called_once_with(args, ANY)
        mock_anno.assert_called_once_with(args, ANY, ANY)
        mock_extend.assert_called_once_with(args, ANY)
        mock_collapse.assert_called_once_with(args, ANY)

    @patch('qapa.fasta.main')
    def test_fasta(self, mock_fasta):
        genome = NamedTemporaryFile()
        bed = NamedTemporaryFile()
        output = NamedTemporaryFile()
        args = qapa.getoptions(['fasta', '-f', genome.name, bed.name, output.name])
        self.assertIsInstance(args, argparse.Namespace)
        qapa.fetch_sequences(args)
        mock_fasta.assert_called_once_with(args)


    @patch('os.system')
    def test_quant(self, mock_sys):
        # should call create_merged_data.R and compute_pau.R 
        db = NamedTemporaryFile()
        args = qapa.getoptions(['quant', '--db', db.name, 'test_1.sf',
                                'test_2.sf'])
        qapa.quant(args)
        self.assertRegexpMatches(mock_sys.call_args_list[0][0][0],
            r'create_merged_data.R')
        self.assertRegexpMatches(mock_sys.call_args_list[1][0][0],
            r'compute_pau.R')


if __name__ == '__main__':
    unittest.main()

