import unittest
from io import StringIO
import pandas as pd

from qapa import extend

class ExtendTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Last exons have matching 5' prime
        input_shared = """
seqnames	start	end	name	utr_length	strand	lastexon_cds_start	lastexon_cds_end	name2	exonStarts	exonEnds
chr1	1781220	1781554	ENSRNOT00000073528_ENSRNOG00000049505.2	110	+	1781220	1781444	ENSRNOG00000049505.2	1771720,1772195,1774162,1776350,1777791,1779038,1781220	1771889,1772318,1774343,1776486,1777943,1779161,1781554
chr1	1781220	1782091	ENSRNOT00000080138_ENSRNOG00000049505.2	647	+	1781220	1781444	ENSRNOG00000049505.2	1771709,1772195,1773278,1774162,1776350,1777791,1779038,1781220	1771889,1772318,1773356,1774343,1776486,1777943,1779161,1782091
        """
        cls.df_shared = pd.read_table(StringIO(input_shared))

        # Last exons don't have matching 5' prime
        input_non_shared = """
seqnames	start	end	name	utr_length	strand	lastexon_cds_start	lastexon_cds_end	name2	exonStarts	exonEnds
chr1	1777791	1781554	ENSRNOT00000073528_ENSRNOG00000049505.2	110	+	1781220	1781444	ENSRNOG00000049505.2	1771720,1772195,1774162,1776350,1777791,1779038,1781220	1771889,1772318,1774343,1776486,1777943,1779161,1781554
chr1	1781220	1782091	ENSRNOT00000080138_ENSRNOG00000049505.2	647	+	1781220	1781444	ENSRNOG00000049505.2	1771709,1772195,1773278,1774162,1776350,1777791,1779038,1781220	1771889,1772318,1773356,1774343,1776486,1777943,1779161,1782091
        """
        cls.df_non_shared = pd.read_table(StringIO(input_non_shared))

    def test_extend_0(self):
        result = extend.extend_5prime(self.df_shared, 0)
        self.assertEqual(result.start.tolist(), self.df_shared.start.tolist())
        self.assertEqual(result.end.tolist(), self.df_shared.end.tolist())
        
        result = extend.extend_5prime(self.df_non_shared, 0)
        self.assertEqual(result.start.tolist(), [1781220]*2)
        self.assertEqual(result.end.tolist(), self.df_non_shared.end.tolist())
        
    def test_extend_1(self):
        result = extend.extend_5prime(self.df_shared, 1)
        self.assertEqual(result.start.tolist(), self.df_shared.start.tolist())
        self.assertEqual(result.end.tolist(), self.df_shared.end.tolist())

        result = extend.extend_5prime(self.df_non_shared, 1)
        self.assertEqual(result.start.tolist(), [1781220]*2)
        self.assertEqual(result.end.tolist(), self.df_non_shared.end.tolist())
        
    def test_extend_2(self):
        result = extend.extend_5prime(self.df_shared, 1)
        self.assertEqual(result.start.tolist(), self.df_shared.start.tolist())
        self.assertEqual(result.end.tolist(), self.df_shared.end.tolist())

        result = extend.extend_5prime(self.df_non_shared, 2)
        self.assertEqual(result.start.tolist(), [1777791]*2)
        self.assertEqual(result.end.tolist(), self.df_non_shared.end.tolist())

if __name__ == '__main__':
    unittest.main()
