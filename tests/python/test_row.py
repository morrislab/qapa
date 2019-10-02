import unittest
import sys
from qapa import extract as ex
from qapa.extract import get_stripped_name

class RowTestCase(unittest.TestCase):

    def setUp(self):
        example_row = '143	ENSMUST00000100750.9	chrX_random	-	74026591	74085669	74035416	74079934	4	74026591,74036981,74079908,74085509,	74036494,74037332,74080032,74085669,	0	Mecp2	cmpl	cmpl	2,2,0,-1,'
        self.row = ex.Row(example_row)

    def test_3utr_length(self):
        target = self.row.get_3utr_length()
        expected = self.row.cdsStart - self.row.txStart
        self.assertEqual(target, expected)


    def test_random_chromosome(self):
        target = self.row.is_on_random_chromosome()
        self.assertTrue(target)

    def test_random_chromosome_no_chr(self):
        my_row = '143	ENSMUST00000100750.9	X	-	74026591	74085669	74035416	74079934	4	74026591,74036981,74079908,74085509,	74036494,74037332,74080032,74085669,	0	Mecp2	cmpl	cmpl	2,2,0,-1,'
        my_row = ex.Row(my_row)
        target = my_row.is_on_random_chromosome()
        self.assertFalse(target)

    def test_chrY(self):
        my_row = '143	ENSMUST00000100750.9	chrY	-	74026591	74085669	74035416	74079934	4	74026591,74036981,74079908,74085509,	74036494,74037332,74080032,74085669,	0	Mecp2	cmpl	cmpl	2,2,0,-1,'
        my_row = ex.Row(my_row)
        target = my_row.is_on_random_chromosome()
        self.assertFalse(target)


    def test_extract_last_exon(self):
        target = self.row.extract_last_exon(n=1, min_utr_length=0)
        self.assertEqual(target[1], 74026591, "Start coord not equal")
        self.assertEqual(target[2], 74036494, "End coord not equal")


    def test_invalid_gene_pred_format(self):
        my_row = '143   ENSMUST00000100750.9    chrY    -   74026591    74085669    74035416    74079934    4   74026591,74036981,74079908,74085509,    74036494,74037332,74080032,74085669,'
        with self.assertRaises(ValueError):
            ex.Row(my_row)

    def test_get_stripped_name(self):
        target = get_stripped_name(self.row.name)
        self.assertEqual(target, "ENSMUST00000100750")

        target = get_stripped_name("ENST00000100750.11")
        self.assertEqual(target, "ENST00000100750")

        target = get_stripped_name("ENSRNOT00000100750.11")
        self.assertEqual(target, "ENSRNOT00000100750")

        target = get_stripped_name("ENST00123")
        self.assertEqual(target, "ENST00123")

if __name__ == '__main__':
    unittest.main()
