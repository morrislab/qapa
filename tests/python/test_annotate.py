import unittest
import sys
import re
import pybedtools
from qapa import annotate as anno

class AnnotateTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        example = "chrX	74026591	74036494	ENSMUST00000100750_Mecp2	8825	-	74035416	74036494	Mecp2	74026591,74036981,74079908,74085509	74036494,74037332,74080032,74085669"
        cls.bed = pybedtools.BedTool(example, from_string=True)

    def test_extend_feature(self):
        l = 1
        feature = self.bed[0]
        target = anno.extend_feature(feature, length=l)
        self.assertEqual(target.start, 74026591 - l)
        self.assertEqual(target.end, 74036494)

    
    def test_gene_at_beginning_of_chr_1(self):
        example = "chr17_KI270861v1_alt	0	5793	ENST00000634102_SLC43A2 5631	-	5631	5793	SLC43A2 0,6624,8277,13231,15940,20829,21296,21593,23214,43222,45006,46589,57742,58821 5793,6748,8351,13364,16079,20976,21499,21727,23307,43299,45062,46797,57948,58864"
        bed = pybedtools.BedTool(example, from_string=True)
        l = 24
        feature = bed[0]
        target = anno.extend_feature(feature, length=l)
        self.assertEqual(target.start, 0)
        self.assertEqual(target.end, 5793)
        
        target = anno.restore_feature(target, length=l)
        self.assertEqual(target.start, 0)

    def test_gene_at_beginning_of_chr_2(self):
        example = "chr17_KI270861v1_alt	24	5793	ENST00000634102_SLC43A2 5631	-	5631	5793	SLC43A2 24,6624,8277,13231,15940,20829,21296,21593,23214,43222,45006,46589,57742,58821 5793,6748,8351,13364,16079,20976,21499,21727,23307,43299,45062,46797,57948,58864"
        bed = pybedtools.BedTool(example, from_string=True)
        l = 24
        feature = bed[0]
        target = anno.extend_feature(feature, length=l)
        self.assertEqual(target.start, 0)
        self.assertEqual(target.end, 5793)

        target = anno.restore_feature(target, length=l)
        self.assertEqual(target.start, 24)

    def test_preprocess_gencode(self):
        result = anno.preprocess_gencode_polya("python/files/gencode.polya.example.bed")

        # Test no polya_signals
        count = 0
        for item in result:
            self.assertEqual(item.name, 'polyA_site')
            count += 1

        # Should be four entries total
        self.assertEqual(count, 4)

    def test_preprocess_polyasite_v1(self):
        with self.assertLogs("qapa", level="INFO") as cm:
            result = anno.preprocess_polyasite("python/files/polyasite_v1.example.bed", 2)
        self.assertIn("Detected PolyASite version 1", cm.output[1])

        count = 0
        for item in result:
            self.assertTrue(re.search(r"(DS|TE)$", item.name))
            self.assertGreaterEqual(int(item[4]), 2)
            count += 1
        self.assertEqual(count, 2)

    def test_preprocess_polyasite_v2(self):
        with self.assertLogs("qapa", level="INFO") as cm:
            result = anno.preprocess_polyasite("python/files/polyasite_v2.example.bed", 2)
        self.assertIn("Detected PolyASite version 2", cm.output[1])

        count = 0
        for item in result:
            self.assertTrue(re.match(r"^chr", item[0]))
            self.assertEqual(item[0], "chrX")
            self.assertGreaterEqual(int(item[4]), 2)
            count += 1
        self.assertEqual(count, 3)


if __name__ == '__main__':
    unittest.main()
