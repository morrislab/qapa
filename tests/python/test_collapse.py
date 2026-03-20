import unittest
from collections import namedtuple

from qapa import collapse
from qapa.collapse import Interval


class CollapseTestCase(unittest.TestCase):
    def setUp(self):
        UTR = namedtuple(
            "UTR",
            [
                "seqnames",
                "start",
                "end",
                "name",
                "utr_length",
                "strand",
                "lastexon_cds_start",
                "lastexon_cds_end",
                "name2",
                "blockCount",
                "blockSizes",
                "blockStarts",
            ],
        )

        utr_a = UTR(
            "chrX",
            "74035282",
            "74036330",
            "ENSMUST00000000001_Mecp2",
            "961",
            "-",
            "74036243",
            "74036330",
            "Mecp2",
            "1",
            "1042",
            "0",
        )
        utr_b = UTR(
            "chrX",
            "74035283",
            "74036330",
            "ENSMUST00000000002_Mecp2",
            "960",
            "-",
            "74036243",
            "74036330",
            "Mecp2",
            "1",
            "1042",
            "0",
        )

        self.interval_a = Interval(utr_a)
        self.interval_b = Interval(utr_b)

    def test_merge(self):
        self.interval_a.merge(self.interval_b)
        self.assertEqual(self.interval_a.start, 74035282)
        self.assertEqual(self.interval_a.end, 74036330)
        self.assertEqual(self.interval_a.start2, 74036243)
        self.assertEqual(self.interval_a.end2, 74036330)
        self.assertEqual(
            self.interval_a.name, "ENSMUST00000000001_Mecp2,ENSMUST00000000002_Mecp2"
        )

    def test_overlaps(self):
        result = collapse.overlaps(self.interval_a, self.interval_b, 10)
        self.assertTrue(result)

        self.interval_a.start = 73036272
        result = collapse.overlaps(self.interval_a, self.interval_b, 10)
        self.assertFalse(result)


if __name__ == "__main__":
    unittest.main()
