from __future__ import print_function
import re
import sys
import fileinput
import numpy as np
import pandas as pd
# import sqlite3


class Row:
    def __init__(self, row):
        l = row.rstrip().split("\t")
        self.name = l[1]
        self.chrom = l[2]
        self.strand = l[3]
        self.txStart = int(l[4])
        self.txEnd = int(l[5])
        self.cdsStart = int(l[6])
        self.cdsEnd = int(l[7])
        self.exonStarts = [int(x) for x in [_f for _f in l[9].split(",") if _f]]
        self.exonEnds = [int(x) for x in [_f for _f in l[10].split(",") if _f]]
        #self.exonStarts = [int(x) for x in filter(None, l[9].split(","))]
        #self.exonEnds = [int(x) for x in filter(None, l[10].split(","))]
        self.name2 = l[12]

        self.utr3 = [self.cdsEnd, self.txEnd]
        self.has_intron_in_3utr = self.cdsEnd < self.exonStarts[-1]

        if self.strand == "-":
            self.utr3 = [self.txStart, self.cdsStart]
            self.has_intron_in_3utr = self.cdsStart > self.exonEnds[0]

    def extract_last_exon(self, n=1, min_utr_length=0):
        bed = None
        name = self.get_stripped_name() + "_" + self.name2

        if not self.has_intron_in_3utr and \
                self.get_3utr_length() >= min_utr_length:
            if self.strand == "+":
                start = np.max([self.exonStarts[-n], self.cdsStart])
                bed = [self.chrom, start, self.txEnd, name,
                       self.get_3utr_length(), self.strand,
                       start, self.cdsEnd]
            else:
                end = np.min([self.exonEnds[n - 1], self.cdsEnd])
                bed = [self.chrom, self.txStart, end, name,
                       self.get_3utr_length(), self.strand, self.cdsStart,
                       end]
            bed.append(self.name2)
            bed.append(','.join([str(x) for x in self.exonStarts]))
            bed.append(','.join([str(x) for x in self.exonEnds]))
        else:
            pass
            # print >> sys.stderr, "Skipping " + self.name + " because it" + \
            #" contains an intron in 3' UTR"
        return bed

    def extract_3utr(self, min_utr_length=0):
        bed = None
        name = self.get_stripped_name() + "_" + self.name2

        if not self.has_intron_in_3utr and \
                self.get_3utr_length() >= min_utr_length:
            if self.strand == "+":
                bed = [self.chrom] + self.utr3 + [name,
                                                  self.get_3utr_length(),
                                                  self.strand]
            else:
                bed = [self.chrom] + self.utr3 + [name,
                                                  self.get_3utr_length(),
                                                  self.strand]
        return bed

    def get_3utr_length(self):
        return self.utr3[1] - self.utr3[0]

    def is_on_random_chromosome(self):
        return not re.match(r'chr[0-9XY]+$', self.chrom)

    def get_block_sizes(self, n):
        sizes = [0] * n
        for i in range(0, n):
            if self.strand == "+":
                sizes[i] = self.exonEnds[i - 1] - self.exonStarts[i - 1]
            else:
                sizes[i] = self.exonEnds[i] - self.exonStarts[i]
        return ",".join([str(x) for x in sizes])

    def get_block_starts(self, n):
        if self.strand == "+":
            return ",".join([str(x - self.exonStarts[-n])
                             for x in self.exonStarts[-n:]])
        return ",".join([str(x - self.txStart) for x in self.exonStarts[0:n]])

    def get_stripped_name(self):
        # If Gencode tables are supplied, the Ensembl transcript ID has a
        # version number appended to the ID. We want to strip this out.
        if re.match('^ENS(MUS)*T', self.name):
            m = re.match('^ENS(MUS)*T\d+', self.name)
            return m.group()
        return self.name


def main(args, fout=sys.stdout):

    # print "\t".join(["seqnames", "start", "end", "name", "utr_length", "strand",
                     #"lastexon_cds_start", "lastexon_cds_end", "name2",
                     #"exonStarts", "exonEnds"])

    # conn = sqlite3.connect(args.db)

    # query = "select gene_biotype, transcript_biotype from ensembl_id where transcript_id = ?"

    conn = pd.read_table(args.db)
    conn = conn.loc[:, ['Transcript stable ID', 'Gene type',
                        'Transcript type']].drop_duplicates()
    conn = conn.set_index(['Transcript stable ID'])

    c = 0
    n = 0
    for row in fileinput.input(args.annotation_file[0],
                               openhook=fileinput.hook_compressed):
        n = n + 1

        if fileinput.isfirstline():
            continue

        if re.match(r"^#", row):
            c = c + 1
            continue

        rowobj = Row(row)

        if rowobj.is_on_random_chromosome():
            c = c + 1
            continue

        # filter for only protein-coding genes
        # result = conn.execute(query, (rowobj.get_stripped_name(),))
        # result = result.fetchone()
        # if result is None or \
        #     not (result[0] == "protein_coding" and \
        #     result[1] == "protein_coding"):
        #         c = c + 1
        #         continue

        # filter for only protein-coding genes
        try:
            result = conn.loc[rowobj.get_stripped_name()]
            if isinstance(result, pd.DataFrame):
                result = result.iloc[0, ]
            if not (result['Gene type'] == "protein_coding" and
                    result['Transcript type'] == "protein_coding"):
                c = c + 1
                continue
        except KeyError:
            c = c + 1
            continue

        bed = rowobj.extract_last_exon()

        if bed is not None:
            fout.write("\t".join([str(x) for x in bed]) + "\n")
        else:
            c = c + 1

    fileinput.close()
    # conn.close()
    if float(c) / float(n) > 0.75:
        print("Warning: More than 75% of entries skipped. Are you using the "
              "correct database?", file=sys.stderr)


if __name__ == '__main__':
    pass
