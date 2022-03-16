import re
import sys
import fileinput
import pandas as pd
import logging

logger = logging.getLogger(__name__)

class Row:
    def __init__(self, row):
        l = row.rstrip().split("\t")
        if len(l) == 15:
            # Assume the format is generated from gtfToGenePred.
            # Add a dummy column in front of list to represent bin column in
            # UCSC genePred tables
            l.insert(0, None)
        if len(l) <= 15:
            raise ValueError("Insufficient number of columns in gene"
                             " prediction file: %s. If gtfToGenePred was used,"
                             " please ensure that the -genePredExt option is"
                             " enabled!" % row)
        try:
            self.name = l[1]
            self.chrom = l[2]
            self.strand = l[3]
            self.txStart = int(l[4])
            self.txEnd = int(l[5])
            self.cdsStart = int(l[6])
            self.cdsEnd = int(l[7])
            self.exonStarts = [int(x) for x in [_f for _f in l[9].split(",")
                            if _f]]
            self.exonEnds = [int(x) for x in [_f for _f in l[10].split(",")
                                    if _f]]
            self.name2 = l[12]
        except Exception as e:
            raise ValueError("Unable to parse the input genePred file."
                         " QAPA expects a genePred file from UCSC Table Browser"
                         " or converted using gtfToGenePred with -genePredExt"
                         " option enabled. \nExample row: %s" % row)

        self.utr3 = [self.cdsEnd, self.txEnd]
        self.has_intron_in_3utr = self.cdsEnd < self.exonStarts[-1]

        if self.strand == "-":
            self.utr3 = [self.txStart, self.cdsStart]
            self.has_intron_in_3utr = self.cdsStart > self.exonEnds[0]

    def extract_last_exon(self, n=1, min_utr_length=0):
        bed = None
        name = self._join_names()

        if not self.has_intron_in_3utr and \
                self.get_3utr_length() >= min_utr_length:
            if self.strand == "+":
                start = max([self.exonStarts[-n], self.cdsStart])
                bed = [self.chrom, start, self.txEnd, name,
                       self.get_3utr_length(), self.strand,
                       start, self.cdsEnd]
            else:
                end = min([self.exonEnds[n - 1], self.cdsEnd])
                bed = [self.chrom, self.txStart, end, name,
                       self.get_3utr_length(), self.strand, self.cdsStart,
                       end]
            bed.append(self.name2)
            bed.append(','.join([str(x) for x in self.exonStarts]))
            bed.append(','.join([str(x) for x in self.exonEnds]))
        #else:
            #logger.debug("Skipping %s because it contains an intron in 3' UTR" %
                    #self.name)
        return bed

    def extract_3utr(self, min_utr_length=0):
        bed = None
        name = self._join_names()

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
        return not re.match(r'^(chr)*[0-9XYM]+$', self.chrom)

    def chromosome_contains_underscore(self):
        return re.search(r'_', self.chrom)

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

    def _join_names(self):
        return get_stripped_name(self.name) + "_" + get_stripped_name(self.name2)


def get_stripped_name(name):
    # If Gencode tables are supplied, the Ensembl transcript ID has a
    # version number appended to the ID. We want to strip this out.
    match = re.match(r'ENS\w*T\d+', name)
    if match:
        return match.group()
    return name

def main(args, fout=sys.stdout):
    if args.debug:
        logger.setLevel(logging.DEBUG)

    conn = pd.read_table(args.db)
    conn = conn.loc[:, ['Transcript stable ID', 'Gene type',
                        'Transcript type']].drop_duplicates()
    conn = conn.set_index(['Transcript stable ID'])

    max_warnings = 10
    w = 0
    c = 0
    n = 0
    bad_chroms = set()
    for row in fileinput.input(args.annotation_file[0],
                               openhook=fileinput.hook_compressed):
        n = n + 1

        if fileinput.isfirstline():
            if row.startswith("#bin") or row.startswith("bin"):
                logger.debug("Header detected in genePred file. Assuming UCSC"
                             " format.")
                continue
            else:
                logger.debug("No header detected. Assuming custom genePred.")


        if row.startswith("#"):
            logger.debug("Skipping line %s", row)
            continue

        rowobj = Row(row)

        if not args.no_skip_random_chromosomes and \
            rowobj.is_on_random_chromosome():
            logger.debug("Skipping line %s (random chromosome)", row)
            c = c + 1
            continue

        if rowobj.chromosome_contains_underscore():
            w = w + 1

            if rowobj.chrom not in bad_chroms:
                logger.warning("Skipping chromosome %s because it contains"
                               " underscores" % rowobj.chrom)
                bad_chroms.add(rowobj.chrom)
            continue

        # filter for only protein-coding genes
        try:
            result = conn.loc[get_stripped_name(rowobj.name)]
            if isinstance(result, pd.DataFrame):
                result = result.iloc[0, ]
            if not (result['Gene type'] == "protein_coding" and
                    result['Transcript type'] == "protein_coding"):
                c = c + 1
                continue
        except KeyError:
            logger.debug("Skipping %s (not a protein coding gene)", rowobj.name)
            c = c + 1
            continue

        bed = rowobj.extract_last_exon()

        if bed is not None:
            fout.write("\t".join([str(x) for x in bed]) + "\n")
        else:
            logger.debug("Skipping %s (could not extract last exon",
                    rowobj.name)
            c = c + 1

    fileinput.close()

    if (skipped := c / n) > 0.99:
        raise RuntimeError("Too many entries, %d/%d (%0.2f%%), were skipped. "
              "Are you using the correct database?" % (c, n, 100 * skipped))
    elif skipped  > 0.75:
        logger.warning("%d/%d (%0.2f%%) were skipped. Are you using the "
              "correct database?" % (c, n, 100 * skipped))

